/**
 * @file decoder.c
 * @brief Duration grammar construction for CRF training and decoding.
 */
#include "include/crf.h"
#include "include/cwt.h"
#include "include/surf.h"
#include "include/utils.h"

#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

/**
 * @brief Add an incoming edge to a grammar state.
 * @param[in,out] grammar Grammar to modify.
 * @param[in] from Source state index.
 * @param[in] to Destination state index.
 * @return None.
 */
static void grammar_edge(Grammar *grammar, int from, int to) {
  int edge = grammar->backpointers[to]++;
  if (edge >= GRAMMAR_MAX_BACKPOINTERS)
    abort();
  grammar->from[to][edge] = from;
}

/**
 * @brief Add a minimum-duration state chain and its self-loop.
 * @param[in,out] grammar Grammar to modify.
 * @param[in] class Class emitted by each state in the chain.
 * @param[in] length Number of states in the chain.
 * @return Index of the first state in the chain.
 */
static int grammar_chain(Grammar *grammar, uint8_t class, int length) {
  int head = grammar->count;
  for (int i = 0; i < length; ++i) {
    grammar->cls[grammar->count] = class;
    if (i)
      grammar_edge(grammar, grammar->count - 1, grammar->count);
    grammar->count++;
  }
  grammar_edge(grammar, grammar->count - 1, grammar->count - 1);
  return head;
}

/**
 * @brief Build the duration-expanded gene grammar.
 *
 * The intergenic state connects to start and end exons. Exons connect through
 * strand-specific introns, and internal exons may repeat before reaching an
 * end exon. Minus-strand paths use the reverse transcript order.
 * @param[out] grammar Grammar to initialize.
 * @return None.
 */
void build_grammar(Grammar *grammar) {
  memset(grammar, 0, sizeof(*grammar));
  int intergenic = grammar_chain(grammar, INTERGENIC, 1);
  for (int strand = 0; strand < STRANDS; ++strand) {
    uint8_t start_class =
        strand ? MINUS_START_SINGLE_EXON : PLUS_START_SINGLE_EXON;
    uint8_t end_class = strand ? MINUS_END_EXON : PLUS_END_EXON;
    uint8_t internal_class = strand ? MINUS_INTERNAL_EXON : PLUS_INTERNAL_EXON;
    uint8_t intron_class = strand ? MINUS_INTRON : PLUS_INTRON;
    int start = grammar_chain(grammar, start_class, MIN_EXON_LENGTH);
    int end = grammar_chain(grammar, end_class, MIN_EXON_LENGTH);
    int internal = grammar_chain(grammar, internal_class, MIN_EXON_LENGTH);
    int intron = grammar_chain(grammar, intron_class, MIN_INTRON_LENGTH);
    int start_tail = start + MIN_EXON_LENGTH - 1;
    int end_tail = end + MIN_EXON_LENGTH - 1;
    int internal_tail = internal + MIN_EXON_LENGTH - 1;
    int intron_tail = intron + MIN_INTRON_LENGTH - 1;

    grammar_edge(grammar, intergenic, start);
    grammar_edge(grammar, start_tail, intergenic);
    if (strand)
      grammar_edge(grammar, intergenic, end);
    grammar_edge(grammar, strand ? end_tail : start_tail, intron);
    grammar_edge(grammar, intron_tail, internal);
    grammar_edge(grammar, internal_tail, intron);
    grammar_edge(grammar, intron_tail, strand ? start : end);
    grammar_edge(grammar, end_tail, intergenic);
  }
}

struct Decoder {
  Grammar grammar;
  int64_t length;
  int64_t position;
  double *previous;
  double *current;
  uint8_t *back;
};

Decoder *decoder_create(int64_t length) {
  if (length < 1)
    return NULL;
  Decoder *decoder = calloc(1, sizeof(*decoder));
  if (!decoder)
    abort();
  build_grammar(&decoder->grammar);
  decoder->length = length;
  decoder->previous = calloc((size_t)decoder->grammar.count, sizeof(double));
  decoder->current = calloc((size_t)decoder->grammar.count, sizeof(double));
  decoder->back = calloc((size_t)(length > 1 ? length - 1 : 1) *
                             (size_t)decoder->grammar.count,
                         sizeof(*decoder->back));
  if (!decoder->previous || !decoder->current || !decoder->back)
    abort();
  return decoder;
}

void decoder_destroy(Decoder *decoder) {
  if (!decoder)
    return;
  free(decoder->previous);
  free(decoder->current);
  free(decoder->back);
  free(decoder);
}

void decoder_step(Decoder *decoder, const double *emission,
                  const double *transition) {
  const Grammar *grammar = &decoder->grammar;
  int64_t position = decoder->position++;
  if (!position) {
    for (int state = 0; state < grammar->count; ++state)
      decoder->previous[state] = -INFINITY;
    decoder->previous[0] = emission[INTERGENIC];
    return;
  }
  double top = -INFINITY;
  for (int state = 0; state < grammar->count; ++state) {
    double best = -INFINITY;
    int choice = 0;
    uint8_t to_class = grammar->cls[state];

    for (int edge = 0; edge < grammar->backpointers[state]; ++edge) {
      int previous = grammar->from[state][edge];
      uint8_t from_class = grammar->cls[previous];
      double score = decoder->previous[previous];
      if (transition && from_class != to_class)
        score += transition[(size_t)from_class * CLASSES + to_class];
      if (score > best) {
        best = score;
        choice = edge;
      }
    }
    decoder->current[state] = best + emission[to_class];
    if (decoder->current[state] > top)
      top = decoder->current[state];
    decoder->back[(size_t)(position - 1) * (size_t)grammar->count +
                  (size_t)state] = (uint8_t)choice;
  }
  for (int state = 0; state < grammar->count; ++state)
    decoder->current[state] -= top;
  double *swap = decoder->previous;
  decoder->previous = decoder->current;
  decoder->current = swap;
}

int decoder_finish(Decoder *decoder, uint8_t *path) {
  const Grammar *grammar = &decoder->grammar;
  if (decoder->position != decoder->length || !isfinite(decoder->previous[0]))
    return ERANGE;
  int state = 0;
  for (int64_t position = decoder->length - 1; position > 0; --position) {
    path[position] = grammar->cls[state];
    int choice = decoder->back[(size_t)(position - 1) * (size_t)grammar->count +
                               (size_t)state];
    if (choice >= grammar->backpointers[state])
      return ERANGE;
    state = grammar->from[state][choice];
  }
  path[0] = grammar->cls[state];
  return state == 0 ? 0 : ERANGE;
}

typedef struct {
  char *parent;
  int sequence;
  char strand;
  IntervalList exons;
  int64_t orf_length;
} LabelTranscript;
static void label_fail(const char *path, long line, const char *message) {
  bye("%s:%ld: %s", path, line, message);
}
static LabelTranscript *label_get_transcript(LabelTranscript **items,
                                             int *count, int *capacity,
                                             const char *parent, int sequence,
                                             char strand) {
  for (int i = 0; i < *count; i++)
    if ((*items)[i].sequence == sequence && (*items)[i].strand == strand &&
        !strcmp((*items)[i].parent, parent))
      return &(*items)[i];
  if (*count == *capacity) {
    *capacity = *capacity ? *capacity * 2 : 64;
    *items = srealloc(*items, (size_t)*capacity, sizeof(**items));
  }
  LabelTranscript *t = &(*items)[(*count)++];
  *t = (LabelTranscript){
      .parent = sstrdup(parent), .sequence = sequence, .strand = strand};
  return t;
}
static int label_transcript_compare(const void *a, const void *b) {
  const LabelTranscript *x = a, *y = b;
  return x->orf_length < y->orf_length   ? 1
         : x->orf_length > y->orf_length ? -1
                                         : 0;
}
static void label_set_orf_length(LabelTranscript *t, const Sequence *s) {
  int64_t length;
  char *transcript = splice_sequence(s, &t->exons, t->strand, &length);
  int64_t start = INT64_MAX, end = 0;
  t->orf_length =
      longest_orf(transcript, length, &start, &end) ? end - start : 0;
  free(transcript);
}
static char *label_parent_id(char *attributes) {
  for (char *field = strtok(attributes, ";"); field;
       field = strtok(NULL, ";")) {
    while (*field == ' ')
      ++field;
    if (!strncmp(field, "Parent=", 7) && field[7])
      return field + 7;
  }
  return NULL;
}
void read_gff3(const char *path, Sequence *sequences, int count) {
  gzFile input = gzopen(path, "rb");
  if (!input)
    bye("cannot open annotation %s", path);
  LabelTranscript *transcripts = NULL;
  int transcript_count = 0, transcript_capacity = 0;
  char line[65536];
  long line_number = 0;
  while (gzgets(input, line, sizeof(line))) {
    ++line_number;
    if (line[0] == '#' || line[0] == '\n' || line[0] == '\r')
      continue;
    line[strcspn(line, "\r\n")] = '\0';
    char *fields[9], *cursor = line;
    int columns = 0;
    while (columns < 9 && cursor) {
      fields[columns++] = cursor;
      cursor = strchr(cursor, '\t');
      if (cursor)
        *cursor++ = '\0';
    }
    if (columns != 9 || cursor)
      label_fail(path, line_number, "expected 9 columns");
    if (strcmp(fields[2], "exon"))
      continue;
    int sequence = -1;
    for (int i = 0; i < count; ++i)
      if (!strcmp(sequences[i].name, fields[0]))
        sequence = i;
    if (sequence < 0)
      label_fail(path, line_number, "sequence absent from FASTA");
    char *end = NULL;
    errno = 0;
    int64_t start = strtoll(fields[3], &end, 10);
    if (errno || end == fields[3] || *end || start < 1)
      label_fail(path, line_number, "invalid start coordinate");
    errno = 0;
    int64_t finish = strtoll(fields[4], &end, 10);
    if (errno || end == fields[4] || *end || finish < start ||
        finish > sequences[sequence].length)
      label_fail(path, line_number, "invalid end coordinate");
    char strand = fields[6][0];
    if ((strand != '+' && strand != '-') || fields[6][1])
      label_fail(path, line_number, "unsupported strand");
    char *parent = label_parent_id(fields[8]);
    if (!parent)
      label_fail(path, line_number, "exon without Parent");
    LabelTranscript *t =
        label_get_transcript(&transcripts, &transcript_count,
                             &transcript_capacity, parent, sequence, strand);
    push_interval(&t->exons, start - 1, finish);
  }
  if (gzclose(input) != Z_OK)
    bye("cannot read annotation %s", path);
  for (int i = 0; i < transcript_count; ++i) {
    LabelTranscript *t = &transcripts[i];
    qsort(t->exons.items, (size_t)t->exons.count, sizeof(*t->exons.items),
          compare_interval);
    label_set_orf_length(t, &sequences[t->sequence]);
  }
  qsort(transcripts, (size_t)transcript_count, sizeof(*transcripts),
        label_transcript_compare);
  for (int i = 0; i < transcript_count; ++i) {
    LabelTranscript *t = &transcripts[i];
    Sequence *s = &sequences[t->sequence];
    if (!s->labels)
      s->labels = scalloc((size_t)s->length, sizeof(*s->labels));
    int conflict = 0;
    for (int e = 0; e < t->exons.count && !conflict; ++e)
      for (int64_t p = t->exons.items[e].start; p < t->exons.items[e].end; ++p)
        if (s->labels[p]) {
          conflict = 1;
          break;
        }
    if (conflict) {
      free(t->parent);
      free(t->exons.items);
      continue;
    }
    int plus = t->strand == '+';
    int first = plus ? 0 : t->exons.count - 1;
    int last = plus ? t->exons.count - 1 : 0;
    for (int e = 0; e < t->exons.count; ++e) {
      Interval exon = t->exons.items[e];
      uint8_t class =
          t->exons.count == 1 ? (plus ? PLUS_START_SINGLE_EXON : MINUS_END_EXON)
          : e == first
              ? (plus ? PLUS_START_SINGLE_EXON : MINUS_START_SINGLE_EXON)
          : e == last ? (plus ? PLUS_END_EXON : MINUS_END_EXON)
                      : (plus ? PLUS_INTERNAL_EXON : MINUS_INTERNAL_EXON);
      for (int64_t p = exon.start; p < exon.end; ++p)
        s->labels[p] = class;
    }
    uint8_t intron = plus ? PLUS_INTRON : MINUS_INTRON;
    for (int e = 0; e + 1 < t->exons.count; ++e)
      for (int64_t p = t->exons.items[e].end; p < t->exons.items[e + 1].start;
           ++p)
        s->labels[p] = intron;
    free(t->parent);
    free(t->exons.items);
  }
  free(transcripts);
}
static const char *const prediction_class_names[CLASSES] = {
    "intergenic",   "plus_start_single_exon",  "plus_end_exon",
    "plus_intron",  "minus_start_single_exon", "minus_end_exon",
    "minus_intron", "plus_internal_exon",      "minus_internal_exon"};
typedef struct {
  FILE *files[CLASSES];
  int value[CLASSES];
  int64_t start[CLASSES];
} PredictionBedgraph;
static void prediction_emit(PredictionBedgraph *b, const char *name, int class,
                            int64_t position) {
  if (b->value[class] >= 0 && position > b->start[class] &&
      fprintf(b->files[class], "%s\t%lld\t%lld\t%d.%02d\n", name,
              (long long)b->start[class], (long long)position,
              b->value[class] / 100, b->value[class] % 100) < 0)
    bye("cannot write bedGraph");
}
static void prediction_add(PredictionBedgraph *b, const char *name,
                           int64_t position, const double *scores) {
  if (!b->files[0])
    return;
  double top = -INFINITY, total = 0;
  for (int c = 0; c < CLASSES; c++)
    if (scores[c] > top)
      top = scores[c];
  for (int c = 0; c < CLASSES; c++)
    total += exp(scores[c] - top);
  for (int c = 0; c < CLASSES; c++) {
    int value = (int)lrint(exp(scores[c] - top) / total * 100.0);
    if (value == b->value[c])
      continue;
    prediction_emit(b, name, c, position);
    b->value[c] = value;
    b->start[c] = position;
  }
}
unsigned long predict(const Model *model, const Sequence *sequences, int count,
                      FILE *gff, const char *prefix) {
  PredictionBedgraph b = {0};
  for (int c = 0; prefix && c < CLASSES; c++) {
    char path[PATH_BUFFER];
    if (snprintf(path, sizeof(path), "%s.%s.bedgraph", prefix,
                 prediction_class_names[c]) >= (int)sizeof(path))
      bye("output prefix too long");
    if (!(b.files[c] = fopen(path, "w")))
      bye("cannot create %s", path);
    fprintf(b.files[c], "track type=bedGraph name=\"%s\"\n",
            prediction_class_names[c]);
  }
  unsigned long genes = 0;
  for (int si = 0; si < count; si++) {
    const Sequence *s = &sequences[si];
    Decoder *d = decoder_create(s->length);
    /* Context scoring finishes before traceback reuses this buffer as labels.
     */
    uint8_t *path = scalloc((size_t)s->length, sizeof(*path));
    for (int64_t p = 0; p < s->length; ++p) {
      int base = base_index(s->bases[p]);
      path[p] = (uint8_t)(base < 0 ? NUCLEOTIDES : base);
    }
    Sample sample = {s->length, path, NULL, NULL};
    for (int c = 0; c < CLASSES; c++)
      b.value[c] = -1;
    float *power =
        scalloc(smul((size_t)s->length, CWT_CHANNELS), sizeof(*power));
    cwt_power(s->bases, s->length, 0, s->length, power);
    for (int64_t p = 0; p < s->length; p++) {
      float row[INPUT_CHANNELS];
      double scores[CLASSES], transition[TRANSITIONS];
      compute_features(s->bases[p], power + p * CWT_CHANNELS, row);
      compute_emissions(model, row, scores);
      if (!compute_transitions(model->transition_weights,
                               model->boundary_weights, &sample, (size_t)p,
                               transition))
        bye("non-finite transition score for %s", s->name);
      prediction_add(&b, s->name, p, scores);
      decoder_step(d, scores, transition);
    }
    free(power);
    for (int c = 0; c < CLASSES; c++)
      if (b.files[c])
        prediction_emit(&b, s->name, c, s->length);
    if (decoder_finish(d, path))
      bye("grammar decoding failed for %s", s->name);
    size_t coding = 0;
    for (int64_t p = 0; p < s->length; p++)
      if (path[p] != INTERGENIC)
        coding++;
    cry("%s: %zu non-intergenic bases", s->name, coding);
    decoder_destroy(d);
    long written = write_gff(gff, s, path, &genes);
    free(path);
    if (written < 0)
      bye("cannot write GFF3");
    cry("%s: %ld genes", s->name, written);
  }
  for (int c = 0; c < CLASSES; c++)
    if (b.files[c] && fclose(b.files[c]))
      bye("cannot close bedGraph");
  return genes;
}
