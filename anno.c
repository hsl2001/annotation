/**
 * @file anno.c
 * @brief anno: FASTA/GFF3 input, CWT, pure-C 1D UNet (training + inference),
 *        9-state grammar Viterbi, longest ORF and GFF3 output.
 *
 * Pipeline: bases -> one-hot + CWT power -> UNet -> per-base probability of
 * ANNO_CLASSES classes -> log emissions -> grammar Viterbi -> transcripts ->
 * longest ORF -> gene/mRNA/exon/UTR/CDS records.
 */
#define _POSIX_C_SOURCE 200809L

#include "anno.h"

#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>
#include <zlib.h>

/** pi, since C11 does not define M_PI. */
#define ANNO_PI 3.14159265358979323846

/* ------------------------------------------------------------------------ */
/* Basic utilities                                                          */
/* ------------------------------------------------------------------------ */

/** Class names; also used for bedGraph file names. */
static const char *const class_names[ANNO_CLASSES] = {
  "intergenic", "plus_five_prime_UTR", "plus_exon", "plus_intron",
  "plus_three_prime_UTR", "minus_five_prime_UTR", "minus_exon",
  "minus_intron", "minus_three_prime_UTR"
};
/** Morlet scales of the CWT channels. */
static const int cwt_scales[ANNO_CWT_CHANNELS] = {4, 5, 6};

/** @brief Print an error to stderr and exit. */
static _Noreturn void fail(const char *format, ...) {
  va_list arguments;
  va_start(arguments, format);
  fputs("anno: error: ", stderr);
  vfprintf(stderr, format, arguments);
  fputc('\n', stderr);
  va_end(arguments);
  exit(EXIT_FAILURE);
}

/** @brief Print a log line to stderr. */
static void say(const char *format, ...) {
  va_list arguments;
  va_start(arguments, format);
  fputs("anno: ", stderr);
  vfprintf(stderr, format, arguments);
  fputc('\n', stderr);
  va_end(arguments);
}

/** @brief Multiply two sizes, failing on overflow. */
static size_t mul_size(size_t a, size_t b) {
  if (a && b > SIZE_MAX / a) fail("size overflow");
  return a * b;
}

/** @brief calloc() that fails on overflow or exhaustion. Caller frees. */
static void *xcalloc(size_t count, size_t size) {
  void *memory = calloc(count ? count : 1, size ? size : 1);
  if (!memory || (size && count > SIZE_MAX / size)) fail("out of memory");
  return memory;
}

/** @brief realloc() that fails on overflow or exhaustion. */
static void *xrealloc(void *old, size_t count, size_t size) {
  void *memory = realloc(old, mul_size(count ? count : 1, size));
  if (!memory) fail("out of memory");
  return memory;
}

/** @brief strdup() that fails on exhaustion. Caller frees. */
static char *xstrdup(const char *text) {
  char *copy = strdup(text);
  if (!copy) fail("out of memory");
  return copy;
}

/** @brief Monotonic wall-clock seconds. */
static double now_seconds(void) {
  struct timespec time;
  clock_gettime(CLOCK_MONOTONIC, &time);
  return (double)time.tv_sec + 1e-9 * (double)time.tv_nsec;
}

/** @brief A/C/G/T -> 0..3, anything else -> -1. */
static int base_index(char base) {
  switch (base) {
  case 'A': return 0;
  case 'C': return 1;
  case 'G': return 2;
  case 'T': return 3;
  default: return -1;
  }
}

/** @brief Watson-Crick complement (N for unknown bases). */
static char complement(char base) {
  int index = base_index(base);
  return index < 0 ? 'N' : "TGCA"[index];
}

/** @brief xorshift64* generator; deterministic for a given seed. */
static uint64_t random_next(uint64_t *state) {
  *state ^= *state >> 12;
  *state ^= *state << 25;
  *state ^= *state >> 27;
  return *state * 2685821657736338717ULL;
}

/** @brief Uniform float in [0, 1). */
static float random_unit(uint64_t *state) {
  return (float)(random_next(state) >> 40) / 16777216.0f;
}

/** @brief Uniform integer in [0, bound). bound must be positive. */
static int64_t random_below(uint64_t *state, int64_t bound) {
  return (int64_t)(random_next(state) % (uint64_t)bound);
}

/* ------------------------------------------------------------------------ */
/* Threads                                                                  */
/* ------------------------------------------------------------------------ */

/** Job callback: (context, job index, worker index). */
typedef void (*JobFunction)(void *context, int64_t job, int worker);

/** Arguments of one worker thread. */
typedef struct {
  JobFunction function;
  void *context;
  int64_t jobs;
  int worker;
  int workers;
} WorkerArgs;

/** @brief Worker body: runs jobs worker, worker+workers, ... */
static void *worker_main(void *pointer) {
  WorkerArgs *args = pointer;
  for (int64_t job = args->worker; job < args->jobs; job += args->workers)
    args->function(args->context, job, args->worker);
  return NULL;
}

/**
 * @brief Run jobs [0,jobs) on up to `threads` threads with a static,
 *        deterministic job-to-worker assignment.
 */
static void run_parallel(int threads, int64_t jobs, JobFunction function,
                         void *context) {
  int workers = threads < 1 ? 1 : threads;
  if (jobs < workers) workers = (int)(jobs > 0 ? jobs : 1);
  pthread_t *ids = xcalloc((size_t)workers, sizeof(*ids));
  WorkerArgs *args = xcalloc((size_t)workers, sizeof(*args));
  for (int i = 0; i < workers; i++) {
    args[i] = (WorkerArgs){function, context, jobs, i, workers};
    if (i && pthread_create(&ids[i], NULL, worker_main, &args[i]))
      fail("cannot create thread");
  }
  worker_main(&args[0]);
  for (int i = 1; i < workers; i++) pthread_join(ids[i], NULL);
  free(args);
  free(ids);
}

/* ------------------------------------------------------------------------ */
/* FASTA                                                                    */
/* ------------------------------------------------------------------------ */

AnnoSequence *anno_read_fasta(const char *path, int *count_out) {
  gzFile input = gzopen(path, "rb");
  if (!input) fail("cannot open FASTA %s", path);
  AnnoSequence *sequences = NULL;
  int count = 0;
  size_t capacity = 0, bases_capacity = 0;
  AnnoSequence *current = NULL;
  char buffer[65536];
  while (gzgets(input, buffer, sizeof(buffer))) {
    size_t length = strlen(buffer);
    if (buffer[0] == '>') {
      char *name = buffer + 1;
      name[strcspn(name, " \t\r\n")] = '\0';
      if (!*name) fail("empty FASTA sequence ID in %s", path);
      for (int i = 0; i < count; i++)
        if (!strcmp(sequences[i].name, name))
          fail("duplicate FASTA sequence ID %s", name);
      if ((size_t)count == capacity) {
        capacity = capacity ? capacity * 2 : 8;
        sequences = xrealloc(sequences, capacity, sizeof(*sequences));
      }
      current = &sequences[count++];
      *current = (AnnoSequence){.name = xstrdup(name)};
      bases_capacity = 0;
      continue;
    }
    for (size_t i = 0; i < length; i++) {
      if (isspace((unsigned char)buffer[i])) continue;
      if (!current) fail("FASTA sequence before header in %s", path);
      if ((size_t)current->length + 2 > bases_capacity) {
        bases_capacity = bases_capacity ? bases_capacity * 2 : 1 << 20;
        current->bases = xrealloc(current->bases, bases_capacity, 1);
      }
      current->bases[current->length++] = (char)toupper((unsigned char)buffer[i]);
    }
  }
  int eof = gzeof(input);
  if (gzclose(input) != Z_OK || !eof) fail("cannot read FASTA %s", path);
  if (!count) fail("no FASTA records in %s", path);
  for (int i = 0; i < count; i++) {
    if (!sequences[i].length) fail("empty FASTA sequence %s", sequences[i].name);
    sequences[i].bases[sequences[i].length] = '\0';
  }
  *count_out = count;
  return sequences;
}

void anno_sequences_free(AnnoSequence *sequences, int count) {
  for (int i = 0; sequences && i < count; i++) {
    free(sequences[i].name);
    free(sequences[i].bases);
    free(sequences[i].labels);
  }
  free(sequences);
}

/* ------------------------------------------------------------------------ */
/* GFF3 -> label vector                                                     */
/* ------------------------------------------------------------------------ */

enum { KIND_EXON, KIND_CDS, KIND_5UTR, KIND_3UTR, KIND_INTRON };

/** One GFF feature attached to one transcript. Owns `parent`. */
typedef struct {
  char *parent;
  int sequence;
  int kind;
  char strand;
  int64_t start, end; /* 0-based half-open */
} GffFeature;

/** Half-open genomic interval. */
typedef struct {
  int64_t start, end;
} Segment;

/** Growable segment list. */
typedef struct {
  Segment *items;
  int count, capacity;
} SegmentList;

/** @brief Append a segment. */
static void segments_push(SegmentList *list, int64_t start, int64_t end) {
  if (list->count == list->capacity) {
    list->capacity = list->capacity ? list->capacity * 2 : 8;
    list->items = xrealloc(list->items, (size_t)list->capacity,
                           sizeof(*list->items));
  }
  list->items[list->count++] = (Segment){start, end};
}

/** @brief Order segments by start. */
static int segment_compare(const void *a, const void *b) {
  const Segment *x = a, *y = b;
  return x->start < y->start ? -1 : x->start > y->start ? 1
       : x->end < y->end ? -1 : x->end > y->end;
}

/** @brief Order features by parent, kind and start. */
static int feature_compare(const void *a, const void *b) {
  const GffFeature *x = a, *y = b;
  int order = strcmp(x->parent, y->parent);
  if (order) return order;
  if (x->kind != y->kind) return x->kind - y->kind;
  return x->start < y->start ? -1 : x->start > y->start;
}

/** @brief GFF feature type -> KIND_*, or -1 when ignored. */
static int feature_kind(const char *type) {
  if (!strcmp(type, "exon")) return KIND_EXON;
  if (!strcmp(type, "CDS")) return KIND_CDS;
  if (!strcmp(type, "five_prime_UTR")) return KIND_5UTR;
  if (!strcmp(type, "three_prime_UTR")) return KIND_3UTR;
  if (!strcmp(type, "intron")) return KIND_INTRON;
  if (!strcmp(type, "gene") || !strcmp(type, "mRNA") ||
      !strcmp(type, "transcript"))
    return -2;
  return -1;
}

/** @brief Parse a non-negative decimal integer. @return 1 on success. */
static int parse_int64(const char *text, int64_t *value) {
  char *end = NULL;
  errno = 0;
  long long parsed = strtoll(text, &end, 10);
  if (errno || end == text || *end || parsed < 0) return 0;
  *value = parsed;
  return 1;
}

/** One coding transcript ready to be painted. Owns its lists. */
typedef struct {
  int sequence;
  char strand;
  int64_t start, end;
  SegmentList exons, explicit_utr5, explicit_utr3, explicit_introns;
  int64_t cds_start, cds_end;
} LabelTranscript;

/** @brief Order transcripts by decreasing span (longer transcripts win). */
static int transcript_compare(const void *a, const void *b) {
  const LabelTranscript *x = a, *y = b;
  int64_t lx = x->end - x->start, ly = y->end - y->start;
  if (lx != ly) return lx > ly ? -1 : 1;
  if (x->sequence != y->sequence) return x->sequence - y->sequence;
  return x->start < y->start ? -1 : x->start > y->start;
}

/** @brief Merge sorted touching/overlapping segments in place. */
static void segments_merge(SegmentList *list) {
  if (!list->count) return;
  qsort(list->items, (size_t)list->count, sizeof(Segment), segment_compare);
  int out = 0;
  for (int i = 1; i < list->count; i++) {
    if (list->items[i].start <= list->items[out].end) {
      if (list->items[i].end > list->items[out].end)
        list->items[out].end = list->items[i].end;
    } else {
      list->items[++out] = list->items[i];
    }
  }
  list->count = out + 1;
}

/** @brief Read GFF features that belong to transcripts. */
static GffFeature *read_gff_features(const char *path,
                                     const AnnoSequence *sequences, int count,
                                     size_t *feature_count) {
  gzFile input = gzopen(path, "rb");
  if (!input) fail("cannot open annotation %s", path);
  GffFeature *features = NULL;
  size_t used = 0, capacity = 0;
  char line[65536];
  long line_number = 0;
  while (gzgets(input, line, sizeof(line))) {
    line_number++;
    if (!strncmp(line, "##FASTA", 7)) break;
    if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;
    line[strcspn(line, "\r\n")] = '\0';
    char *fields[9];
    char *cursor = line;
    int columns = 0;
    while (columns < 9 && cursor) {
      fields[columns++] = cursor;
      cursor = strchr(cursor, '\t');
      if (cursor) *cursor++ = '\0';
    }
    if (columns != 9 || cursor) fail("%s:%ld: expected 9 columns", path, line_number);
    int kind = feature_kind(fields[2]);
    if (kind == -1) continue;
    int sequence = -1;
    for (int i = 0; i < count; i++)
      if (!strcmp(sequences[i].name, fields[0])) sequence = i;
    if (sequence < 0) fail("%s:%ld: sequence %s absent from FASTA", path,
                           line_number, fields[0]);
    int64_t start = 0, end = 0;
    if (!parse_int64(fields[3], &start) || !parse_int64(fields[4], &end) ||
        start < 1 || end < start)
      fail("%s:%ld: invalid coordinates", path, line_number);
    if (end > sequences[sequence].length)
      fail("%s:%ld: coordinate beyond %s length %lld", path, line_number,
           fields[0], (long long)sequences[sequence].length);
    char strand = fields[6][0];
    if ((strand != '+' && strand != '-') || fields[6][1])
      fail("%s:%ld: unsupported strand %s", path, line_number, fields[6]);
    if (kind < 0) continue;
    char *parents = NULL;
    for (char *field = strtok(fields[8], ";"); field; field = strtok(NULL, ";")) {
      while (*field == ' ') field++;
      if (!strncmp(field, "Parent=", 7)) parents = field + 7;
    }
    if (!parents || !*parents)
      fail("%s:%ld: %s without Parent", path, line_number, fields[2]);
    for (char *save = NULL, *parent = strtok_r(parents, ",", &save); parent;
         parent = strtok_r(NULL, ",", &save)) {
      if (used == capacity) {
        capacity = capacity ? capacity * 2 : 1024;
        features = xrealloc(features, capacity, sizeof(*features));
      }
      features[used++] = (GffFeature){xstrdup(parent), sequence, kind, strand,
                                      start - 1, end};
    }
  }
  if (gzclose(input) != Z_OK) fail("cannot read annotation %s", path);
  *feature_count = used;
  return features;
}

/** @brief Paint one transcript into `labels` (caller checked overlap). */
static void paint_transcript(const LabelTranscript *t, uint8_t *labels) {
  int plus = t->strand == '+';
  uint8_t intron = plus ? ANNO_PLUS_INTRON : ANNO_MINUS_INTRON;
  uint8_t coding = plus ? ANNO_PLUS_EXON : ANNO_MINUS_EXON;
  uint8_t utr5 = plus ? ANNO_PLUS_5UTR : ANNO_MINUS_5UTR;
  uint8_t utr3 = plus ? ANNO_PLUS_3UTR : ANNO_MINUS_3UTR;
  for (int64_t i = t->start; i < t->end; i++) labels[i] = intron;
  for (int e = 0; e < t->exons.count; e++)
    for (int64_t i = t->exons.items[e].start; i < t->exons.items[e].end; i++)
      labels[i] = i < t->cds_start ? (plus ? utr5 : utr3)
                : i >= t->cds_end ? (plus ? utr3 : utr5) : coding;
  /* explicit UTR/intron features outrank the inferred labels */
  for (int e = 0; e < t->explicit_utr5.count; e++)
    for (int64_t i = t->explicit_utr5.items[e].start;
         i < t->explicit_utr5.items[e].end; i++)
      labels[i] = utr5;
  for (int e = 0; e < t->explicit_utr3.count; e++)
    for (int64_t i = t->explicit_utr3.items[e].start;
         i < t->explicit_utr3.items[e].end; i++)
      labels[i] = utr3;
  for (int e = 0; e < t->explicit_introns.count; e++)
    for (int64_t i = t->explicit_introns.items[e].start;
         i < t->explicit_introns.items[e].end; i++)
      labels[i] = intron;
}

void anno_read_labels(const char *path, AnnoSequence *sequences, int count) {
  size_t feature_count = 0;
  GffFeature *features = read_gff_features(path, sequences, count, &feature_count);
  qsort(features, feature_count, sizeof(*features), feature_compare);
  LabelTranscript *transcripts = NULL;
  size_t transcript_count = 0, skipped = 0;
  for (size_t first = 0; first < feature_count;) {
    size_t last = first;
    while (last < feature_count &&
           !strcmp(features[last].parent, features[first].parent))
      last++;
    LabelTranscript t = {.sequence = features[first].sequence,
                         .strand = features[first].strand,
                         .cds_start = INT64_MAX, .cds_end = -1};
    SegmentList pieces = {0};
    for (size_t i = first; i < last; i++) {
      const GffFeature *f = &features[i];
      if (f->sequence != t.sequence || f->strand != t.strand)
        fail("transcript %s mixes sequences or strands", f->parent);
      if (f->kind == KIND_EXON) segments_push(&t.exons, f->start, f->end);
      if (f->kind == KIND_5UTR) segments_push(&t.explicit_utr5, f->start, f->end);
      if (f->kind == KIND_3UTR) segments_push(&t.explicit_utr3, f->start, f->end);
      if (f->kind == KIND_INTRON)
        segments_push(&t.explicit_introns, f->start, f->end);
      if (f->kind == KIND_CDS) {
        if (f->start < t.cds_start) t.cds_start = f->start;
        if (f->end > t.cds_end) t.cds_end = f->end;
      }
      if (f->kind == KIND_CDS || f->kind == KIND_5UTR || f->kind == KIND_3UTR)
        segments_push(&pieces, f->start, f->end);
    }
    if (!t.exons.count) {
      t.exons = pieces;
      pieces = (SegmentList){0};
    }
    free(pieces.items);
    segments_merge(&t.exons);
    if (t.cds_end < 0 || !t.exons.count) { /* non-coding transcript */
      free(t.exons.items);
      free(t.explicit_utr5.items);
      free(t.explicit_utr3.items);
      free(t.explicit_introns.items);
      first = last;
      continue;
    }
    t.start = t.exons.items[0].start;
    t.end = t.exons.items[t.exons.count - 1].end;
    transcripts = xrealloc(transcripts, transcript_count + 1, sizeof(*transcripts));
    transcripts[transcript_count++] = t;
    first = last;
  }
  for (size_t i = 0; i < feature_count; i++) free(features[i].parent);
  free(features);
  qsort(transcripts, transcript_count, sizeof(*transcripts), transcript_compare);
  for (int i = 0; i < count; i++) {
    free(sequences[i].labels);
    sequences[i].labels = xcalloc((size_t)sequences[i].length, 1);
  }
  for (size_t k = 0; k < transcript_count; k++) {
    LabelTranscript *t = &transcripts[k];
    uint8_t *labels = sequences[t->sequence].labels;
    int free_span = 1;
    for (int64_t i = t->start; i < t->end && free_span; i++)
      free_span = !labels[i];
    if (free_span) paint_transcript(t, labels);
    else skipped++;
    free(t->exons.items);
    free(t->explicit_utr5.items);
    free(t->explicit_utr3.items);
    free(t->explicit_introns.items);
  }
  say("labelled %zu coding transcripts (%zu overlapping skipped)",
      transcript_count - skipped, skipped);
  free(transcripts);
}

/* ------------------------------------------------------------------------ */
/* CWT                                                                      */
/* ------------------------------------------------------------------------ */

void anno_cwt(const char *bases, int64_t length, int64_t start, int64_t end,
              float *out) {
  static const double signal[4][2] = {{1, 1}, {-1, 1}, {-1, -1}, {1, -1}};
  for (int s = 0; s < ANNO_CWT_CHANNELS; s++) {
    int radius = 4 * cwt_scales[s];
    double g[64], c[64], sn[64];
    for (int o = -radius; o <= radius; o++) {
      double x = (double)o / cwt_scales[s];
      g[o + radius] = exp(-0.5 * x * x);
      c[o + radius] = g[o + radius] * cos(6.0 * x);
      sn[o + radius] = g[o + radius] * sin(6.0 * x);
    }
    for (int64_t p = start; p < end; p++) {
      double real = 0, imaginary = 0, norm = 0;
      for (int o = -radius; o <= radius; o++) {
        int64_t q = p + o;
        int b = q >= 0 && q < length ? base_index(bases[q]) : -1;
        if (b < 0) continue;
        real += signal[b][0] * c[o + radius] - signal[b][1] * sn[o + radius];
        imaginary += signal[b][0] * sn[o + radius] + signal[b][1] * c[o + radius];
        norm += g[o + radius];
      }
      float *value = out + ((size_t)(p - start) * ANNO_CWT_CHANNELS + (size_t)s) * 2;
      value[0] = norm > 0 ? (float)(real / norm) : 0.0f;
      value[1] = norm > 0 ? (float)(imaginary / norm) : 0.0f;
    }
  }
}

/** @brief CWT power (|W|^2) for [start,end) into out[(p-start)*C + s]. */
static void cwt_power(const char *bases, int64_t length, int64_t start,
                      int64_t end, float *out) {
  for (int64_t block = start; block < end; block += ANNO_CWT_WINDOW) {
    float values[ANNO_CWT_WINDOW * ANNO_CWT_CHANNELS * 2];
    int64_t stop = block + ANNO_CWT_WINDOW < end ? block + ANNO_CWT_WINDOW : end;
    anno_cwt(bases, length, block, stop, values);
    for (int64_t i = 0; i < (stop - block) * ANNO_CWT_CHANNELS; i++)
      out[(size_t)(block - start) * ANNO_CWT_CHANNELS + (size_t)i] =
        values[2 * i] * values[2 * i] + values[2 * i + 1] * values[2 * i + 1];
  }
}

/** Context of the parallel whole-sequence CWT power job. */
typedef struct {
  const AnnoSequence *sequence;
  float *power;
} PowerJob;

/** @brief Job: CWT power of one ANNO_CWT_WINDOW block. */
static void power_job(void *pointer, int64_t job, int worker) {
  (void)worker;
  PowerJob *context = pointer;
  int64_t start = job * ANNO_CWT_WINDOW, end = start + ANNO_CWT_WINDOW;
  if (end > context->sequence->length) end = context->sequence->length;
  cwt_power(context->sequence->bases, context->sequence->length, start, end,
            context->power + (size_t)start * ANNO_CWT_CHANNELS);
}

/** @brief Whole-sequence CWT power computed in parallel blocks. Caller frees. */
static float *sequence_power(const AnnoSequence *sequence, int threads) {
  PowerJob job = {sequence, xcalloc(mul_size((size_t)sequence->length,
                                             ANNO_CWT_CHANNELS), sizeof(float))};
  run_parallel(threads, (sequence->length + ANNO_CWT_WINDOW - 1) / ANNO_CWT_WINDOW,
               power_job, &job);
  return job.power;
}

/** Context of the cwt-subcommand block job. */
typedef struct {
  const AnnoSequence *sequence;
  int64_t start;
  float *out;
} ComplexJob;

/** @brief Job: complex CWT of one block for the cwt subcommand. */
static void complex_job(void *pointer, int64_t job, int worker) {
  (void)worker;
  ComplexJob *context = pointer;
  int64_t start = context->start + job * ANNO_CWT_WINDOW;
  int64_t end = start + ANNO_CWT_WINDOW;
  if (end > context->sequence->length) end = context->sequence->length;
  if (start < end)
    anno_cwt(context->sequence->bases, context->sequence->length, start, end,
             context->out + (size_t)(start - context->start) * ANNO_CWT_CHANNELS * 2);
}

/**
 * @brief `anno cwt`: write {dir}/matrix.bin (complex64, plus rows then the
 *        reverse-complement rows) and {dir}/contigs.tsv for plot.py.
 *        Minus rows are the conjugated, mirrored plus rows (no recompute).
 */
static void run_cwt(const char *fasta, const char *directory, int threads) {
  int count = 0;
  AnnoSequence *sequences = anno_read_fasta(fasta, &count);
  char command_path[4096];
  if (snprintf(command_path, sizeof(command_path), "%s/matrix.bin", directory) >=
      (int)sizeof(command_path))
    fail("CWT directory path too long");
  if (mkdir(directory, 0777) && errno != EEXIST) fail("cannot create %s", directory);
  FILE *matrix = fopen(command_path, "w+b");
  snprintf(command_path, sizeof(command_path), "%s/contigs.tsv", directory);
  FILE *index = fopen(command_path, "w");
  if (!matrix || !index) fail("cannot create CWT output in %s", directory);
  fprintf(index, "# anno-cwt-v2\n# scales\t%d,%d,%d\nseqid\tlength\tplus_offset\t"
          "minus_offset\n", cwt_scales[0], cwt_scales[1], cwt_scales[2]);
  const int64_t chunk = (int64_t)ANNO_CWT_WINDOW * 256;
  const size_t row = ANNO_CWT_CHANNELS * 2;
  float *buffer = xcalloc(mul_size((size_t)chunk, row), sizeof(float));
  int64_t rows = 0;
  for (int i = 0; i < count; i++) {
    const AnnoSequence *s = &sequences[i];
    fprintf(index, "%s\t%lld\t%lld\t%lld\n", s->name, (long long)s->length,
            (long long)rows, (long long)(rows + s->length));
    long base = (long)(rows * (int64_t)row * (int64_t)sizeof(float));
    for (int64_t start = 0; start < s->length; start += chunk) {
      int64_t n = s->length - start < chunk ? s->length - start : chunk;
      ComplexJob job = {s, start, buffer};
      run_parallel(threads, (n + ANNO_CWT_WINDOW - 1) / ANNO_CWT_WINDOW,
                   complex_job, &job);
      if (fwrite(buffer, sizeof(float) * row, (size_t)n, matrix) != (size_t)n)
        fail("cannot write CWT matrix");
    }
    for (int64_t end = s->length; end > 0; end -= chunk) {
      int64_t n = end < chunk ? end : chunk, start = end - n;
      if (fseek(matrix, base + (long)(start * (int64_t)row * (int64_t)sizeof(float)),
                SEEK_SET) ||
          fread(buffer, sizeof(float) * row, (size_t)n, matrix) != (size_t)n ||
          fseek(matrix, 0, SEEK_END))
        fail("cannot reread CWT matrix");
      for (int64_t a = 0, b = n - 1; a <= b; a++, b--)
        for (size_t k = 0; k < row; k++) {
          float x = buffer[(size_t)a * row + k], y = buffer[(size_t)b * row + k];
          buffer[(size_t)a * row + k] = (k & 1) ? -y : y;
          buffer[(size_t)b * row + k] = (k & 1) ? -x : x;
        }
      if (fwrite(buffer, sizeof(float) * row, (size_t)n, matrix) != (size_t)n)
        fail("cannot write CWT matrix");
    }
    rows += 2 * s->length;
  }
  free(buffer);
  if (fclose(matrix) || fclose(index)) fail("cannot close CWT output");
  anno_sequences_free(sequences, count);
}

/* ------------------------------------------------------------------------ */
/* 1D UNet                                                                  */
/* ------------------------------------------------------------------------ */

/** Convolutions of one level: in->C, C->C, (C+C')->C, C->C. */
enum { CONV_IN, CONV_MID, CONV_UP, CONV_OUT, CONVS_PER_LEVEL };

struct Unet {
  UnetConfig config;
  float cwt_mean[ANNO_CWT_CHANNELS];
  float cwt_std[ANNO_CWT_CHANNELS];
  size_t parameters;
  size_t offset[ANNO_MAX_DEPTH][CONVS_PER_LEVEL];
  size_t head;
  float *weights; /* owned */
};

/**
 * Per-thread buffers for a window of up to `capacity` positions, laid out
 * position-major ([position][channel]). Owned by one worker only.
 */
typedef struct {
  int64_t capacity;
  float *input, *logits, *probability, *glogits;
  float *a[ANNO_MAX_DEPTH], *b[ANNO_MAX_DEPTH], *pooled[ANNO_MAX_DEPTH];
  float *cat[ANNO_MAX_DEPTH], *c[ANNO_MAX_DEPTH], *d[ANNO_MAX_DEPTH];
  float *ga[ANNO_MAX_DEPTH], *gb[ANNO_MAX_DEPTH], *gin[ANNO_MAX_DEPTH];
  float *gcat[ANNO_MAX_DEPTH], *gc[ANNO_MAX_DEPTH], *gout[ANNO_MAX_DEPTH];
  float *gradient; /* parameter gradient */
  int threads;      /* threads used inside each convolution */
  float *partial;   /* threads x largest convolution: weight-gradient partials */
  uint8_t *labels;
} Workspace;

UnetConfig unet_default_config(void) {
  return (UnetConfig){.depth = 4, .kernel = 9, .pool = 4,
                      .channels = {16, 32, 48, 64}};
}

/** @brief Channels entering level l. */
static int level_input(const UnetConfig *config, int level) {
  return level ? config->channels[level - 1] : ANNO_INPUT_CHANNELS;
}

/** @brief Number of weights + biases of one convolution. */
static size_t conv_size(int in, int out, int kernel) {
  return (size_t)out * ((size_t)in * (size_t)kernel + 1);
}

/** @brief Pooling multiple every window length must have. */
static int64_t unet_multiple(const UnetConfig *config) {
  int64_t multiple = 1;
  for (int l = 1; l < config->depth; l++) multiple *= config->pool;
  return multiple;
}

/** @brief Validate a configuration. */
static int config_valid(const UnetConfig *config) {
  if (!config || config->depth < 1 || config->depth > ANNO_MAX_DEPTH ||
      config->kernel < 1 || config->kernel > 31 || !(config->kernel & 1) ||
      config->pool < 2 || config->pool > 8)
    return 0;
  for (int l = 0; l < config->depth; l++)
    if (config->channels[l] < 1 || config->channels[l] > 256) return 0;
  return unet_multiple(config) <= ANNO_UNET_WINDOW;
}

/** @brief Compute parameter offsets for every convolution. */
static void unet_layout(Unet *unet) {
  const UnetConfig *config = &unet->config;
  size_t offset = 0;
  for (int l = 0; l < config->depth; l++) {
    int in = level_input(config, l), c = config->channels[l];
    unet->offset[l][CONV_IN] = offset;
    offset += conv_size(in, c, config->kernel);
    unet->offset[l][CONV_MID] = offset;
    offset += conv_size(c, c, config->kernel);
    if (l + 1 < config->depth) {
      unet->offset[l][CONV_UP] = offset;
      offset += conv_size(c + config->channels[l + 1], c, config->kernel);
      unet->offset[l][CONV_OUT] = offset;
      offset += conv_size(c, c, config->kernel);
    }
  }
  unet->head = offset;
  unet->parameters = offset + conv_size(config->channels[0], ANNO_CLASSES, 1);
}

/** @brief He-uniform initialisation of one convolution, zero bias. */
static void init_conv(float *weights, int in, int out, int kernel,
                      uint64_t *seed) {
  float bound = sqrtf(6.0f / (float)(in * kernel));
  for (size_t i = 0; i < (size_t)in * (size_t)kernel * (size_t)out; i++)
    weights[i] = (2.0f * random_unit(seed) - 1.0f) * bound;
}

Unet *unet_create(const UnetConfig *config) {
  if (!config_valid(config)) return NULL;
  Unet *unet = xcalloc(1, sizeof(*unet));
  unet->config = *config;
  for (int s = 0; s < ANNO_CWT_CHANNELS; s++) unet->cwt_std[s] = 1.0f;
  unet_layout(unet);
  unet->weights = xcalloc(unet->parameters, sizeof(float));
  uint64_t seed = ANNO_SEED * 0x9E3779B97F4A7C15ULL + 1;
  for (int l = 0; l < config->depth; l++) {
    int in = level_input(config, l), c = config->channels[l], k = config->kernel;
    init_conv(unet->weights + unet->offset[l][CONV_IN], in, c, k, &seed);
    init_conv(unet->weights + unet->offset[l][CONV_MID], c, c, k, &seed);
    if (l + 1 < config->depth) {
      init_conv(unet->weights + unet->offset[l][CONV_UP],
                c + config->channels[l + 1], c, k, &seed);
      init_conv(unet->weights + unet->offset[l][CONV_OUT], c, c, k, &seed);
    }
  }
  init_conv(unet->weights + unet->head, config->channels[0], ANNO_CLASSES, 1, &seed);
  return unet;
}

void unet_destroy(Unet *unet) {
  if (!unet) return;
  free(unet->weights);
  free(unet);
}

/**
 * @brief Allocate a workspace; `training` adds gradient buffers and
 *        `threads` workers share each convolution of this workspace.
 */
static Workspace *workspace_create(const Unet *unet, int64_t capacity,
                                   int training, int threads) {
  const UnetConfig *config = &unet->config;
  Workspace *w = xcalloc(1, sizeof(*w));
  w->capacity = capacity;
  w->threads = threads < 1 ? 1 : threads;
  size_t n = (size_t)capacity;
  w->input = xcalloc(mul_size(n, ANNO_INPUT_CHANNELS), sizeof(float));
  w->logits = xcalloc(mul_size(n, ANNO_CLASSES), sizeof(float));
  w->probability = xcalloc(mul_size(n, ANNO_CLASSES), sizeof(float));
  size_t largest = conv_size(config->channels[0], ANNO_CLASSES, 1);
  for (int l = 0; l < config->depth; l++) {
    int c = config->channels[l], next = l + 1 < config->depth ? config->channels[l + 1] : 0;
    size_t sizes[3] = {conv_size(level_input(config, l), c, config->kernel),
                       conv_size(c, c, config->kernel),
                       conv_size(c + next, c, config->kernel)};
    for (int i = 0; i < 3; i++) if (sizes[i] > largest) largest = sizes[i];
    w->a[l] = xcalloc(mul_size(n, (size_t)c), sizeof(float));
    w->b[l] = xcalloc(mul_size(n, (size_t)c), sizeof(float));
    if (next) {
      w->pooled[l] = xcalloc(mul_size(n / (size_t)config->pool, (size_t)c), sizeof(float));
      w->cat[l] = xcalloc(mul_size(n, (size_t)(c + next)), sizeof(float));
      w->c[l] = xcalloc(mul_size(n, (size_t)c), sizeof(float));
      w->d[l] = xcalloc(mul_size(n, (size_t)c), sizeof(float));
    }
    if (training) {
      w->ga[l] = xcalloc(mul_size(n, (size_t)c), sizeof(float));
      w->gb[l] = xcalloc(mul_size(n, (size_t)c), sizeof(float));
      w->gout[l] = xcalloc(mul_size(n, (size_t)c), sizeof(float));
      if (l) w->gin[l] = xcalloc(mul_size(n, (size_t)level_input(config, l)), sizeof(float));
      if (next) {
        w->gcat[l] = xcalloc(mul_size(n, (size_t)(c + next)), sizeof(float));
        w->gc[l] = xcalloc(mul_size(n, (size_t)c), sizeof(float));
      }
    }
    n /= (size_t)config->pool;
  }
  if (training) {
    w->glogits = xcalloc(mul_size((size_t)capacity, ANNO_CLASSES), sizeof(float));
    w->gradient = xcalloc(unet->parameters, sizeof(float));
    w->labels = xcalloc((size_t)capacity, 1);
    w->partial = xcalloc(mul_size(largest, (size_t)w->threads), sizeof(float));
  }
  return w;
}

/** @brief Free a workspace and every buffer it owns. */
static void workspace_destroy(Workspace *w) {
  if (!w) return;
  for (int l = 0; l < ANNO_MAX_DEPTH; l++) {
    free(w->a[l]); free(w->b[l]); free(w->pooled[l]); free(w->cat[l]);
    free(w->c[l]); free(w->d[l]); free(w->ga[l]); free(w->gb[l]);
    free(w->gin[l]); free(w->gcat[l]); free(w->gc[l]); free(w->gout[l]);
  }
  free(w->input); free(w->logits); free(w->probability); free(w->glogits);
  free(w->gradient); free(w->partial); free(w->labels);
  free(w);
}

/** Arguments of one convolution split into position ranges. */
typedef struct {
  const float *input, *output, *weights;
  float *result, *grad, *gweights, *ginput, *partial;
  int64_t n, chunk;
  int in, out, kernel, relu;
} Conv;

/** @brief Number of position chunks a convolution is split into. */
static int64_t conv_chunks(const Workspace *w, int64_t n, int64_t *chunk) {
  int64_t jobs = w->threads > 1 && n >= 4096 ? w->threads : 1;
  *chunk = (n + jobs - 1) / jobs;
  return jobs;
}

/** @brief Job: forward convolution of positions [job*chunk, +chunk). */
static void conv_forward_job(void *pointer, int64_t job, int worker) {
  (void)worker;
  const Conv *v = pointer;
  int radius = v->kernel / 2, in = v->in, out = v->out;
  const float *bias = v->weights + (size_t)v->kernel * (size_t)in * (size_t)out;
  int64_t end = (job + 1) * v->chunk < v->n ? (job + 1) * v->chunk : v->n;
  for (int64_t p = job * v->chunk; p < end; p++) {
    float *o = v->result + (size_t)p * (size_t)out;
    memcpy(o, bias, (size_t)out * sizeof(float));
    for (int k = 0; k < v->kernel; k++) {
      int64_t q = p + k - radius;
      if (q < 0 || q >= v->n) continue; /* explicit zero padding */
      const float *x = v->input + (size_t)q * (size_t)in;
      const float *wk = v->weights + (size_t)k * (size_t)in * (size_t)out;
      for (int c = 0; c < in; c++) {
        float value = x[c];
        if (value == 0.0f) continue;
        const float *wr = wk + (size_t)c * (size_t)out;
        for (int j = 0; j < out; j++) o[j] += value * wr[j];
      }
    }
    if (v->relu)
      for (int j = 0; j < out; j++) o[j] = o[j] > 0.0f ? o[j] : 0.0f;
  }
}

/**
 * @brief Same-padded 1D convolution, weights [k][in][out] then bias[out];
 *        positions outside [0,n) are zero. Split over w->threads.
 */
static void conv_forward(const Workspace *w, const float *input, int64_t n, int in,
                         int out, int kernel, const float *weights, float *output,
                         int relu) {
  Conv v = {.input = input, .weights = weights, .result = output, .n = n,
            .in = in, .out = out, .kernel = kernel, .relu = relu};
  run_parallel(w->threads, conv_chunks(w, n, &v.chunk), conv_forward_job, &v);
}

/**
 * @brief Job: for positions [job*chunk, +chunk) mask the ReLU gradient,
 *        accumulate weight gradients into this job's partial buffer and
 *        gather input gradients (no write conflicts between jobs).
 */
static void conv_backward_job(void *pointer, int64_t job, int worker) {
  (void)worker;
  const Conv *v = pointer;
  int radius = v->kernel / 2, in = v->in, out = v->out;
  size_t per = (size_t)in * (size_t)out;
  int64_t start = job * v->chunk, end = start + v->chunk < v->n ? start + v->chunk : v->n;
  float *gweights = v->partial + (size_t)job * (per * (size_t)v->kernel + (size_t)out);
  float *gbias = gweights + (size_t)v->kernel * per;
  if (v->relu)
    for (size_t i = (size_t)start * (size_t)out; i < (size_t)end * (size_t)out; i++)
      if (v->output[i] <= 0.0f) v->grad[i] = 0.0f;
  for (int64_t p = start; p < end; p++) {
    const float *g = v->grad + (size_t)p * (size_t)out;
    int any = 0;
    for (int j = 0; j < out; j++) {
      gbias[j] += g[j];
      any |= g[j] != 0.0f;
    }
    if (!any) continue;
    for (int k = 0; k < v->kernel; k++) {
      int64_t q = p + k - radius;
      if (q < 0 || q >= v->n) continue;
      const float *x = v->input + (size_t)q * (size_t)in;
      float *gw = gweights + (size_t)k * per;
      for (int c = 0; c < in; c++) {
        float value = x[c];
        if (value == 0.0f) continue;
        float *row = gw + (size_t)c * (size_t)out;
        for (int j = 0; j < out; j++) row[j] += value * g[j];
      }
    }
  }
}

/** @brief Job: gather dL/dinput for positions [job*chunk, +chunk). */
static void conv_input_job(void *pointer, int64_t job, int worker) {
  (void)worker;
  const Conv *v = pointer;
  int radius = v->kernel / 2, in = v->in, out = v->out;
  int64_t start = job * v->chunk, end = start + v->chunk < v->n ? start + v->chunk : v->n;
  for (int64_t q = start; q < end; q++) {
    float *gi = v->ginput + (size_t)q * (size_t)in;
    for (int k = 0; k < v->kernel; k++) {
      int64_t p = q - k + radius; /* output position that read q through tap k */
      if (p < 0 || p >= v->n) continue;
      const float *g = v->grad + (size_t)p * (size_t)out;
      const float *wk = v->weights + (size_t)k * (size_t)in * (size_t)out;
      for (int c = 0; c < in; c++) {
        const float *wr = wk + (size_t)c * (size_t)out;
        float sum = 0.0f;
        for (int j = 0; j < out; j++) sum += g[j] * wr[j];
        gi[c] += sum;
      }
    }
  }
}

/**
 * @brief Backward pass of conv_forward(). `grad` (dL/doutput) is masked in
 *        place by the ReLU; parameter gradients accumulate into `gweights`
 *        and input gradients accumulate into `ginput` when non-NULL.
 */
static void conv_backward(const Workspace *w, const float *input, const float *output,
                          float *grad, int64_t n, int in, int out, int kernel,
                          int relu, const float *weights, float *gweights,
                          float *ginput) {
  size_t size = conv_size(in, out, kernel);
  Conv v = {.input = input, .output = output, .weights = weights, .grad = grad,
            .ginput = ginput, .partial = w->partial, .n = n, .in = in, .out = out,
            .kernel = kernel, .relu = relu};
  int64_t jobs = conv_chunks(w, n, &v.chunk);
  memset(w->partial, 0, (size_t)jobs * size * sizeof(float));
  run_parallel(w->threads, jobs, conv_backward_job, &v);
  for (int64_t job = 0; job < jobs; job++) /* fixed order: deterministic sums */
    for (size_t i = 0; i < size; i++) gweights[i] += w->partial[(size_t)job * size + i];
  if (ginput) run_parallel(w->threads, jobs, conv_input_job, &v);
}

/** @brief Forward pass of level l on n positions; returns its output. */
static float *level_forward(const Unet *unet, Workspace *w, int level,
                            const float *x, int64_t n) {
  const UnetConfig *config = &unet->config;
  int in = level_input(config, level), c = config->channels[level];
  int k = config->kernel, p = config->pool;
  conv_forward(w, x, n, in, c, k, unet->weights + unet->offset[level][CONV_IN],
               w->a[level], 1);
  conv_forward(w, w->a[level], n, c, c, k,
               unet->weights + unet->offset[level][CONV_MID], w->b[level], 1);
  if (level + 1 >= config->depth || level + 1 >= ANNO_MAX_DEPTH) return w->b[level];
  int next = config->channels[level + 1];
  int64_t m = n / p;
  for (int64_t i = 0; i < m; i++)
    for (int j = 0; j < c; j++) {
      float sum = 0.0f;
      for (int t = 0; t < p; t++) sum += w->b[level][(size_t)(i * p + t) * (size_t)c + (size_t)j];
      w->pooled[level][(size_t)i * (size_t)c + (size_t)j] = sum / (float)p;
    }
  const float *u = level_forward(unet, w, level + 1, w->pooled[level], m);
  for (int64_t i = 0; i < n; i++) {
    float *row = w->cat[level] + (size_t)i * (size_t)(c + next);
    memcpy(row, w->b[level] + (size_t)i * (size_t)c, (size_t)c * sizeof(float));
    memcpy(row + c, u + (size_t)(i / p) * (size_t)next, (size_t)next * sizeof(float));
  }
  conv_forward(w, w->cat[level], n, c + next, c, k,
               unet->weights + unet->offset[level][CONV_UP], w->c[level], 1);
  conv_forward(w, w->c[level], n, c, c, k,
               unet->weights + unet->offset[level][CONV_OUT], w->d[level], 1);
  return w->d[level];
}

/** @brief Numerically stable softmax of one row. */
static void softmax_row(const float *logits, float *probability) {
  float top = logits[0];
  for (int j = 1; j < ANNO_CLASSES; j++) if (logits[j] > top) top = logits[j];
  double sum = 0.0;
  for (int j = 0; j < ANNO_CLASSES; j++) sum += exp((double)(logits[j] - top));
  for (int j = 0; j < ANNO_CLASSES; j++)
    probability[j] = (float)(exp((double)(logits[j] - top)) / sum);
}

/**
 * @brief Full forward pass on w->input (n positions, n a pooling multiple)
 *        producing w->logits and w->probability.
 * @return 0, or ERANGE when a non-finite value appears.
 */
static int unet_forward(const Unet *unet, Workspace *w, int64_t n) {
  const float *top = level_forward(unet, w, 0, w->input, n);
  conv_forward(w, top, n, unet->config.channels[0], ANNO_CLASSES, 1,
               unet->weights + unet->head, w->logits, 0);
  for (int64_t i = 0; i < n; i++) {
    const float *row = w->logits + (size_t)i * ANNO_CLASSES;
    for (int j = 0; j < ANNO_CLASSES; j++)
      if (!isfinite(row[j])) return ERANGE;
    softmax_row(row, w->probability + (size_t)i * ANNO_CLASSES);
  }
  return 0;
}

/** @brief Backward pass of level l; gout[level] holds dL/doutput. */
static void level_backward(const Unet *unet, Workspace *w, int level,
                           const float *x, int64_t n, float *gx) {
  const UnetConfig *config = &unet->config;
  int in = level_input(config, level), c = config->channels[level];
  int k = config->kernel, p = config->pool;
  float *gw = w->gradient;
  float *gb = w->gout[level];
  if (level + 1 < config->depth && level + 1 < ANNO_MAX_DEPTH) {
    int next = config->channels[level + 1];
    int64_t m = n / p;
    memset(w->gc[level], 0, (size_t)n * (size_t)c * sizeof(float));
    conv_backward(w, w->c[level], w->d[level], w->gout[level], n, c, c, k, 1,
                  unet->weights + unet->offset[level][CONV_OUT],
                  gw + unet->offset[level][CONV_OUT], w->gc[level]);
    memset(w->gcat[level], 0, (size_t)n * (size_t)(c + next) * sizeof(float));
    conv_backward(w, w->cat[level], w->c[level], w->gc[level], n, c + next, c, k, 1,
                  unet->weights + unet->offset[level][CONV_UP],
                  gw + unet->offset[level][CONV_UP], w->gcat[level]);
    float *gu = w->gout[level + 1];
    memset(gu, 0, (size_t)m * (size_t)next * sizeof(float));
    gb = w->gb[level];
    for (int64_t i = 0; i < n; i++) {
      const float *row = w->gcat[level] + (size_t)i * (size_t)(c + next);
      memcpy(gb + (size_t)i * (size_t)c, row, (size_t)c * sizeof(float));
      float *u = gu + (size_t)(i / p) * (size_t)next;
      for (int j = 0; j < next; j++) u[j] += row[c + j];
    }
    memset(w->gin[level + 1], 0, (size_t)m * (size_t)c * sizeof(float));
    level_backward(unet, w, level + 1, w->pooled[level], m, w->gin[level + 1]);
    for (int64_t i = 0; i < n; i++) {
      const float *gp = w->gin[level + 1] + (size_t)(i / p) * (size_t)c;
      for (int j = 0; j < c; j++) gb[(size_t)i * (size_t)c + (size_t)j] += gp[j] / (float)p;
    }
  }
  memset(w->ga[level], 0, (size_t)n * (size_t)c * sizeof(float));
  conv_backward(w, w->a[level], w->b[level], gb, n, c, c, k, 1,
                unet->weights + unet->offset[level][CONV_MID],
                gw + unet->offset[level][CONV_MID], w->ga[level]);
  conv_backward(w, x, w->a[level], w->ga[level], n, in, c, k, 1,
                unet->weights + unet->offset[level][CONV_IN],
                gw + unet->offset[level][CONV_IN], gx);
}

/**
 * @brief ANNEVO-style loss of one window: focal cross entropy (gamma 2)
 *        plus per-class soft-F1 (dice) of every non-intergenic class; a
 *        class absent from the window contributes its mean false-positive
 *        probability. Writes scale * dL/dlogits to w->glogits.
 * @param valid labelled positions; padding beyond it gets zero gradient.
 * @return unscaled loss.
 */
static double window_loss(Workspace *w, int64_t valid, int64_t n, float scale) {
  const float *prob = w->probability;
  const uint8_t *y = w->labels;
  double intersection[ANNO_CLASSES] = {0}, sum_p[ANNO_CLASSES] = {0};
  double sum_y[ANNO_CLASSES] = {0}, fp_sum[ANNO_CLASSES] = {0};
  double fp_count[ANNO_CLASSES] = {0};
  for (int64_t i = 0; i < valid; i++) {
    const float *row = prob + (size_t)i * ANNO_CLASSES;
    int top = 0;
    for (int j = 0; j < ANNO_CLASSES; j++) {
      sum_p[j] += row[j];
      if (row[j] > row[top]) top = j;
    }
    intersection[y[i]] += row[y[i]];
    sum_y[y[i]] += 1.0;
    if (top != y[i]) {
      fp_sum[top] += row[top];
      fp_count[top] += 1.0;
    }
  }
  double loss = 0.0, slope[ANNO_CLASSES] = {0}, offset[ANNO_CLASSES] = {0};
  double fp_weight[ANNO_CLASSES] = {0};
  for (int j = 1; j < ANNO_CLASSES; j++) {
    if (sum_y[j] > 0) {
      double u = sum_p[j] + sum_y[j] + ANNO_DICE_EPSILON;
      double top = 2.0 * intersection[j] + ANNO_DICE_EPSILON;
      loss += ANNO_DICE_WEIGHT * (1.0 - top / u);
      slope[j] = -2.0 / u;           /* times y_i */
      offset[j] = top / (u * u);     /* every position */
    } else {
      double count = fp_count[j] > 1 ? fp_count[j] : 1.0;
      loss += ANNO_DICE_WEIGHT * fp_sum[j] / count;
      fp_weight[j] = 1.0 / count;
    }
  }
  memset(w->glogits, 0, (size_t)n * ANNO_CLASSES * sizeof(float));
  for (int64_t i = 0; i < valid; i++) {
    const float *row = prob + (size_t)i * ANNO_CLASSES;
    float *g = w->glogits + (size_t)i * ANNO_CLASSES;
    int t = y[i], top = 0;
    for (int j = 1; j < ANNO_CLASSES; j++) if (row[j] > row[top]) top = j;
    double pt = row[t] > 1e-12f ? row[t] : 1e-12;
    double miss = 1.0 - pt;
    loss += -miss * miss * log(pt) / (double)valid;
    /* d/dz of -(1-pt)^2 log pt: (2(1-pt) pt log pt - (1-pt)^2)(delta - p) */
    double a = (ANNO_FOCAL_GAMMA * miss * pt * log(pt) - miss * miss) / (double)valid;
    double gp[ANNO_CLASSES], dot = 0.0;
    for (int j = 0; j < ANNO_CLASSES; j++) {
      gp[j] = ANNO_DICE_WEIGHT * (slope[j] * (t == j) + offset[j] +
                                  (top == j && top != t ? fp_weight[j] : 0.0));
      dot += gp[j] * row[j];
    }
    for (int j = 0; j < ANNO_CLASSES; j++)
      g[j] = scale * (float)(a * ((t == j) - row[j]) + row[j] * (gp[j] - dot));
  }
  return loss;
}

/** @brief Backward pass from w->glogits into w->gradient (accumulating). */
static void unet_backward(const Unet *unet, Workspace *w, int64_t n) {
  int c = unet->config.channels[0];
  const float *top = unet->config.depth > 1 ? w->d[0] : w->b[0];
  memset(w->gout[0], 0, (size_t)n * (size_t)c * sizeof(float));
  conv_backward(w, top, w->logits, w->glogits, n, c, ANNO_CLASSES, 1, 0,
                unet->weights + unet->head, w->gradient + unet->head,
                w->gout[0]);
  level_backward(unet, w, 0, w->input, n, NULL);
}

int unet_predict(const Unet *unet, const float *input, int64_t length,
                 float *probability) {
  if (!unet || !input || !probability || length < 1) return EINVAL;
  int64_t multiple = unet_multiple(&unet->config);
  if (length > INT64_MAX - multiple) return EOVERFLOW;
  int64_t n = (length + multiple - 1) / multiple * multiple;
  Workspace *w = workspace_create(unet, n, 0, 1);
  memcpy(w->input, input, (size_t)length * ANNO_INPUT_CHANNELS * sizeof(float));
  int status = unet_forward(unet, w, n);
  if (!status)
    memcpy(probability, w->probability, (size_t)length * ANNO_CLASSES * sizeof(float));
  workspace_destroy(w);
  return status;
}

/* ---- model file ---- */

static const char model_magic[8] = {'A', 'N', 'N', 'O', 'U', 'N', 'E', 'T'};
static const uint32_t model_version = 1;

/** @brief FNV-1a checksum. */
static uint32_t checksum(uint32_t hash, const void *data, size_t size) {
  const unsigned char *bytes = data;
  for (size_t i = 0; i < size; i++) hash = (hash ^ bytes[i]) * 16777619u;
  return hash;
}

/** @brief fwrite() a block and fold it into the checksum. */
static int write_block(FILE *file, const void *data, size_t size, uint32_t *hash) {
  *hash = checksum(*hash, data, size);
  return fwrite(data, 1, size, file) == size;
}

/** @brief fread() a block and fold it into the checksum. */
static int read_block(FILE *file, void *data, size_t size, uint32_t *hash) {
  if (fread(data, 1, size, file) != size) return 0;
  *hash = checksum(*hash, data, size);
  return 1;
}

int unet_save(const char *path, const Unet *unet) {
  if (!path || !unet) return EINVAL;
  FILE *file = fopen(path, "wb");
  if (!file) return errno ? errno : EIO;
  uint32_t hash = 2166136261u, endian = 0x01020304u;
  int32_t shape[3 + ANNO_MAX_DEPTH] = {unet->config.depth, unet->config.kernel,
                                       unet->config.pool};
  for (int l = 0; l < ANNO_MAX_DEPTH; l++) shape[3 + l] = unet->config.channels[l];
  uint64_t count = unet->parameters;
  int ok = write_block(file, model_magic, sizeof(model_magic), &hash) &&
           write_block(file, &model_version, sizeof(model_version), &hash) &&
           write_block(file, &endian, sizeof(endian), &hash) &&
           write_block(file, shape, sizeof(shape), &hash) &&
           write_block(file, unet->cwt_mean, sizeof(unet->cwt_mean), &hash) &&
           write_block(file, unet->cwt_std, sizeof(unet->cwt_std), &hash) &&
           write_block(file, &count, sizeof(count), &hash) &&
           write_block(file, unet->weights, unet->parameters * sizeof(float), &hash);
  ok = ok && fwrite(&hash, sizeof(hash), 1, file) == 1;
  if (fclose(file)) ok = 0;
  return ok ? 0 : EIO;
}

int unet_load(const char *path, Unet **out) {
  if (!path || !out) return EINVAL;
  *out = NULL;
  FILE *file = fopen(path, "rb");
  if (!file) return errno ? errno : EIO;
  uint32_t hash = 2166136261u, version = 0, endian = 0, stored = 0;
  char magic[8];
  int32_t shape[3 + ANNO_MAX_DEPTH];
  float mean[ANNO_CWT_CHANNELS], deviation[ANNO_CWT_CHANNELS];
  uint64_t count = 0;
  int ok = read_block(file, magic, sizeof(magic), &hash) &&
           !memcmp(magic, model_magic, sizeof(magic)) &&
           read_block(file, &version, sizeof(version), &hash) &&
           version == model_version &&
           read_block(file, &endian, sizeof(endian), &hash) &&
           endian == 0x01020304u &&
           read_block(file, shape, sizeof(shape), &hash) &&
           read_block(file, mean, sizeof(mean), &hash) &&
           read_block(file, deviation, sizeof(deviation), &hash) &&
           read_block(file, &count, sizeof(count), &hash);
  UnetConfig config = {0};
  if (ok) {
    config.depth = shape[0];
    config.kernel = shape[1];
    config.pool = shape[2];
    for (int l = 0; l < ANNO_MAX_DEPTH; l++)
      config.channels[l] = l < config.depth ? shape[3 + l] : 0;
  }
  Unet *unet = ok ? unet_create(&config) : NULL;
  ok = unet && count == unet->parameters &&
       read_block(file, unet->weights, unet->parameters * sizeof(float), &hash) &&
       fread(&stored, sizeof(stored), 1, file) == 1 && stored == hash;
  for (int s = 0; ok && s < ANNO_CWT_CHANNELS; s++) {
    ok = isfinite(mean[s]) && isfinite(deviation[s]) && deviation[s] > 0;
    if (ok) {
      unet->cwt_mean[s] = mean[s];
      unet->cwt_std[s] = deviation[s];
    }
  }
  for (size_t i = 0; ok && i < unet->parameters; i++) ok = isfinite(unet->weights[i]);
  fclose(file);
  if (!ok) {
    unet_destroy(unet);
    return EINVAL;
  }
  *out = unet;
  return 0;
}

/* ------------------------------------------------------------------------ */
/* Training                                                                 */
/* ------------------------------------------------------------------------ */

/** @brief Map a class to the same feature on the opposite strand. */
static uint8_t swap_strand(uint8_t label) {
  return label == ANNO_INTERGENIC ? label : label <= 4 ? label + 4 : label - 4;
}

/**
 * @brief Write normalised features of bases[start, start+count) into
 *        `input`; with `flip` the window is reverse-complemented using the
 *        mirrored (strand-symmetric) CWT power instead of a new transform.
 */
static void fill_input(const Unet *unet, const AnnoSequence *sequence,
                       const float *power, int64_t start, int64_t count,
                       int flip, float *input) {
  for (int64_t q = 0; q < count; q++) {
    int64_t p = flip ? start + count - 1 - q : start + q;
    float *row = input + (size_t)q * ANNO_INPUT_CHANNELS;
    int b = base_index(sequence->bases[p]);
    if (flip && b >= 0) b = 3 - b;
    for (int j = 0; j < 4; j++) row[j] = (float)(b == j);
    for (int s = 0; s < ANNO_CWT_CHANNELS; s++)
      row[4 + s] = (power[(size_t)p * ANNO_CWT_CHANNELS + (size_t)s] -
                    unet->cwt_mean[s]) / unet->cwt_std[s];
  }
}

/** One sampled training window. */
typedef struct {
  int sequence;
  int64_t start, valid;
  int flip;
} Sample;

/** Shared context of a training step. */
typedef struct {
  const Unet *unet;
  const AnnoSequence *sequences;
  float *const *power;
  const Sample *samples;
  Workspace **workspaces;
  double *losses;
  int failed;
} TrainStep;

/** @brief Job: forward + loss + backward of one sampled window. */
static void train_job(void *pointer, int64_t job, int worker) {
  TrainStep *step = pointer;
  Workspace *w = step->workspaces[worker];
  const Sample *s = &step->samples[job];
  const AnnoSequence *sequence = &step->sequences[s->sequence];
  int64_t multiple = unet_multiple(&step->unet->config);
  int64_t n = (s->valid + multiple - 1) / multiple * multiple;
  memset(w->input, 0, (size_t)n * ANNO_INPUT_CHANNELS * sizeof(float));
  fill_input(step->unet, sequence, step->power[s->sequence], s->start, s->valid,
             s->flip, w->input);
  for (int64_t q = 0; q < s->valid; q++) {
    uint8_t label = sequence->labels[s->flip ? s->start + s->valid - 1 - q : s->start + q];
    w->labels[q] = s->flip ? swap_strand(label) : label;
  }
  if (unet_forward(step->unet, w, n)) {
    step->failed = 1;
    return;
  }
  step->losses[job] = window_loss(w, s->valid, n, 1.0f / ANNO_BATCH_WINDOWS);
  unet_backward(step->unet, w, n);
}

/** @brief Mean/std of CWT power and class priors (head bias) of the training data. */
static void compute_normalisation(Unet *unet, const AnnoSequence *sequences,
                                  float *const *power, int count) {
  double sum[ANNO_CWT_CHANNELS] = {0}, square[ANNO_CWT_CHANNELS] = {0}, total = 0;
  for (int i = 0; i < count; i++)
    for (int64_t p = 0; p < sequences[i].length; p++, total++)
      for (int s = 0; s < ANNO_CWT_CHANNELS; s++) {
        double v = power[i][(size_t)p * ANNO_CWT_CHANNELS + (size_t)s];
        sum[s] += v;
        square[s] += v * v;
      }
  for (int s = 0; s < ANNO_CWT_CHANNELS; s++) {
    double mean = sum[s] / total, variance = square[s] / total - mean * mean;
    unet->cwt_mean[s] = (float)mean;
    unet->cwt_std[s] = variance > 1e-12 ? (float)sqrt(variance) : 1.0f;
  }
  /* head bias = log class frequency (strand-symmetrised for the flips) */
  double frequency[ANNO_CLASSES] = {0};
  for (int i = 0; i < count; i++)
    for (int64_t p = 0; p < sequences[i].length; p++) {
      frequency[sequences[i].labels[p]] += 0.5;
      frequency[swap_strand(sequences[i].labels[p])] += 0.5;
    }
  float *bias = unet->weights + unet->head + (size_t)unet->config.channels[0] * ANNO_CLASSES;
  for (int j = 0; j < ANNO_CLASSES; j++) bias[j] = (float)log((frequency[j] + 1.0) / (total + 1.0));
}

/**
 * @brief Train a UNet with mini-batches of random (optionally reverse-
 *        complemented) windows and Adam. Returns the checkpoint with the
 *        lowest smoothed training loss. Caller owns the result.
 */
static Unet *train_unet(const AnnoSequence *sequences, int count, int iterations,
                        int threads) {
  UnetConfig config = unet_default_config();
  Unet *unet = unet_create(&config);
  if (!unet) fail("invalid UNet configuration");
  double started = now_seconds();
  float **power = xcalloc((size_t)count, sizeof(*power));
  int64_t total = 0;
  for (int i = 0; i < count; i++) {
    power[i] = sequence_power(&sequences[i], threads);
    total += sequences[i].length;
  }
  compute_normalisation(unet, sequences, power, count);
  say("CWT features ready in %.1fs (%lld bp)", now_seconds() - started,
      (long long)total);
  /* one workspace; every convolution is split over all threads */
  Workspace *workspace = workspace_create(unet, ANNO_UNET_WINDOW, 1, threads);
  size_t parameters = unet->parameters;
  float *m = xcalloc(parameters, sizeof(float)), *v = xcalloc(parameters, sizeof(float));
  float *gradient = xcalloc(parameters, sizeof(float));
  float *best = xcalloc(parameters, sizeof(float));
  memcpy(best, unet->weights, parameters * sizeof(float));
  double best_loss = INFINITY, smooth = NAN;
  uint64_t seed = ANNO_SEED;
  Sample samples[ANNO_BATCH_WINDOWS];
  double losses[ANNO_BATCH_WINDOWS];
  for (int iteration = 1; iteration <= iterations; iteration++) {
    for (int b = 0; b < ANNO_BATCH_WINDOWS; b++) {
      int64_t pick = random_below(&seed, total);
      int sequence = 0;
      while (pick >= sequences[sequence].length) pick -= sequences[sequence++].length;
      int64_t valid = sequences[sequence].length < ANNO_UNET_WINDOW
        ? sequences[sequence].length : ANNO_UNET_WINDOW;
      samples[b] = (Sample){sequence,
                            random_below(&seed, sequences[sequence].length - valid + 1),
                            valid, (int)(random_next(&seed) >> 63)};
    }
    memset(workspace->gradient, 0, parameters * sizeof(float));
    TrainStep step = {unet, sequences, power, samples, &workspace, losses, 0};
    for (int b = 0; b < ANNO_BATCH_WINDOWS && !step.failed; b++) train_job(&step, b, 0);
    if (step.failed) fail("non-finite UNet activation during training");
    double loss = 0.0, norm = 0.0;
    for (int b = 0; b < ANNO_BATCH_WINDOWS; b++) loss += losses[b] / ANNO_BATCH_WINDOWS;
    for (size_t k = 0; k < parameters; k++) {
      float g = workspace->gradient[k];
      gradient[k] = g;
      norm += (double)g * g;
    }
    norm = sqrt(norm);
    if (!isfinite(norm) || !isfinite(loss)) fail("non-finite gradient during training");
    float clip = norm > ANNO_GRADIENT_CLIP ? (float)(ANNO_GRADIENT_CLIP / norm) : 1.0f;
    double progress = (double)(iteration - 1) / (double)iterations;
    float rate = ANNO_LEARNING_RATE *
                 (float)(iteration <= 20 ? iteration / 20.0
                         : 0.05 + 0.95 * 0.5 * (1.0 + cos(ANNO_PI * progress)));
    float correction1 = 1.0f - powf(0.9f, (float)iteration);
    float correction2 = 1.0f - powf(0.999f, (float)iteration);
    for (size_t k = 0; k < parameters; k++) {
      float g = gradient[k] * clip;
      m[k] = 0.9f * m[k] + 0.1f * g;
      v[k] = 0.999f * v[k] + 0.001f * g * g;
      unet->weights[k] -= rate * (m[k] / correction1) /
                          (sqrtf(v[k] / correction2) + 1e-8f);
    }
    smooth = isnan(smooth) ? loss : 0.9 * smooth + 0.1 * loss;
    if (iteration >= iterations / 10 && smooth < best_loss) {
      best_loss = smooth;
      memcpy(best, unet->weights, parameters * sizeof(float));
    }
    if (iteration == 1 || iteration % 10 == 0 || iteration == iterations)
      say("iteration %d/%d loss %.4f smooth %.4f |g| %.3f lr %.5f %.1fs", iteration,
          iterations, loss, smooth, norm, rate, now_seconds() - started);
  }
  memcpy(unet->weights, best, parameters * sizeof(float));
  say("best smoothed training loss %.4f", best_loss);
  workspace_destroy(workspace);
  for (int i = 0; i < count; i++) free(power[i]);
  free(power); free(m); free(v); free(gradient); free(best);
  return unet;
}

/* ------------------------------------------------------------------------ */
/* Grammar Viterbi                                                          */
/* ------------------------------------------------------------------------ */

#define GRAMMAR_MAX_STATES 256
#define GRAMMAR_MAX_PREDECESSORS 5

enum { MOTIF_NONE, MOTIF_PLUS_DONOR, MOTIF_PLUS_ACCEPTOR, MOTIF_MINUS_ENTRY,
       MOTIF_MINUS_EXIT };

/**
 * Expanded grammar: every class state is a chain of sub-states enforcing a
 * minimum length; only the chain tail self-loops.
 */
typedef struct {
  int count;
  uint8_t cls[GRAMMAR_MAX_STATES];
  int predecessors[GRAMMAR_MAX_STATES];
  int from[GRAMMAR_MAX_STATES][GRAMMAR_MAX_PREDECESSORS];
  double score[GRAMMAR_MAX_STATES][GRAMMAR_MAX_PREDECESSORS];
  uint8_t motif[GRAMMAR_MAX_STATES][GRAMMAR_MAX_PREDECESSORS];
  int shift[GRAMMAR_MAX_STATES], bits[GRAMMAR_MAX_STATES];
} Grammar;

/** @brief Add a transition from -> to. */
static void grammar_edge(Grammar *g, int from, int to, double score, int motif) {
  int k = g->predecessors[to]++;
  if (k >= GRAMMAR_MAX_PREDECESSORS) fail("grammar has too many predecessors");
  g->from[to][k] = from;
  g->score[to][k] = score;
  g->motif[to][k] = (uint8_t)motif;
}

/** @brief Add a chain of `length` sub-states; returns the head index. */
static int grammar_chain(Grammar *g, uint8_t cls, int length) {
  int head = g->count;
  for (int i = 0; i < length; i++) {
    g->cls[g->count] = cls;
    if (i) grammar_edge(g, g->count - 1, g->count, 0.0, MOTIF_NONE);
    g->count++;
  }
  grammar_edge(g, g->count - 1, g->count - 1, 0.0, MOTIF_NONE);
  return head;
}

/**
 * @brief Build the strand-aware gene grammar. Plus strand in genomic order:
 *   intergenic -> [5'UTR (intron 5'UTR)*] -> exon (intron exon)* ->
 *   [3'UTR (intron 3'UTR)*] -> intergenic; minus strand is the mirror image.
 *   Introns are never entered from intergenic and UTRs never touch directly.
 */
static void grammar_build(Grammar *g) {
  memset(g, 0, sizeof(*g));
  int inter = grammar_chain(g, ANNO_INTERGENIC, 1);
  for (int strand = 0; strand < 2; strand++) {
    int o = strand ? 4 : 0, donor = strand ? MOTIF_MINUS_ENTRY : MOTIF_PLUS_DONOR;
    int acceptor = strand ? MOTIF_MINUS_EXIT : MOTIF_PLUS_ACCEPTOR;
    int u5 = grammar_chain(g, (uint8_t)(ANNO_PLUS_5UTR + o), 1);
    int ex = grammar_chain(g, (uint8_t)(ANNO_PLUS_EXON + o), ANNO_MIN_EXON_LENGTH);
    int ex_tail = ex + ANNO_MIN_EXON_LENGTH - 1;
    int u3 = grammar_chain(g, (uint8_t)(ANNO_PLUS_3UTR + o), 1);
    int introns[3];
    for (int i = 0; i < 3; i++)
      introns[i] = grammar_chain(g, (uint8_t)(ANNO_PLUS_INTRON + o),
                                 ANNO_MIN_INTRON_LENGTH);
    int tail = ANNO_MIN_INTRON_LENGTH - 1;
    /* exonic pieces in genomic order: first, middle (coding), last */
    int first = strand ? u3 : u5, last = strand ? u5 : u3;
    grammar_edge(g, inter, first, -ANNO_TRANSCRIPT_PENALTY, MOTIF_NONE);
    grammar_edge(g, inter, ex, -ANNO_TRANSCRIPT_PENALTY, MOTIF_NONE);
    grammar_edge(g, first, ex, 0.0, MOTIF_NONE);
    grammar_edge(g, ex_tail, last, 0.0, MOTIF_NONE);
    grammar_edge(g, ex_tail, inter, 0.0, MOTIF_NONE);
    grammar_edge(g, last, inter, 0.0, MOTIF_NONE);
    int pieces[3] = {first, ex, last}, exits[3] = {first, ex_tail, last};
    for (int i = 0; i < 3; i++) {
      grammar_edge(g, exits[i], introns[i], 0.0, donor);
      grammar_edge(g, introns[i] + tail, pieces[i], 0.0, acceptor);
    }
  }
  int shift = 0;
  for (int s = 0; s < g->count; s++) {
    int bits = 0;
    while ((1 << bits) < g->predecessors[s]) bits++;
    g->shift[s] = shift;
    g->bits[s] = bits;
    shift += bits;
  }
  if (shift > 32 || g->count > GRAMMAR_MAX_STATES) fail("grammar too large");
}

/** @brief Whether the transition motif holds when entering position i. */
static int motif_matches(int motif, const char *bases, int64_t length, int64_t i) {
  switch (motif) {
  case MOTIF_PLUS_DONOR: return i + 1 < length && bases[i] == 'G' && bases[i + 1] == 'T';
  case MOTIF_PLUS_ACCEPTOR: return i >= 2 && bases[i - 2] == 'A' && bases[i - 1] == 'G';
  case MOTIF_MINUS_ENTRY: return i + 1 < length && bases[i] == 'C' && bases[i + 1] == 'T';
  case MOTIF_MINUS_EXIT: return i >= 2 && bases[i - 2] == 'A' && bases[i - 1] == 'C';
  default: return 1;
  }
}

/** Streaming Viterbi decoder. Owns its buffers. */
typedef struct {
  Grammar grammar;
  const char *bases;
  int64_t length, position;
  double *previous, *current;
  uint32_t *back;
} Decoder;

/** @brief Create a decoder for one sequence. */
static Decoder *decoder_create(const char *bases, int64_t length) {
  Decoder *d = xcalloc(1, sizeof(*d));
  grammar_build(&d->grammar);
  d->bases = bases;
  d->length = length;
  d->previous = xcalloc((size_t)d->grammar.count, sizeof(double));
  d->current = xcalloc((size_t)d->grammar.count, sizeof(double));
  d->back = xcalloc((size_t)length, sizeof(uint32_t));
  return d;
}

/** @brief Free a decoder. */
static void decoder_destroy(Decoder *d) {
  if (!d) return;
  free(d->previous);
  free(d->current);
  free(d->back);
  free(d);
}

/** @brief Consume the probability row of the next position. */
static void decoder_step(Decoder *d, const float *row) {
  const Grammar *g = &d->grammar;
  double emission[ANNO_CLASSES];
  for (int j = 0; j < ANNO_CLASSES; j++) {
    float p = isfinite(row[j]) && row[j] > 0.0f ? row[j] : 0.0f;
    emission[j] = log((double)p + ANNO_EMISSION_EPSILON);
  }
  int64_t i = d->position++;
  if (!i) {
    for (int s = 0; s < g->count; s++) d->previous[s] = -INFINITY;
    d->previous[0] = emission[ANNO_INTERGENIC];
    return;
  }
  uint32_t back = 0;
  double top = -INFINITY;
  for (int s = 0; s < g->count; s++) {
    double best = -INFINITY;
    int choice = 0;
    for (int k = 0; k < g->predecessors[s]; k++) {
      double score = d->previous[g->from[s][k]] + g->score[s][k];
      if (g->motif[s][k] && !motif_matches(g->motif[s][k], d->bases, d->length, i))
        score -= ANNO_SPLICE_BONUS; /* canonical motifs keep the bonus */
      if (score > best) {
        best = score;
        choice = k;
      }
    }
    d->current[s] = best + emission[g->cls[s]];
    if (d->current[s] > top) top = d->current[s];
    back |= (uint32_t)choice << g->shift[s];
  }
  for (int s = 0; s < g->count; s++) d->current[s] -= top; /* keep scores small */
  d->back[i] = back;
  double *swap = d->previous;
  d->previous = d->current;
  d->current = swap;
}

/** @brief Trace back the best path ending in intergenic. @return 0 on success. */
static int decoder_finish(Decoder *d, uint8_t *path) {
  const Grammar *g = &d->grammar;
  if (d->position != d->length || !isfinite(d->previous[0])) return ERANGE;
  int state = 0;
  for (int64_t i = d->length - 1; i > 0; i--) {
    path[i] = g->cls[state];
    int choice = (int)((d->back[i] >> g->shift[state]) & ((1u << g->bits[state]) - 1u));
    state = g->from[state][choice];
  }
  path[0] = g->cls[state];
  return state == 0 ? 0 : ERANGE;
}

int anno_decode(const char *bases, int64_t length, const float *probability,
                uint8_t *path) {
  if (!bases || !probability || !path || length < 1) return EINVAL;
  Decoder *d = decoder_create(bases, length);
  for (int64_t i = 0; i < length; i++)
    decoder_step(d, probability + (size_t)i * ANNO_CLASSES);
  int status = decoder_finish(d, path);
  decoder_destroy(d);
  return status;
}

/* ------------------------------------------------------------------------ */
/* ORF and GFF3 output                                                      */
/* ------------------------------------------------------------------------ */

/** @brief Whether the 3 bases form a stop codon. */
static int is_stop(const char *c) {
  return c[0] == 'T' && ((c[1] == 'A' && (c[2] == 'A' || c[2] == 'G')) ||
                         (c[1] == 'G' && c[2] == 'A'));
}

/**
 * @brief Longest ATG..stop ORF (stop included) in any frame of a spliced
 *        transcript; ties keep the 5'-most ORF.
 * @return 1 when found, filling [*start,*end) in transcript coordinates.
 */
static int longest_orf(const char *transcript, int64_t length, int64_t *start,
                       int64_t *end) {
  int64_t best = 0;
  for (int frame = 0; frame < 3; frame++) {
    int64_t open = -1;
    for (int64_t p = frame; p + 3 <= length; p += 3) {
      const char *c = transcript + p;
      if (open < 0 && c[0] == 'A' && c[1] == 'T' && c[2] == 'G') open = p;
      if (!is_stop(c)) continue;
      if (open >= 0 && (p + 3 - open > best ||
                        (p + 3 - open == best && open < *start))) {
        best = p + 3 - open;
        *start = open;
        *end = p + 3;
      }
      open = -1;
    }
  }
  return best > 0;
}

/** @brief Map transcript range [from,to) onto genomic exon pieces. */
static void map_range(const SegmentList *exons, char strand, int64_t from,
                      int64_t to, SegmentList *out) {
  int64_t offset = 0;
  for (int order = 0; order < exons->count; order++) {
    const Segment *e = &exons->items[strand == '+' ? order : exons->count - 1 - order];
    int64_t length = e->end - e->start;
    int64_t a = from > offset ? from - offset : 0;
    int64_t b = to < offset + length ? to - offset : length;
    if (a < b)
      segments_push(out, strand == '+' ? e->start + a : e->end - b,
                    strand == '+' ? e->start + b : e->end - a);
    offset += length;
  }
  qsort(out->items, (size_t)out->count, sizeof(Segment), segment_compare);
}

/** One GFF3 record waiting for validation. */
typedef struct {
  const char *type;
  int64_t start, end;
  int phase;
  char id[96], parent[96];
} Record;

/** @brief Whether [start,end) lies within one of the segments. */
static int inside(const SegmentList *list, int64_t start, int64_t end) {
  for (int i = 0; i < list->count; i++)
    if (start >= list->items[i].start && end <= list->items[i].end) return 1;
  return 0;
}

/**
 * @brief Check ID/Parent links and containment of the records of one gene.
 * @return 1 when consistent.
 */
static int records_valid(const Record *records, int count, const SegmentList *exons) {
  for (int i = 0; i < count; i++) {
    const Record *r = &records[i];
    if (r->start < 0 || r->end <= r->start) return 0;
    if (i == 0) {
      if (strcmp(r->type, "gene") || r->parent[0]) return 0;
      continue;
    }
    int parent = -1;
    for (int j = 0; j < count; j++)
      if (r->parent[0] && !strcmp(records[j].id, r->parent)) parent = j;
    if (parent < 0 || r->start < records[parent].start || r->end > records[parent].end)
      return 0;
    if (strcmp(r->type, "mRNA") && strcmp(r->type, "exon") &&
        !inside(exons, r->start, r->end))
      return 0;
    for (int j = 0; j < i; j++)
      if (r->id[0] && !strcmp(records[j].id, r->id)) return 0;
  }
  return 1;
}

/** @brief Append a record. */
static void record_add(Record **records, int *count, const char *type,
                       int64_t start, int64_t end, int phase, const char *id,
                       const char *parent) {
  *records = xrealloc(*records, (size_t)*count + 1, sizeof(**records));
  Record *r = &(*records)[(*count)++];
  r->type = type;
  r->start = start;
  r->end = end;
  r->phase = phase;
  snprintf(r->id, sizeof(r->id), "%s", id);
  snprintf(r->parent, sizeof(r->parent), "%s", parent);
}

/**
 * @brief Emit one transcript (exons on one strand) with its longest ORF.
 * @return 1 written, 0 skipped (no ORF), -1 on error.
 */
static int write_transcript(FILE *output, const AnnoSequence *sequence,
                            const SegmentList *exons, char strand,
                            unsigned long gene_number) {
  int64_t length = 0;
  for (int i = 0; i < exons->count; i++) length += exons->items[i].end - exons->items[i].start;
  char *transcript = xcalloc((size_t)length + 1, 1);
  int64_t at = 0;
  for (int order = 0; order < exons->count; order++) {
    const Segment *e = &exons->items[strand == '+' ? order : exons->count - 1 - order];
    for (int64_t p = 0; p < e->end - e->start; p++)
      transcript[at++] = strand == '+' ? sequence->bases[e->start + p]
                                       : complement(sequence->bases[e->end - 1 - p]);
  }
  int64_t orf_start = INT64_MAX, orf_end = 0;
  int found = longest_orf(transcript, length, &orf_start, &orf_end);
  free(transcript);
  if (!found) return 0;
  SegmentList utr5 = {0}, cds = {0}, utr3 = {0};
  map_range(exons, strand, 0, orf_start, &utr5);
  map_range(exons, strand, orf_start, orf_end, &cds);
  map_range(exons, strand, orf_end, length, &utr3);
  char gene[48], mrna[64], id[96];
  snprintf(gene, sizeof(gene), "gene%lu", gene_number);
  snprintf(mrna, sizeof(mrna), "%s.t1", gene);
  Record *records = NULL;
  int count = 0;
  int64_t start = exons->items[0].start, end = exons->items[exons->count - 1].end;
  record_add(&records, &count, "gene", start, end, -1, gene, "");
  record_add(&records, &count, "mRNA", start, end, -1, mrna, gene);
  for (int i = 0; i < exons->count; i++) {
    snprintf(id, sizeof(id), "%s.exon%d", mrna, i + 1);
    record_add(&records, &count, "exon", exons->items[i].start, exons->items[i].end,
               -1, id, mrna);
  }
  for (int i = 0; i < utr5.count; i++) {
    snprintf(id, sizeof(id), "%s.utr5p%d", mrna, i + 1);
    record_add(&records, &count, "five_prime_UTR", utr5.items[i].start,
               utr5.items[i].end, -1, id, mrna);
  }
  for (int i = 0; i < cds.count; i++) {
    int64_t before = 0; /* coding bases upstream in transcript orientation */
    for (int j = 0; j < cds.count; j++)
      if (strand == '+' ? j < i : j > i) before += cds.items[j].end - cds.items[j].start;
    snprintf(id, sizeof(id), "%s.cds%d", mrna, i + 1);
    record_add(&records, &count, "CDS", cds.items[i].start, cds.items[i].end,
               (int)((3 - before % 3) % 3), id, mrna);
  }
  for (int i = 0; i < utr3.count; i++) {
    snprintf(id, sizeof(id), "%s.utr3p%d", mrna, i + 1);
    record_add(&records, &count, "three_prime_UTR", utr3.items[i].start,
               utr3.items[i].end, -1, id, mrna);
  }
  int ok = records_valid(records, count, exons);
  if (!ok) fail("internal error: inconsistent GFF3 records for %s", gene);
  for (int i = 0; ok && i < count; i++) {
    const Record *r = &records[i];
    char phase[12] = ".";
    if (r->phase >= 0) snprintf(phase, sizeof(phase), "%d", r->phase);
    ok = fprintf(output, "%s\tanno\t%s\t%lld\t%lld\t.\t%c\t%s\tID=%s%s%s\n",
                 sequence->name, r->type, (long long)r->start + 1,
                 (long long)r->end, strand, phase, r->id,
                 r->parent[0] ? ";Parent=" : "", r->parent) > 0;
  }
  free(records);
  free(utr5.items);
  free(cds.items);
  free(utr3.items);
  return ok ? 1 : -1;
}

/** @brief Strand of a class: '+', '-' or 0 for intergenic. */
static char class_strand(uint8_t cls) {
  return cls == ANNO_INTERGENIC ? 0 : cls <= 4 ? '+' : '-';
}

long anno_write_gff(FILE *output, const AnnoSequence *sequence,
                    const uint8_t *path, unsigned long *gene_counter) {
  long written = 0;
  for (int64_t p = 0; p < sequence->length;) {
    char strand = class_strand(path[p]);
    if (!strand) {
      p++;
      continue;
    }
    SegmentList exons = {0};
    while (p < sequence->length && class_strand(path[p]) == strand) {
      int intron = path[p] == ANNO_PLUS_INTRON || path[p] == ANNO_MINUS_INTRON;
      int64_t start = p;
      while (p < sequence->length && class_strand(path[p]) == strand &&
             (path[p] == ANNO_PLUS_INTRON || path[p] == ANNO_MINUS_INTRON) == intron)
        p++;
      if (!intron) segments_push(&exons, start, p);
    }
    int result = exons.count
      ? write_transcript(output, sequence, &exons, strand, *gene_counter + 1) : 0;
    free(exons.items);
    if (result < 0) return -1;
    if (result) {
      ++*gene_counter;
      written++;
    }
  }
  return written;
}

/* ------------------------------------------------------------------------ */
/* Prediction pipeline                                                      */
/* ------------------------------------------------------------------------ */

/** Shared context of one batch of inference windows. */
typedef struct {
  const Unet *unet;
  const AnnoSequence *sequence;
  Workspace **workspaces;
  int64_t first_window, windows, core;
  float *rows; /* windows x core x ANNO_CLASSES */
  int failed;
} PredictBatch;

/** @brief Job: predict one window, keep its core rows. */
static void predict_job(void *pointer, int64_t job, int worker) {
  PredictBatch *batch = pointer;
  Workspace *w = batch->workspaces[worker];
  const AnnoSequence *s = batch->sequence;
  int64_t core_start = (batch->first_window + job) * batch->core;
  int64_t core_end = core_start + batch->core < s->length ? core_start + batch->core : s->length;
  if (core_start >= core_end) return;
  int64_t start = core_start > ANNO_UNET_MARGIN ? core_start - ANNO_UNET_MARGIN : 0;
  int64_t end = core_end + ANNO_UNET_MARGIN < s->length ? core_end + ANNO_UNET_MARGIN : s->length;
  int64_t multiple = unet_multiple(&batch->unet->config);
  int64_t n = (end - start + multiple - 1) / multiple * multiple;
  float *power = xcalloc(mul_size((size_t)(end - start), ANNO_CWT_CHANNELS), sizeof(float));
  cwt_power(s->bases, s->length, start, end, power);
  memset(w->input, 0, (size_t)n * ANNO_INPUT_CHANNELS * sizeof(float));
  AnnoSequence view = {.bases = s->bases + start, .length = end - start};
  fill_input(batch->unet, &view, power, 0, end - start, 0, w->input);
  free(power);
  if (unet_forward(batch->unet, w, n)) {
    batch->failed = 1;
    return;
  }
  memcpy(batch->rows + (size_t)job * (size_t)batch->core * ANNO_CLASSES,
         w->probability + (size_t)(core_start - start) * ANNO_CLASSES,
         (size_t)(core_end - core_start) * ANNO_CLASSES * sizeof(float));
}

/** Run-length bedGraph writers, one file per class (values in 1/1000). */
typedef struct {
  FILE *files[ANNO_CLASSES];
  int value[ANNO_CLASSES];
  int64_t start[ANNO_CLASSES];
} Bedgraph;

/** @brief Write the pending run of class j ending at `position`. */
static void bedgraph_emit(Bedgraph *b, const char *name, int j, int64_t position) {
  if (b->value[j] >= 0 && position > b->start[j] &&
      fprintf(b->files[j], "%s\t%lld\t%lld\t%d.%02d\n", name, (long long)b->start[j],
              (long long)position, b->value[j] / 100, b->value[j] % 100) < 0)
    fail("cannot write bedGraph");
}

/** @brief Add one probability row at `position` (rows arrive in order). */
static void bedgraph_add(Bedgraph *b, const char *name, int64_t position,
                         const float *row) {
  if (!b->files[0]) return;
  for (int j = 0; j < ANNO_CLASSES; j++) {
    int value = (int)lrintf(row[j] * 100.0f); /* hundredths keep runs long */
    if (value == b->value[j]) continue;
    bedgraph_emit(b, name, j, position);
    b->value[j] = value;
    b->start[j] = position;
  }
}

/**
 * @brief Predict, decode and write GFF3 for every sequence.
 * @return number of genes written.
 */
static unsigned long predict_genome(const Unet *unet, const AnnoSequence *sequences,
                                    int count, int threads, FILE *gff,
                                    const char *prefix) {
  Bedgraph bedgraph = {0};
  for (int j = 0; prefix && j < ANNO_CLASSES; j++) {
    char path[4096];
    if (snprintf(path, sizeof(path), "%s.%s.bedgraph", prefix, class_names[j]) >=
        (int)sizeof(path))
      fail("output prefix too long");
    if (!(bedgraph.files[j] = fopen(path, "w"))) fail("cannot create %s", path);
    fprintf(bedgraph.files[j], "track type=bedGraph name=\"%s\"\n", class_names[j]);
  }
  int workers = threads < 1 ? 1 : threads;
  Workspace **workspaces = xcalloc((size_t)workers, sizeof(*workspaces));
  for (int i = 0; i < workers; i++)
    workspaces[i] = workspace_create(unet, ANNO_UNET_WINDOW, 0, 1);
  const int64_t core = ANNO_UNET_WINDOW - 2 * ANNO_UNET_MARGIN;
  float *rows = xcalloc(mul_size(mul_size((size_t)workers, (size_t)core), ANNO_CLASSES),
                        sizeof(float));
  unsigned long genes = 0;
  for (int i = 0; i < count; i++) {
    const AnnoSequence *s = &sequences[i];
    double started = now_seconds();
    Decoder *decoder = decoder_create(s->bases, s->length);
    int64_t windows = (s->length + core - 1) / core;
    for (int j = 0; j < ANNO_CLASSES; j++) bedgraph.value[j] = -1;
    for (int64_t first = 0; first < windows; first += workers) {
      int64_t batch_windows = windows - first < workers ? windows - first : workers;
      PredictBatch batch = {unet, s, workspaces, first, batch_windows, core, rows, 0};
      run_parallel(workers, batch_windows, predict_job, &batch);
      if (batch.failed) fail("non-finite UNet output for %s", s->name);
      int64_t begin = first * core;
      int64_t stop = (first + batch_windows) * core < s->length
        ? (first + batch_windows) * core : s->length;
      for (int64_t p = begin; p < stop; p++) {
        const float *row = rows + (size_t)(p - begin) * ANNO_CLASSES;
        bedgraph_add(&bedgraph, s->name, p, row);
        decoder_step(decoder, row);
      }
    }
    for (int j = 0; j < ANNO_CLASSES && bedgraph.files[0]; j++)
      bedgraph_emit(&bedgraph, s->name, j, s->length);
    uint8_t *path = xcalloc((size_t)s->length, 1);
    if (decoder_finish(decoder, path)) fail("grammar decoding failed for %s", s->name);
    decoder_destroy(decoder);
    long written = anno_write_gff(gff, s, path, &genes);
    free(path);
    if (written < 0) fail("cannot write GFF3");
    say("%s: %ld genes (%.1fs)", s->name, written, now_seconds() - started);
  }
  for (int j = 0; j < ANNO_CLASSES; j++)
    if (bedgraph.files[j] && fclose(bedgraph.files[j])) fail("cannot close bedGraph");
  for (int i = 0; i < workers; i++) workspace_destroy(workspaces[i]);
  free(workspaces);
  free(rows);
  return genes;
}

/* ------------------------------------------------------------------------ */
/* Command line                                                             */
/* ------------------------------------------------------------------------ */

#ifndef ANNO_NO_MAIN
/** @brief Print usage to stderr. */
static void usage(void) {
  fprintf(stderr,
    "Usage: anno [options] genome.fasta[.gz] [annotation.gff3[.gz]]\n"
    "       anno [-p N] cwt genome.fasta[.gz] directory\n"
    "  with an annotation: train a UNet and write it to -m\n"
    "  without:            predict genes (GFF3 on stdout unless -o)\n"
    "Options:\n"
    "  -p N       threads [1]\n"
    "  -i N       training iterations [%d]\n"
    "  -m FILE    model file [anno.ann]\n"
    "  -o PREFIX  write PREFIX.gff3 and PREFIX.{class}.bedgraph\n"
    "  -h         this help\n", ANNO_DEFAULT_ITERATIONS);
}

/** @brief Parse a positive int option value. */
static int positive(const char *text, const char *name) {
  char *end = NULL;
  errno = 0;
  long value = strtol(text, &end, 10);
  if (errno || end == text || *end || value < 1 || value > INT_MAX)
    fail("invalid %s: %s", name, text);
  return (int)value;
}

/** @brief Entry point. */
int main(int argc, char **argv) {
  int threads = 1, iterations = ANNO_DEFAULT_ITERATIONS;
  const char *model_path = "anno.ann", *prefix = NULL;
  int a = 1;
  for (; a < argc && argv[a][0] == '-' && argv[a][1]; a++) {
    const char *option = argv[a];
    if (!strcmp(option, "-h")) {
      usage();
      return EXIT_SUCCESS;
    }
    if (strlen(option) != 2 || !strchr("pimo", option[1]) || a + 1 >= argc) {
      usage();
      return EXIT_FAILURE;
    }
    const char *value = argv[++a];
    if (option[1] == 'p') threads = positive(value, "thread count");
    if (option[1] == 'i') iterations = positive(value, "iteration count");
    if (option[1] == 'm') model_path = value;
    if (option[1] == 'o') prefix = value;
  }
  if (a < argc && !strcmp(argv[a], "cwt")) {
    if (argc - a != 3) {
      usage();
      return EXIT_FAILURE;
    }
    run_cwt(argv[a + 1], argv[a + 2], threads);
    return EXIT_SUCCESS;
  }
  if (argc - a < 1 || argc - a > 2) {
    usage();
    return EXIT_FAILURE;
  }
  int count = 0;
  AnnoSequence *sequences = anno_read_fasta(argv[a], &count);
  if (argc - a == 2) {
    double started = now_seconds();
    anno_read_labels(argv[a + 1], sequences, count);
    Unet *unet = train_unet(sequences, count, iterations, threads);
    int status = unet_save(model_path, unet);
    if (status) fail("cannot write model %s: %s", model_path, strerror(status));
    say("trained in %.1fs; model written to %s", now_seconds() - started, model_path);
    unet_destroy(unet);
  } else {
    Unet *unet = NULL;
    int status = unet_load(model_path, &unet);
    if (status) fail("cannot read model %s: %s", model_path, strerror(status));
    FILE *gff = stdout;
    if (prefix) {
      char path[4096];
      if (snprintf(path, sizeof(path), "%s.gff3", prefix) >= (int)sizeof(path))
        fail("output prefix too long");
      if (!(gff = fopen(path, "w"))) fail("cannot create %s", path);
    }
    fputs("##gff-version 3\n", gff);
    unsigned long genes = predict_genome(unet, sequences, count, threads, gff, prefix);
    if (fflush(gff) || ferror(gff) || (gff != stdout && fclose(gff)))
      fail("cannot write GFF3");
    say("predicted %lu genes", genes);
    unet_destroy(unet);
  }
  anno_sequences_free(sequences, count);
  return EXIT_SUCCESS;
}
#endif
