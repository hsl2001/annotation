/**
 * @file anno.c
 * @brief anno: FASTA/GFF3 input, CWT, KANN 1D convolution model,
 *        9-state grammar Viterbi, longest ORF and GFF3 output.
 *
 * Pipeline: bases -> one-hot + CWT power -> UNet -> per-base probability of
 * ANNO_CLASSES classes -> log emissions -> grammar Viterbi -> transcripts ->
 * longest ORF -> gene/mRNA/exon/UTR/CDS records.
 */
#define _POSIX_C_SOURCE 200809L

#include "anno.h"
#include "klib/kann.h"
#include "klib/kthread.h"
#include "klib/ketopt.h"

#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>
#include <zlib.h>

#include "klib/kseq.h"
KSEQ_INIT(gzFile, gzread)

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

/**
 * @brief Print an error to stderr and terminate the process.
 * @param format printf-style error format string.
 * @param ... Values referenced by @p format.
 * @return This function does not return.
 */
static _Noreturn void fail(const char *format, ...) {
  va_list arguments;
  va_start(arguments, format);
  fputs("anno: error: ", stderr);
  vfprintf(stderr, format, arguments);
  fputc('\n', stderr);
  va_end(arguments);
  exit(EXIT_FAILURE);
}

/**
 * @brief Print a log line to stderr.
 * @param format printf-style log format string.
 * @param ... Values referenced by @p format.
 * @return None.
 */
static void say(const char *format, ...) {
  va_list arguments;
  va_start(arguments, format);
  fputs("anno: ", stderr);
  vfprintf(stderr, format, arguments);
  fputc('\n', stderr);
  va_end(arguments);
}

/**
 * @brief Multiply two sizes, failing on overflow.
 * @param a First factor.
 * @param b Second factor.
 * @return Product of @p a and @p b.
 */
static size_t mul_size(size_t a, size_t b) {
  if (a && b > SIZE_MAX / a) fail("size overflow");
  return a * b;
}

/**
 * @brief Allocate zero-initialized memory or terminate on failure.
 * @param count Number of elements.
 * @param size Size of each element.
 * @return Allocated zero-initialized memory.
 */
static void *xcalloc(size_t count, size_t size) {
  void *memory = calloc(count ? count : 1, size ? size : 1);
  if (!memory || (size && count > SIZE_MAX / size)) fail("out of memory");
  return memory;
}

/**
 * @brief Resize an allocation or terminate on failure.
 * @param old Existing allocation, or NULL.
 * @param count Number of elements.
 * @param size Size of each element.
 * @return Resized allocation.
 */
static void *xrealloc(void *old, size_t count, size_t size) {
  void *memory = realloc(old, mul_size(count ? count : 1, size));
  if (!memory) fail("out of memory");
  return memory;
}

/**
 * @brief Duplicate a string or terminate on allocation failure.
 * @param text String to duplicate.
 * @return Newly allocated string.
 */
static char *xstrdup(const char *text) {
  char *copy = strdup(text);
  if (!copy) fail("out of memory");
  return copy;
}

/**
 * @brief Return monotonic elapsed time in seconds.
 * @return Monotonic clock value in seconds.
 */
static double now_seconds(void) {
  struct timespec time;
  clock_gettime(CLOCK_MONOTONIC, &time);
  return (double)time.tv_sec + 1e-9 * (double)time.tv_nsec;
}

/**
 * @brief Convert a nucleotide to its one-hot channel index.
 * @param base Nucleotide character.
 * @return A/C/G/T index in the range 0..3, or -1 for another character.
 */
static int base_index(char base) {
  switch (base) {
  case 'A': return 0;
  case 'C': return 1;
  case 'G': return 2;
  case 'T': return 3;
  default: return -1;
  }
}

/**
 * @brief Return the complement of a nucleotide.
 * @param base Nucleotide character.
 * @return Complementary nucleotide, or N for an unknown base.
 */
static char complement(char base) {
  int index = base_index(base);
  return index < 0 ? 'N' : "TGCA"[index];
}

/* ------------------------------------------------------------------------ */
/* Threads                                                                  */
/* ------------------------------------------------------------------------ */

typedef void (*JobFunction)(void *context, long job, int worker);

/**
 * @brief Run independent jobs through klib's static worker pool.
 * @param threads Number of worker threads.
 * @param jobs Number of jobs.
 * @param function Job callback.
 * @param context User context passed to @p function.
 * @return None.
 */
static void run_parallel(int threads, int64_t jobs, JobFunction function,
                         void *context) {
  kt_for(threads, function, context, (long)jobs);
}

/* ------------------------------------------------------------------------ */
/* FASTA                                                                    */
/* ------------------------------------------------------------------------ */

/**
 * @brief Read FASTA records from a possibly compressed file.
 * @param path FASTA file path.
 * @param count_out Receives the number of records.
 * @return Caller-owned sequence array.
 */
AnnoSequence *anno_read_fasta(const char *path, int *count_out) {
  gzFile input = gzopen(path, "rb");
  if (!input) fail("cannot open FASTA %s", path);
  AnnoSequence *sequences = NULL;
  int count = 0;
  kseq_t *reader = kseq_init(input);
  int length;
  while ((length = kseq_read(reader)) >= 0) {
    if (!reader->name.l) fail("empty FASTA sequence ID in %s", path);
    if (!length) fail("empty FASTA sequence %s", reader->name.s);
    for (int i = 0; i < count; i++)
      if (!strcmp(sequences[i].name, reader->name.s))
        fail("duplicate FASTA sequence ID %s", reader->name.s);
    sequences = xrealloc(sequences, (size_t)count + 1, sizeof(*sequences));
    AnnoSequence *sequence = &sequences[count++];
    *sequence = (AnnoSequence){.name = xstrdup(reader->name.s), .length = length,
                               .bases = xcalloc((size_t)length + 1, 1)};
    for (int i = 0; i < length; i++)
      sequence->bases[i] = (char)toupper((unsigned char)reader->seq.s[i]);
  }
  int error = length < -1;
  kseq_destroy(reader);
  if (gzclose(input) != Z_OK || error) fail("cannot read FASTA %s", path);
  if (!count) fail("no FASTA records in %s", path);
  *count_out = count;
  return sequences;
}

/**
 * @brief Release an array of sequence records and their owned data.
 * @param sequences Sequence records, or NULL.
 * @param count Number of records.
 * @return None.
 */
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

/**
 * @brief Append a half-open segment to a growable list.
 * @param list Segment list to extend.
 * @param start Inclusive segment start.
 * @param end Exclusive segment end.
 * @return None.
 */
static void segments_push(SegmentList *list, int64_t start, int64_t end) {
  if (list->count == list->capacity) {
    list->capacity = list->capacity ? list->capacity * 2 : 8;
    list->items = xrealloc(list->items, (size_t)list->capacity,
                           sizeof(*list->items));
  }
  list->items[list->count++] = (Segment){start, end};
}

/**
 * @brief Compare segments by start and then end.
 * @param a First segment.
 * @param b Second segment.
 * @return Negative, zero, or positive comparison result.
 */
static int segment_compare(const void *a, const void *b) {
  const Segment *x = a, *y = b;
  return x->start < y->start ? -1 : x->start > y->start ? 1
       : x->end < y->end ? -1 : x->end > y->end;
}

/**
 * @brief Compare GFF features by parent, kind, and start.
 * @param a First feature.
 * @param b Second feature.
 * @return Negative, zero, or positive comparison result.
 */
static int feature_compare(const void *a, const void *b) {
  const GffFeature *x = a, *y = b;
  int order = strcmp(x->parent, y->parent);
  if (order) return order;
  if (x->kind != y->kind) return x->kind - y->kind;
  return x->start < y->start ? -1 : x->start > y->start;
}

/**
 * @brief Convert a GFF feature type to an internal kind.
 * @param type GFF feature type.
 * @return KIND_* value, -2 for container features, or -1 when ignored.
 */
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

/**
 * @brief Parse a non-negative decimal integer.
 * @param text Decimal text to parse.
 * @param value Receives the parsed integer.
 * @return 1 on success, otherwise 0.
 */
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

/**
 * @brief Compare transcripts by decreasing genomic span.
 * @param a First transcript.
 * @param b Second transcript.
 * @return Negative, zero, or positive comparison result.
 */
static int transcript_compare(const void *a, const void *b) {
  const LabelTranscript *x = a, *y = b;
  int64_t lx = x->end - x->start, ly = y->end - y->start;
  if (lx != ly) return lx > ly ? -1 : 1;
  if (x->sequence != y->sequence) return x->sequence - y->sequence;
  return x->start < y->start ? -1 : x->start > y->start;
}

/**
 * @brief Merge touching or overlapping segments in place.
 * @param list Segment list to sort and merge.
 * @return None.
 */
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

/**
 * @brief Read GFF features that belong to transcripts.
 * @param path GFF3 file path.
 * @param sequences FASTA sequences used to resolve sequence IDs.
 * @param count Number of FASTA sequences.
 * @param feature_count Receives the number of features.
 * @return Caller-owned feature array.
 */
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

/**
 * @brief Paint one transcript into a per-base label array.
 * @param t Transcript to paint.
 * @param labels Destination labels; overlapping transcripts are prechecked.
 * @return None.
 */
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

/**
 * @brief Read GFF3 annotations and fill sequence label arrays.
 * @param path GFF3 file path.
 * @param sequences FASTA sequence records to label.
 * @param count Number of sequence records.
 * @return None.
 */
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

/**
 * @brief Compute the complex Morlet CWT for a sequence interval.
 * @param bases Input sequence bases.
 * @param length Total sequence length.
 * @param start Inclusive interval start.
 * @param end Exclusive interval end.
 * @param out Destination complex CWT values.
 * @return None.
 */
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

/**
 * @brief Compute CWT power for a sequence interval.
 * @param bases Input sequence bases.
 * @param length Total sequence length.
 * @param start Inclusive interval start.
 * @param end Exclusive interval end.
 * @param out Destination power matrix.
 * @return None.
 */
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

/** Context of the cwt-subcommand block job. */
typedef struct {
  const AnnoSequence *sequence;
  int64_t start;
  float *out;
} ComplexJob;

/**
 * @brief Compute complex CWT for one cwt-subcommand block.
 * @param pointer ComplexJob context.
 * @param job Window number.
 * @param worker Worker number, unused.
 * @return None.
 */
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
 * @brief Write CWT matrices and contig offsets for the cwt subcommand.
 *        reverse-complement rows) and {dir}/contigs.tsv for plot.py.
 *        Minus rows are the conjugated, mirrored plus rows (no recompute).
 * @param fasta Input FASTA file.
 * @param directory Output directory.
 * @param threads Number of worker threads.
 * @return None.
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


#define KANN_LENGTH ANNO_UNET_WINDOW

/** @brief Opaque public model backed by a KANN graph. */
struct Unet {
  UnetConfig config;
  kann_t *ann;
  int output;
};

/**
 * @brief Build the fixed-length KANN convolution graph and Soft F1 cost.
 * @param config Convolution configuration.
 * @return Newly allocated KANN graph, or NULL on graph construction failure.
 */
static kann_t *kann_unet(const UnetConfig *config) {
  kad_node_t *input = kad_feed(3, ANNO_BATCH_WINDOWS, ANNO_INPUT_CHANNELS, KANN_LENGTH);
  input->ext_flag |= KANN_F_IN;
  kad_node_t *skip[ANNO_MAX_DEPTH] = {0};
  kad_node_t *x = input;
  for (int level = 0; level < config->depth; level++) {
    x = kad_relu(kann_layer_conv1d(x, config->channels[level], config->kernel,
                                   1, KAD_PAD_SAME));
    x = kad_relu(kann_layer_conv1d(x, config->channels[level], config->kernel,
                                   1, KAD_PAD_SAME));
    skip[level] = x;
    if (level + 1 < config->depth)
      x = kad_avg1d(x, config->pool, config->pool, 0);
  }
  for (int level = config->depth - 2; level >= 0; level--) {
    x = kad_upsample1d(x, config->pool);
    x = kad_concat(1, 2, x, skip[level]);
    x = kad_relu(kann_layer_conv1d(x, config->channels[level], config->kernel,
                                   1, KAD_PAD_SAME));
    x = kad_relu(kann_layer_conv1d(x, config->channels[level], config->kernel,
                                   1, KAD_PAD_SAME));
  }
  x = kann_layer_conv1d(x, ANNO_CLASSES, 1, 1, 0);
  x->ext_flag |= KANN_F_OUT;
  kad_node_t *truth = kad_feed(3, ANNO_BATCH_WINDOWS, ANNO_CLASSES, KANN_LENGTH);
  truth->ext_flag |= KANN_F_TRUTH;
  truth->ext_label = 2;
  kad_node_t *mask = kad_feed(3, ANNO_BATCH_WINDOWS, 1, KANN_LENGTH);
  mask->ext_flag |= KANN_F_TRUTH;
  mask->ext_label = 1;
  kad_node_t *minus_one = kann_new_scalar(KAD_CONST, -1.0f);
  kad_node_t *epsilon = kann_new_scalar(KAD_CONST, ANNO_EPSILON);
  kad_node_t *probabilities = kad_softmax_axis(x, 1);
  kad_node_t *loss[ANNO_CLASSES];
  kad_node_t *cross_entropy[ANNO_CLASSES];
  if (!probabilities || !epsilon)
    fail("invalid Soft F1 probability graph");
  for (int j = 0; j < ANNO_CLASSES; j++) {
    kad_node_t *p = kad_slice(probabilities, 1, j, j + 1);
    kad_node_t *y = kad_slice(truth, 1, j, j + 1);
    kad_node_t *masked_p = kad_mul(p, mask);
    kad_node_t *masked_y = kad_mul(y, mask);
    kad_node_t *intersection = kad_reduce_sum(kad_mul(masked_p, masked_y), 2);
    kad_node_t *sum_p = kad_reduce_sum(masked_p, 2);
    kad_node_t *sum_y = kad_reduce_sum(masked_y, 2);
    kad_node_t *numerator = kad_add(kad_mul(intersection,
                                             kann_new_scalar(KAD_CONST, 2.0f)),
                                    epsilon);
    kad_node_t *denominator = kad_add(kad_add(sum_p, sum_y), epsilon);
    kad_node_t *dice = kad_mul(numerator,
                               kad_exp(kad_mul(kad_log(denominator),
                                               minus_one)));
    kad_node_t *class_loss = kad_1minus(dice);
    cross_entropy[j] = kad_ce_bin(masked_p, masked_y);
    if (!p || !y || !masked_p || !masked_y || !intersection || !sum_p || !sum_y ||
        !cross_entropy[j] || !numerator ||
        !denominator || !dice || !class_loss)
      fail("invalid Soft F1 loss graph");
    class_loss = kad_reduce_mean(class_loss, 0);
    loss[j] = kad_reduce_mean(class_loss, 0);
    if (!loss[j]) fail("invalid Soft F1 reduction");
  }
  kad_node_t *dice_cost = kad_avg(ANNO_CLASSES, loss);
  kad_node_t *ce_cost = kad_avg(ANNO_CLASSES, cross_entropy);
  kad_node_t *cost = kad_add(dice_cost,
                             kad_mul(kann_new_scalar(KAD_CONST, 0.2f), ce_cost));
  cost->ext_flag |= KANN_F_COST;
  kann_t *ann = kann_new(cost, 0);
  kann_set_batch_size(ann, 1);
  return ann;
}

/**
 * @brief Map a class label to its reverse-complement strand.
 * @param label Class label on either strand.
 * @return Equivalent class on the opposite strand.
 */
static uint8_t kann_swap_strand(uint8_t label) {
  return label == ANNO_INTERGENIC ? label : label <= 4 ? label + 4 : label - 4;
}

/**
 * @brief Return the default KANN model configuration.
 * @return Default model configuration.
 */
UnetConfig unet_default_config(void) {
  return (UnetConfig){.depth = 4, .kernel = 9, .pool = 2,
                      .channels = {16, 32, 64, 128}};
}

/**
 * @brief Create a KANN model.
 * @param config Model configuration.
 * @return New model, or NULL when @p config is NULL or allocation fails.
 */
Unet *unet_create(const UnetConfig *config) {
  if (!config) return NULL;
  Unet *unet = xcalloc(1, sizeof(*unet));
  unet->config = *config;
  unet->ann = kann_unet(&unet->config);
  unet->output = kann_find(unet->ann, KANN_F_OUT, 0);
  if (!unet->output) {
    unet_destroy(unet);
    return NULL;
  }
  return unet;
}

/**
 * @brief Release a KANN model and its graph.
 * @param unet Model to release, or NULL.
 * @return None.
 */
void unet_destroy(Unet *unet) {
  if (!unet) return;
  kann_delete(unet->ann);
  free(unet);
}

/**
 * @brief Fill channel-major KANN input and one-hot target tensors.
 * @param sequence Source sequence.
 * @param power Precomputed CWT power channels.
 * @param power_start First genomic position represented by @p power.
 * @param start First source position.
 * @param valid Number of positions to fill.
 * @param flip Whether to reverse-complement the window.
 * @param input Channel-major input tensor.
 * @param labels Optional source labels.
 * @param truth Optional channel-major one-hot target tensor.
 * @return None.
 */
static void features(const AnnoSequence *sequence, const float *power,
                     int64_t power_start, int64_t start, int64_t valid,
                     int flip, float *input, const uint8_t *labels,
                     float *truth) {
  float mean[ANNO_CWT_CHANNELS] = {0};
  float variance[ANNO_CWT_CHANNELS] = {0};
  for (int64_t q = 0; q < valid; q++)
    for (int s = 0; s < ANNO_CWT_CHANNELS; s++)
      mean[s] += power[(size_t)(start + q - power_start) *
                       ANNO_CWT_CHANNELS + s];
  for (int s = 0; s < ANNO_CWT_CHANNELS; s++) mean[s] /= (float)valid;
  for (int64_t q = 0; q < valid; q++)
    for (int s = 0; s < ANNO_CWT_CHANNELS; s++) {
      float delta = power[(size_t)(start + q - power_start) *
                         ANNO_CWT_CHANNELS + s] - mean[s];
      variance[s] += delta * delta;
    }
  for (int s = 0; s < ANNO_CWT_CHANNELS; s++)
    variance[s] = sqrtf(variance[s] / (float)valid + 1e-6f);
  for (int64_t q = 0; q < valid; q++) {
    int64_t p = flip ? start + valid - 1 - q : start + q;
    int base = base_index(sequence->bases[p]);
    if (flip && base >= 0) base = 3 - base;
    if (base >= 0) input[(size_t)base * KANN_LENGTH + q] = 1.0f;
    for (int s = 0; s < ANNO_CWT_CHANNELS; s++)
      input[(size_t)(4 + s) * KANN_LENGTH + q] =
        (power[(size_t)(p - power_start) * ANNO_CWT_CHANNELS + s] - mean[s]) /
        variance[s];
    if (truth) {
      uint8_t label = labels[flip ? start + valid - 1 - q : start + q];
      if (flip) label = kann_swap_strand(label);
      truth[(size_t)label * KANN_LENGTH + q] = 1.0f;
    }
  }
}

/**
 * @brief Evaluate a sequence and return per-base class probabilities.
 * @param unet Model to evaluate.
 * @param input Position-major input features.
 * @param length Number of positions.
 * @param probability Position-major output probabilities.
 * @return 0 on success, or an errno value otherwise.
 */
static int unet_predict_channels(const Unet *unet, const float *input,
                                  int64_t length, float *probability) {
  if (!unet || !input || !probability || length < 1 || length > KANN_LENGTH)
    return EINVAL;
  float *bound = (float *)input;
  kann_feed_bind(unet->ann, KANN_F_IN, 0, &bound);
  kann_set_batch_size(unet->ann, 1);
  kann_eval_out(unet->ann);
  const float *raw = unet->ann->v[unet->output]->x;
  for (int64_t p = 0; p < length; p++) {
    float top = raw[p], sum = 0.0f;
    for (int j = 1; j < ANNO_CLASSES; j++) if (raw[(size_t)j * KANN_LENGTH + p] > top)
      top = raw[(size_t)j * KANN_LENGTH + p];
    for (int j = 0; j < ANNO_CLASSES; j++) {
      float value = expf(raw[(size_t)j * KANN_LENGTH + p] - top);
      probability[(size_t)p * ANNO_CLASSES + j] = value;
      sum += value;
    }
    for (int j = 0; j < ANNO_CLASSES; j++) probability[(size_t)p * ANNO_CLASSES + j] /= sum;
  }
  return 0;
}

int unet_predict(const Unet *unet, const float *input, int64_t length,
                 float *probability) {
  if (!unet || !input || !probability || length < 1 || length > KANN_LENGTH)
    return EINVAL;
  float *channel_input = xcalloc((size_t)ANNO_INPUT_CHANNELS * KANN_LENGTH,
                                 sizeof(float));
  for (int64_t p = 0; p < length; p++)
    for (int c = 0; c < ANNO_INPUT_CHANNELS; c++)
      channel_input[(size_t)c * KANN_LENGTH + p] =
        input[(size_t)p * ANNO_INPUT_CHANNELS + c];
  int status = unet_predict_channels(unet, channel_input, length, probability);
  free(channel_input);
  return status;
}

/**
 * @brief Save a KANN graph and its parameters.
 * @param path Destination model path.
 * @param unet Model to save.
 * @return 0 on success, or EINVAL for invalid arguments.
 */
int unet_save(const char *path, const Unet *unet) {
  if (!path || !unet) return EINVAL;
  kann_save(path, unet->ann);
  return 0;
}

/**
 * @brief Load a KANN graph and wrap it as an Unet handle.
 * @param path Source model path.
 * @param out Receives the loaded model.
 * @return 0 on success, or EINVAL when loading fails.
 */
int unet_load(const char *path, Unet **out) {
  if (!path || !out) return EINVAL;
  *out = NULL;
  kann_t *ann = kann_load(path);
  if (!ann) return EINVAL;
  Unet *unet = xcalloc(1, sizeof(*unet));
  unet->config = unet_default_config();
  unet->ann = ann;
  unet->output = kann_find(ann, KANN_F_OUT, 0);
  if (unet->output < 0) {
    unet_destroy(unet);
    return EINVAL;
  }
  *out = unet;
  return 0;
}

typedef struct {
  const AnnoSequence *sequence;
  int64_t start;
  int64_t valid;
} TrainingWindow;

static uint64_t training_random_next(uint64_t *state) {
  uint64_t value = *state;
  value ^= value >> 12;
  value ^= value << 25;
  value ^= value >> 27;
  *state = value;
  return value * UINT64_C(2685821657736338717);
}

typedef struct {
  const TrainingWindow *windows;
  size_t batch_start;
  float *input;
  float *truth;
  float *mask;
  int flip;
} TrainingBatchJob;

static void training_batch_job(void *pointer, int64_t job, int worker) {
  (void)worker;
  TrainingBatchJob *context = pointer;
  const TrainingWindow *window =
    &context->windows[context->batch_start + (size_t)job];
  size_t slot = (size_t)job;
  float *power = xcalloc(mul_size((size_t)window->valid,
                                  ANNO_CWT_CHANNELS), sizeof(float));
  cwt_power(window->sequence->bases, window->sequence->length,
            window->start, window->start + window->valid, power);
  features(window->sequence, power, window->start, window->start,
           window->valid, context->flip,
           context->input + slot * ANNO_INPUT_CHANNELS * KANN_LENGTH,
           window->sequence->labels,
           context->truth + slot * ANNO_CLASSES * KANN_LENGTH);
  for (int64_t q = 0; q < window->valid; q++)
    context->mask[slot * KANN_LENGTH + q] = 1.0f;
  free(power);
}

/**
 * @brief Train the KANN graph on every labelled sequence window.
 * @param sequences Labelled training sequences.
 * @param count Number of sequences.
 * @param iterations Number of optimizer iterations.
 * @param threads Number of KANN worker threads.
 * @return Newly trained model.
 */
static Unet *train_unet(const AnnoSequence *sequences, int count, int iterations,
                        int threads) {
  UnetConfig config = unet_default_config();
  Unet *unet = unet_create(&config);
  if (!unet) fail("cannot create KANN network");
  kann_mt(unet->ann, threads, ANNO_BATCH_WINDOWS);
  int nvar = kann_size_var(unet->ann);
  float *x = xcalloc((size_t)ANNO_BATCH_WINDOWS * ANNO_INPUT_CHANNELS * KANN_LENGTH, sizeof(float));
  float *y = xcalloc((size_t)ANNO_BATCH_WINDOWS * ANNO_CLASSES * KANN_LENGTH, sizeof(float));
  float *mask = xcalloc((size_t)ANNO_BATCH_WINDOWS * KANN_LENGTH, sizeof(float));
  float *r = xcalloc((size_t)nvar, sizeof(float));
  size_t total_windows = 0;
  for (int i = 0; i < count; i++)
    total_windows += (size_t)((sequences[i].length + KANN_LENGTH - 1) /
                              KANN_LENGTH);
  if (!total_windows) fail("no training windows");
  TrainingWindow *all_windows = xcalloc(total_windows, sizeof(*all_windows));
  TrainingWindow *batch_windows =
    xcalloc(ANNO_BATCH_WINDOWS, sizeof(*batch_windows));
  size_t all_index = 0;
  for (int sequence = 0; sequence < count; sequence++) {
    int64_t sequence_windows =
      (sequences[sequence].length + KANN_LENGTH - 1) / KANN_LENGTH;
    for (int64_t window = 0; window < sequence_windows; window++) {
      int64_t start = window * KANN_LENGTH;
      int64_t valid = sequences[sequence].length - start;
      if (valid > KANN_LENGTH) valid = KANN_LENGTH;
      all_windows[all_index++] =
        (TrainingWindow){&sequences[sequence], start, valid};
    }
  }
  uint64_t shuffle_state = UINT64_C(0x9e3779b97f4a7c15);
  double started = now_seconds();
  for (int iteration = 0; iteration < iterations; iteration++) {
    for (size_t i = total_windows; i > 1; i--) {
      size_t j = (size_t)(training_random_next(&shuffle_state) % i);
      TrainingWindow swap = all_windows[i - 1];
      all_windows[i - 1] = all_windows[j];
      all_windows[j] = swap;
    }
    double epoch_loss = 0.0;
    size_t epoch_batches = 0;
    for (size_t batch_start = 0; batch_start < total_windows;
         batch_start += ANNO_BATCH_WINDOWS) {
      size_t batch_end = batch_start + ANNO_BATCH_WINDOWS;
      if (batch_end > total_windows) batch_end = total_windows;
      memset(x, 0, (size_t)ANNO_BATCH_WINDOWS * ANNO_INPUT_CHANNELS *
             KANN_LENGTH * sizeof(float));
      memset(y, 0, (size_t)ANNO_BATCH_WINDOWS * ANNO_CLASSES *
             KANN_LENGTH * sizeof(float));
      memset(mask, 0, (size_t)ANNO_BATCH_WINDOWS * KANN_LENGTH * sizeof(float));
      size_t batch_size = batch_end - batch_start;
      memcpy(batch_windows, all_windows + batch_start,
             batch_size * sizeof(*batch_windows));
      TrainingBatchJob batch = {
        batch_windows, 0, x, y, mask,
        (iteration + (int)(batch_start / ANNO_BATCH_WINDOWS)) & 1
      };
      run_parallel(threads, (int64_t)(batch_end - batch_start),
                   training_batch_job, &batch);
      kann_feed_bind(unet->ann, KANN_F_IN, 0, &x);
      kann_feed_bind(unet->ann, KANN_F_TRUTH, 2, &y);
      kann_feed_bind(unet->ann, KANN_F_TRUTH, 1, &mask);
      kann_set_batch_size(unet->ann, (int)(batch_end - batch_start));
      float loss = kann_cost(unet->ann, 0, 1);
      kann_RMSprop(nvar, ANNO_LEARNING_RATE, NULL, 0.9f, unet->ann->g,
                   unet->ann->x, r);
      if (!isfinite(loss)) fail("non-finite KANN loss");
      epoch_loss += loss;
      epoch_batches++;
    }
    if (iteration == 0 || iteration + 1 == iterations ||
        (iteration + 1) % 10 == 0)
      say("iteration %d/%d loss %.4f batches %zu %.1fs",
          iteration + 1, iterations, epoch_loss / (double)epoch_batches,
          epoch_batches, now_seconds() - started);
  }
  kann_mt(unet->ann, 0, 0);
  free(all_windows);
  free(batch_windows);
  free(x); free(y); free(mask); free(r);
  return unet;
}

#define GRAMMAR_MAX_STATES 256
#define GRAMMAR_MAX_PREDECESSORS 5

/** @brief Pending class-wise bedGraph runs. */
typedef struct {
  int count;
  uint8_t cls[GRAMMAR_MAX_STATES];
  int predecessors[GRAMMAR_MAX_STATES];
  int from[GRAMMAR_MAX_STATES][GRAMMAR_MAX_PREDECESSORS];
  double score[GRAMMAR_MAX_STATES][GRAMMAR_MAX_PREDECESSORS];
  int shift[GRAMMAR_MAX_STATES], bits[GRAMMAR_MAX_STATES];
} Grammar;

/**
 * @brief Add a scored predecessor edge to a grammar state.
 * @param g Grammar to modify.
 * @param from Predecessor state index.
 * @param to Destination state index.
 * @param score Transition score.
 * @return None.
 */
static void grammar_edge(Grammar *g, int from, int to, double score) {
  int k = g->predecessors[to]++;
  if (k >= GRAMMAR_MAX_PREDECESSORS) fail("grammar has too many predecessors");
  g->from[to][k] = from;
  g->score[to][k] = score;
}

/**
 * @brief Add a self-looping chain of grammar states.
 * @param g Grammar to modify.
 * @param cls Class assigned to each new state.
 * @param length Number of states in the chain.
 * @return Index of the first new state.
 */
static int grammar_chain(Grammar *g, uint8_t cls, int length) {
  int head = g->count;
  for (int i = 0; i < length; i++) {
    g->cls[g->count] = cls;
    if (i) grammar_edge(g, g->count - 1, g->count, 0.0);
    g->count++;
  }
  grammar_edge(g, g->count - 1, g->count - 1, 0.0);
  return head;
}

/**
 * @brief Build the strand-aware gene grammar. Plus strand in genomic order:
 *   intergenic -> [5'UTR (intron 5'UTR)*] -> exon (intron exon)* ->
 *   [3'UTR (intron 3'UTR)*] -> intergenic; minus strand is the mirror image.
 *   Introns are never entered from intergenic and UTRs never touch directly.
 * @param g Grammar to initialize.
 * @return None.
 */
static void grammar_build(Grammar *g) {
  memset(g, 0, sizeof(*g));
  int inter = grammar_chain(g, ANNO_INTERGENIC, 1);
  for (int strand = 0; strand < 2; strand++) {
    int o = strand ? 4 : 0;
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
    grammar_edge(g, inter, first, -ANNO_TRANSCRIPT_PENALTY);
    grammar_edge(g, inter, ex, -ANNO_TRANSCRIPT_PENALTY);
    grammar_edge(g, first, ex, 0.0);
    grammar_edge(g, ex_tail, last, 0.0);
    grammar_edge(g, ex_tail, inter, 0.0);
    grammar_edge(g, last, inter, 0.0);
    int pieces[3] = {first, ex, last}, exits[3] = {first, ex_tail, last};
    for (int i = 0; i < 3; i++) {
      grammar_edge(g, exits[i], introns[i], 0.0);
      grammar_edge(g, introns[i] + tail, pieces[i], 0.0);
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

/** Streaming Viterbi decoder. Owns its buffers. */
typedef struct {
  Grammar grammar;
  int64_t length, position;
  double *previous, *current;
  uint32_t *back;
} Decoder;

/**
 * @brief Create a decoder for a sequence of the given length.
 * @param length Number of sequence positions.
 * @return Newly allocated decoder.
 */
static Decoder *decoder_create(int64_t length) {
  Decoder *d = xcalloc(1, sizeof(*d));
  grammar_build(&d->grammar);
  d->length = length;
  d->previous = xcalloc((size_t)d->grammar.count, sizeof(double));
  d->current = xcalloc((size_t)d->grammar.count, sizeof(double));
  d->back = xcalloc((size_t)length, sizeof(uint32_t));
  return d;
}

/**
 * @brief Free a decoder and its buffers.
 * @param d Decoder to free, or NULL.
 * @return None.
 */
static void decoder_destroy(Decoder *d) {
  if (!d) return;
  free(d->previous);
  free(d->current);
  free(d->back);
  free(d);
}

/**
 * @brief Consume the probability row of the next position.
 * @param d Decoder to advance.
 * @param row Class probabilities for one position.
 * @return None.
 */
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

/**
 * @brief Trace back the best path ending in intergenic.
 * @param d Completed decoder.
 * @param path Receives one class per sequence position.
 * @return 0 on success, otherwise ERANGE.
 */
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
  Decoder *d = decoder_create(length);
  for (int64_t i = 0; i < length; i++)
    decoder_step(d, probability + (size_t)i * ANNO_CLASSES);
  int status = decoder_finish(d, path);
  decoder_destroy(d);
  return status;
}

/* ------------------------------------------------------------------------ */
/* ORF and GFF3 output                                                      */
/* ------------------------------------------------------------------------ */

/**
 * @brief Test whether three bases form a stop codon.
 * @param c Three-base codon.
 * @return 1 for a stop codon, otherwise 0.
 */
static int is_stop(const char *c) {
  return c[0] == 'T' && ((c[1] == 'A' && (c[2] == 'A' || c[2] == 'G')) ||
                         (c[1] == 'G' && c[2] == 'A'));
}

/**
 * @brief Longest ATG..stop ORF (stop included) in any frame of a spliced
 *        transcript; ties keep the 5'-most ORF.
 * @param transcript Spliced transcript bases.
 * @param length Transcript length.
 * @param start Receives the ORF start in transcript coordinates.
 * @param end Receives the exclusive ORF end in transcript coordinates.
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

/**
 * @brief Map a transcript range onto genomic exon pieces.
 * @param exons Genomic exon segments in genomic order.
 * @param strand Transcript strand.
 * @param from Inclusive transcript start.
 * @param to Exclusive transcript end.
 * @param out Destination segment list.
 * @return None.
 */
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

/**
 * @brief Test whether an interval lies within one segment.
 * @param list Segment list to inspect.
 * @param start Inclusive interval start.
 * @param end Exclusive interval end.
 * @return 1 when contained, otherwise 0.
 */
static int inside(const SegmentList *list, int64_t start, int64_t end) {
  for (int i = 0; i < list->count; i++)
    if (start >= list->items[i].start && end <= list->items[i].end) return 1;
  return 0;
}

/**
 * @brief Check ID/Parent links and containment of the records of one gene.
 * @param records Records to validate.
 * @param count Number of records.
 * @param exons Exon segments used for containment checks.
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

/**
 * @brief Append a GFF3 record to a growable record array.
 * @param records Record array to extend.
 * @param count Number of records, updated in place.
 * @param type GFF3 feature type.
 * @param start Inclusive genomic start.
 * @param end Exclusive genomic end.
 * @param phase CDS phase, or -1.
 * @param id Record ID.
 * @param parent Parent record ID, or an empty string.
 * @return None.
 */
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
 * @param output Destination GFF3 stream.
 * @param sequence Source sequence record.
 * @param exons Transcript exon segments.
 * @param strand Transcript strand.
 * @param gene_number Gene number used in generated IDs.
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

/**
 * @brief Map a decoded class to its strand.
 * @param cls Decoded annotation class.
 * @return '+' or '-', or 0 for intergenic.
 */
static char class_strand(uint8_t cls) {
  return cls == ANNO_INTERGENIC ? 0 : cls <= 4 ? '+' : '-';
}

/**
 * @brief Convert a decoded class path into GFF3 transcript records.
 * @param output Destination GFF3 stream.
 * @param sequence Sequence metadata and bases.
 * @param path Per-base decoded classes.
 * @param gene_counter Gene number counter, updated in place.
 * @return Number of transcripts written, or -1 on output error.
 */
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
  FILE *files[ANNO_CLASSES];
  int value[ANNO_CLASSES];
  int64_t start[ANNO_CLASSES];
} Bedgraph;

/**
 * @brief Flush one pending bedGraph run.
 * @param b BedGraph state containing output files and current runs.
 * @param name Sequence name.
 * @param j Class index.
 * @param position Exclusive genomic end position.
 * @return None.
 */
static void bedgraph_emit(Bedgraph *b, const char *name, int j, int64_t position) {
  if (b->value[j] >= 0 && position > b->start[j] &&
      fprintf(b->files[j], "%s\t%lld\t%lld\t%d.%02d\n", name,
              (long long)b->start[j], (long long)position, b->value[j] / 100,
              b->value[j] % 100) < 0)
    fail("cannot write bedGraph");
}

/**
 * @brief Append a probability row to the class bedGraph runs.
 * @param b BedGraph state to update.
 * @param name Sequence name.
 * @param position Genomic position.
 * @param row Class probability row.
 * @return None.
 */
static void bedgraph_add(Bedgraph *b, const char *name, int64_t position,
                          const float *row) {
  if (!b->files[0]) return;
  for (int j = 0; j < ANNO_CLASSES; j++) {
    int value = (int)lrintf(row[j] * 100.0f);
    if (value == b->value[j]) continue;
    bedgraph_emit(b, name, j, position);
    b->value[j] = value;
    b->start[j] = position;
  }
}

/**
 * @brief Predict, decode, and write GFF3 and optional bedGraph files.
 * @param unet Trained model.
 * @param sequences Genome sequences to predict.
 * @param count Number of sequences.
 * @param threads Number of worker threads.
 * @param gff Destination GFF3 stream.
 * @param prefix Optional bedGraph output prefix.
 * @return Number of genes written.
 */
static unsigned long predict_genome(const Unet *unet, const AnnoSequence *sequences,
                                    int count, int threads, FILE *gff,
                                    const char *prefix) {
  (void)threads;
  Bedgraph bedgraph = {0};
  for (int j = 0; prefix && j < ANNO_CLASSES; j++) {
    char path[4096];
    if (snprintf(path, sizeof(path), "%s.%s.bedgraph", prefix, class_names[j]) >=
        (int)sizeof(path)) fail("output prefix too long");
    if (!(bedgraph.files[j] = fopen(path, "w"))) fail("cannot create %s", path);
    fprintf(bedgraph.files[j], "track type=bedGraph name=\"%s\"\n", class_names[j]);
  }
  unsigned long genes = 0;
  for (int i = 0; i < count; i++) {
    const AnnoSequence *s = &sequences[i];
    Decoder *decoder = decoder_create(s->length);
    for (int j = 0; j < ANNO_CLASSES; j++) bedgraph.value[j] = -1;
    int64_t window_count = (s->length + KANN_LENGTH - 1) / KANN_LENGTH;
    float *input = xcalloc((size_t)ANNO_INPUT_CHANNELS * KANN_LENGTH,
                           sizeof(float));
    float *power = xcalloc((size_t)KANN_LENGTH * ANNO_CWT_CHANNELS,
                           sizeof(float));
    float *probability = xcalloc((size_t)KANN_LENGTH * ANNO_CLASSES,
                                 sizeof(float));
    for (int64_t window = 0; window < window_count; window++) {
      int64_t start = window * KANN_LENGTH;
      int64_t valid = s->length - start < KANN_LENGTH ? s->length - start : KANN_LENGTH;
      memset(input, 0, (size_t)ANNO_INPUT_CHANNELS * KANN_LENGTH * sizeof(float));
      memset(power, 0, (size_t)KANN_LENGTH * ANNO_CWT_CHANNELS * sizeof(float));
      cwt_power(s->bases, s->length, start, start + valid, power);
      features(s, power, start, start, valid, 0, input, NULL, NULL);
      memset(probability, 0, (size_t)valid * ANNO_CLASSES * sizeof(float));
      if (unet_predict_channels(unet, input, valid, probability))
        fail("KANN prediction failed");
      for (int64_t p = 0; p < valid; p++) {
        const float *row = probability + (size_t)p * ANNO_CLASSES;
        bedgraph_add(&bedgraph, s->name, start + p, row);
        decoder_step(decoder, row);
      }
    }
    free(input);
    free(power);
    free(probability);
    for (int j = 0; j < ANNO_CLASSES && bedgraph.files[0]; j++)
      bedgraph_emit(&bedgraph, s->name, j, s->length);
    uint8_t *path = xcalloc((size_t)s->length, 1);
    if (decoder_finish(decoder, path)) fail("grammar decoding failed for %s", s->name);
    decoder_destroy(decoder);
    long written = anno_write_gff(gff, s, path, &genes);
    free(path);
    if (written < 0) fail("cannot write GFF3");
    say("%s: %ld genes", s->name, written);
  }
  for (int j = 0; j < ANNO_CLASSES; j++)
    if (bedgraph.files[j] && fclose(bedgraph.files[j])) fail("cannot close bedGraph");
  return genes;
}

/* ------------------------------------------------------------------------ */
/* Command line                                                             */
/* ------------------------------------------------------------------------ */

#ifndef ANNO_NO_MAIN
/**
 * @brief Print command-line usage to standard error.
 * @return None.
 */
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

/**
 * @brief Parse a positive integer option value.
 * @param text Numeric option text.
 * @param name Option name used in diagnostics.
 * @return Parsed positive integer.
 */
static int positive(const char *text, const char *name) {
  char *end = NULL;
  errno = 0;
  long value = strtol(text, &end, 10);
  if (errno || end == text || *end || value < 1 || value > INT_MAX)
    fail("invalid %s: %s", name, text);
  return (int)value;
}

/**
 * @brief Run the anno command-line application.
 * @param argc Number of command-line arguments.
 * @param argv Command-line argument vector.
 * @return EXIT_SUCCESS on success, otherwise EXIT_FAILURE.
 */
int main(int argc, char **argv) {
  int threads = 1, iterations = ANNO_DEFAULT_ITERATIONS;
  const char *model_path = "anno.ann", *prefix = NULL;
  ketopt_t options = KETOPT_INIT;
  int opt;
  while ((opt = ketopt(&options, argc, argv, 1, "hp:i:m:o:", NULL)) >= 0) {
    if (opt == 'h') {
      usage();
      return EXIT_SUCCESS;
    }
    if (opt == 'p') threads = positive(options.arg, "thread count");
    if (opt == 'i') iterations = positive(options.arg, "iteration count");
    if (opt == 'm') model_path = options.arg;
    if (opt == 'o') prefix = options.arg;
  }
  if (opt == ':' || opt == '?') {
    usage();
    return EXIT_FAILURE;
  }
  int a = options.ind;
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
