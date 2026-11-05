/**
 * @file labels.c
 * @brief GFF3 transcript grouping and per-base label construction.
 */
#define _POSIX_C_SOURCE 200809L

#include "include/anno_internal.h"

#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

/* ------------------------------------------------------------------------ */
/* GFF3 -> label vector                                                     */
/* ------------------------------------------------------------------------ */

enum { KIND_EXON, KIND_INTRON };

/** One GFF feature attached to one transcript. Owns `parent`. */
typedef struct {
  char *parent;
  int sequence;
  int kind;
  char strand;
  int64_t start, end; /* 0-based half-open */
} GffFeature;

/**
 * @brief Compare GFF features by parent, kind, and start.
 * @param[in] a First feature.
 * @param[in] b Second feature.
 * @return Negative, zero, or positive comparison result.
 */
static int feature_compare(const void *a, const void *b) {
  const GffFeature *x = a, *y = b;
  int order = strcmp(x->parent, y->parent);
  if (order)
    return order;
  if (x->kind != y->kind)
    return x->kind - y->kind;
  return x->start < y->start ? -1 : x->start > y->start;
}

/**
 * @brief Convert a GFF feature type to an internal kind.
 * @param[in] type GFF feature type.
 * @return KIND_* value, -2 for container features, or -1 when ignored.
 */
static int feature_kind(const char *type) {
  if (!strcmp(type, "exon"))
    return KIND_EXON;
  if (!strcmp(type, "intron"))
    return KIND_INTRON;
  if (!strcmp(type, "gene") || !strcmp(type, "mRNA") ||
      !strcmp(type, "transcript"))
    return -2;
  return -1;
}

/**
 * @brief Parse a non-negative decimal integer.
 * @param[in] text Decimal text to parse.
 * @param[out] value Receives the parsed integer.
 * @return 1 on success, otherwise 0.
 */
static int parse_int64(const char *text, int64_t *value) {
  char *end = NULL;
  errno = 0;
  long long parsed = strtoll(text, &end, 10);
  if (errno || end == text || *end || parsed < 0)
    return 0;
  *value = parsed;
  return 1;
}

/** One coding transcript ready to be painted. Owns its lists. */
typedef struct {
  int sequence;
  char strand;
  int64_t start, end;
  SegmentList exons, explicit_introns;
} LabelTranscript;

/**
 * @brief Compare transcripts by decreasing genomic span.
 * @param[in] a First transcript.
 * @param[in] b Second transcript.
 * @return Negative, zero, or positive comparison result.
 */
static int transcript_compare(const void *a, const void *b) {
  const LabelTranscript *x = a, *y = b;
  int64_t lx = x->end - x->start, ly = y->end - y->start;
  if (lx != ly)
    return lx > ly ? -1 : 1;
  if (x->sequence != y->sequence)
    return x->sequence - y->sequence;
  return x->start < y->start ? -1 : x->start > y->start;
}

/**
 * @brief Merge touching or overlapping segments in place.
 * @param[in,out] list Segment list to sort and merge.
 * @return None.
 */
static void segments_merge(SegmentList *list) {
  if (!list->count)
    return;
  qsort(list->items, (size_t)list->count, sizeof(Segment), anno_segment_compare);
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
 * @param[in] path GFF3 file path.
 * @param[in] sequences FASTA sequences used to resolve sequence IDs.
 * @param[in] count Number of FASTA sequences.
 * @param[out] feature_count Receives the number of features.
 * @return Caller-owned feature array.
 */
static GffFeature *read_gff_features(const char *path,
                                     const AnnoSequence *sequences, int count,
                                     size_t *feature_count) {
  gzFile input = gzopen(path, "rb");
  if (!input)
    anno_fail("cannot open annotation %s", path);
  GffFeature *features = NULL;
  size_t used = 0, capacity = 0;
  char line[65536];
  long line_number = 0;
  while (gzgets(input, line, sizeof(line))) {
    line_number++;
    if (!strncmp(line, "##FASTA", 7))
      break;
    if (line[0] == '#' || line[0] == '\n' || line[0] == '\r')
      continue;
    line[strcspn(line, "\r\n")] = '\0';
    char *fields[9];
    char *cursor = line;
    int columns = 0;
    while (columns < 9 && cursor) {
      fields[columns++] = cursor;
      cursor = strchr(cursor, '\t');
      if (cursor)
        *cursor++ = '\0';
    }
    if (columns != 9 || cursor)
      anno_fail("%s:%ld: expected 9 columns", path, line_number);
    int kind = feature_kind(fields[2]);
    if (kind == -1)
      continue;
    int sequence = -1;
    for (int i = 0; i < count; i++)
      if (!strcmp(sequences[i].name, fields[0]))
        sequence = i;
    if (sequence < 0)
      anno_fail("%s:%ld: sequence %s absent from FASTA", path, line_number,
           fields[0]);
    int64_t start = 0, end = 0;
    if (!parse_int64(fields[3], &start) || !parse_int64(fields[4], &end) ||
        start < 1 || end < start)
      anno_fail("%s:%ld: invalid coordinates", path, line_number);
    if (end > sequences[sequence].length)
      anno_fail("%s:%ld: coordinate beyond %s length %lld", path, line_number,
           fields[0], (long long)sequences[sequence].length);
    char strand = fields[6][0];
    if ((strand != '+' && strand != '-') || fields[6][1])
      anno_fail("%s:%ld: unsupported strand %s", path, line_number, fields[6]);
    if (kind < 0)
      continue;
    char *parents = NULL;
    for (char *field = strtok(fields[8], ";"); field;
         field = strtok(NULL, ";")) {
      while (*field == ' ')
        field++;
      if (!strncmp(field, "Parent=", 7))
        parents = field + 7;
    }
    if (!parents || !*parents)
      anno_fail("%s:%ld: %s without Parent", path, line_number, fields[2]);
    for (char *save = NULL, *parent = strtok_r(parents, ",", &save); parent;
         parent = strtok_r(NULL, ",", &save)) {
      if (used == capacity) {
        capacity = capacity ? capacity * 2 : 1024;
        features = anno_xrealloc(features, capacity, sizeof(*features));
      }
      features[used++] =
          (GffFeature){anno_xstrdup(parent), sequence, kind, strand, start - 1, end};
    }
  }
  if (gzclose(input) != Z_OK)
    anno_fail("cannot read annotation %s", path);
  *feature_count = used;
  return features;
}

/**
 * @brief Paint one transcript into a per-base label array.
 * @param[in] t Transcript to paint.
 * @param[in,out] labels Destination labels; overlapping transcripts are prechecked.
 * @return None.
 */
static void paint_transcript(const LabelTranscript *t, uint8_t *labels) {
  int plus = t->strand == '+';
  uint8_t intron = plus ? ANNO_PLUS_INTRON : ANNO_MINUS_INTRON;
  for (int64_t i = t->start; i < t->end; i++)
    labels[i] = intron;
  for (int e = 0; e < t->exons.count; e++) {
    int first = e == (plus ? 0 : t->exons.count - 1);
    int last = e == (plus ? t->exons.count - 1 : 0);
    uint8_t exon = first ? (plus ? ANNO_PLUS_START_SINGLE_EXON
                                : ANNO_MINUS_START_SINGLE_EXON)
                   : last ? (plus ? ANNO_PLUS_END_EXON : ANNO_MINUS_END_EXON)
                          : (plus ? ANNO_PLUS_INTERNAL_EXON
                                  : ANNO_MINUS_INTERNAL_EXON);
    for (int64_t i = t->exons.items[e].start; i < t->exons.items[e].end; i++)
      labels[i] = exon;
  }
  for (int e = 0; e < t->explicit_introns.count; e++)
    for (int64_t i = t->explicit_introns.items[e].start;
         i < t->explicit_introns.items[e].end; i++)
      labels[i] = intron;
}

/**
 * @brief Read GFF3 annotations and fill sequence label arrays.
 * @param[in] path GFF3 file path.
 * @param[in,out] sequences FASTA sequence records to label.
 * @param[in] count Number of sequence records.
 * @return None.
 */
void anno_read_labels(const char *path, AnnoSequence *sequences, int count) {
  size_t feature_count = 0;
  GffFeature *features =
      read_gff_features(path, sequences, count, &feature_count);
  qsort(features, feature_count, sizeof(*features), feature_compare);
  LabelTranscript *transcripts = NULL;
  size_t transcript_count = 0, skipped = 0;
  for (size_t first = 0; first < feature_count;) {
    size_t last = first;
    while (last < feature_count &&
           !strcmp(features[last].parent, features[first].parent))
      last++;
    LabelTranscript t = {.sequence = features[first].sequence,
                         .strand = features[first].strand};
    for (size_t i = first; i < last; i++) {
      const GffFeature *f = &features[i];
      if (f->sequence != t.sequence || f->strand != t.strand)
        anno_fail("transcript %s mixes sequences or strands", f->parent);
      if (f->kind == KIND_EXON)
        anno_segments_push(&t.exons, f->start, f->end);
      if (f->kind == KIND_INTRON)
        anno_segments_push(&t.explicit_introns, f->start, f->end);
    }
    segments_merge(&t.exons);
    if (!t.exons.count) {
      free(t.exons.items);
      free(t.explicit_introns.items);
      first = last;
      continue;
    }
    t.start = t.exons.items[0].start;
    t.end = t.exons.items[t.exons.count - 1].end;
    transcripts =
        anno_xrealloc(transcripts, transcript_count + 1, sizeof(*transcripts));
    transcripts[transcript_count++] = t;
    first = last;
  }
  for (size_t i = 0; i < feature_count; i++)
    free(features[i].parent);
  free(features);
  qsort(transcripts, transcript_count, sizeof(*transcripts),
        transcript_compare);
  /* Reserve whole transcript spans, including their introns. */
  uint8_t **occupied = anno_xcalloc((size_t)count, sizeof(*occupied));
  for (int i = 0; i < count; i++) {
    free(sequences[i].labels);
    sequences[i].labels = anno_xcalloc((size_t)sequences[i].length, 1);
    occupied[i] = anno_xcalloc((size_t)sequences[i].length, 1);
  }
  for (size_t k = 0; k < transcript_count; k++) {
    LabelTranscript *t = &transcripts[k];
    uint8_t *labels = sequences[t->sequence].labels;
    int free_span = 1;
    for (int64_t i = t->start; i < t->end && free_span; i++)
      free_span = !occupied[t->sequence][i];
    if (free_span) {
      paint_transcript(t, labels);
      memset(occupied[t->sequence] + t->start, 1, (size_t)(t->end - t->start));
    } else
      skipped++;
    free(t->exons.items);
    free(t->explicit_introns.items);
  }
  anno_log("labelled %zu transcripts (%zu overlapping skipped)",
      transcript_count - skipped, skipped);
  free(transcripts);
  for (int i = 0; i < count; i++)
    free(occupied[i]);
  free(occupied);
}
