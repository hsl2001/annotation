/**
 * @file gff.c
 * @brief ORF selection, transcript-coordinate mapping, and GFF3 emission.
 */
#include "anno_internal.h"

#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------------ */
/* ORF and GFF3 output                                                      */
/* ------------------------------------------------------------------------ */

/**
 * @brief Test whether three bases form a stop codon.
 * @param[in] c Three-base codon.
 * @return 1 for a stop codon, otherwise 0.
 */
int is_stop(const char *c) {
  return c[0] == 'T' && ((c[1] == 'A' && (c[2] == 'A' || c[2] == 'G')) ||
                         (c[1] == 'G' && c[2] == 'A'));
}

/**
 * @brief Longest ATG..stop ORF (stop included) in any frame of a spliced
 *        transcript; ties keep the 5'-most ORF.
 * @param[in] transcript Spliced transcript bases.
 * @param[in] length Transcript length.
 * @param[out] start Receives the ORF start in transcript coordinates.
 * @param[out] end Receives the exclusive ORF end in transcript coordinates.
 * @return 1 when found, filling [*start,*end) in transcript coordinates.
 */
int longest_orf(const char *transcript, int64_t length, int64_t *start,
                       int64_t *end) {
  int64_t best = 0;
  for (int frame = 0; frame < 3; frame++) {
    int64_t open = -1;
    for (int64_t p = frame; p + 3 <= length; p += 3) {
      const char *c = transcript + p;
      if (open < 0 && c[0] == 'A' && c[1] == 'T' && c[2] == 'G')
        open = p;
      if (!is_stop(c))
        continue;
      if (open >= 0 &&
          (p + 3 - open > best || (p + 3 - open == best && open < *start))) {
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
 * @param[in] exons Genomic exon segments in genomic order.
 * @param[in] strand Transcript strand.
 * @param[in] from Inclusive transcript start.
 * @param[in] to Exclusive transcript end.
 * @param[out] out Destination segment list.
 * @return None.
 */
static void map_range(const SegmentList *exons, char strand, int64_t from,
                      int64_t to, SegmentList *out) {
  int64_t offset = 0;
  for (int order = 0; order < exons->count; order++) {
    const Segment *e =
        &exons->items[strand == '+' ? order : exons->count - 1 - order];
    int64_t length = e->end - e->start;
    int64_t a = from > offset ? from - offset : 0;
    int64_t b = to < offset + length ? to - offset : length;
    if (a < b)
      anno_segments_push(out, strand == '+' ? e->start + a : e->end - b,
                    strand == '+' ? e->start + b : e->end - a);
    offset += length;
  }
  qsort(out->items, (size_t)out->count, sizeof(Segment), anno_segment_compare);
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
 * @param[in] list Segment list to inspect.
 * @param[in] start Inclusive interval start.
 * @param[in] end Exclusive interval end.
 * @return 1 when contained, otherwise 0.
 */
static int inside(const SegmentList *list, int64_t start, int64_t end) {
  for (int i = 0; i < list->count; i++)
    if (start >= list->items[i].start && end <= list->items[i].end)
      return 1;
  return 0;
}

/**
 * @brief Check ID/Parent links and containment of the records of one gene.
 * @param[in] records Records to validate.
 * @param[in] count Number of records.
 * @param[in] exons Exon segments used for containment checks.
 * @return 1 when consistent.
 */
static int records_valid(const Record *records, int count,
                         const SegmentList *exons) {
  for (int i = 0; i < count; i++) {
    const Record *r = &records[i];
    if (r->start < 0 || r->end <= r->start)
      return 0;
    if (i == 0) {
      if (strcmp(r->type, "gene") || r->parent[0])
        return 0;
      continue;
    }
    int parent = -1;
    for (int j = 0; j < count; j++)
      if (r->parent[0] && !strcmp(records[j].id, r->parent))
        parent = j;
    if (parent < 0 || r->start < records[parent].start ||
        r->end > records[parent].end)
      return 0;
    if (strcmp(r->type, "mRNA") && strcmp(r->type, "exon") &&
        !inside(exons, r->start, r->end))
      return 0;
    for (int j = 0; j < i; j++)
      if (r->id[0] && !strcmp(records[j].id, r->id))
        return 0;
  }
  return 1;
}

/**
 * @brief Append a GFF3 record to a growable record array.
 * @param[in,out] records Record array to extend.
 * @param[in,out] count Number of records, updated in place.
 * @param[in] type GFF3 feature type.
 * @param[in] start Inclusive genomic start.
 * @param[in] end Exclusive genomic end.
 * @param[in] phase CDS phase, or -1.
 * @param[in] id Record ID.
 * @param[in] parent Parent record ID, or an empty string.
 * @return None.
 */
static void record_add(Record **records, int *count, const char *type,
                       int64_t start, int64_t end, int phase, const char *id,
                       const char *parent) {
  *records = anno_xrealloc(*records, (size_t)*count + 1, sizeof(**records));
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
 * @param[in,out] output Destination GFF3 stream.
 * @param[in] sequence Source sequence record.
 * @param[in] exons Transcript exon segments.
 * @param[in] strand Transcript strand.
 * @param[in] gene_number Gene number used in generated IDs.
 * @return 1 written, 0 skipped (no ORF), -1 on error.
 */
static int write_transcript(FILE *output, const AnnoSequence *sequence,
                            const SegmentList *exons, char strand,
                            unsigned long gene_number) {
  int64_t length = 0;
  for (int i = 0; i < exons->count; i++)
    length += exons->items[i].end - exons->items[i].start;
  char *transcript = anno_xcalloc((size_t)length + 1, 1);
  int64_t at = 0;
  for (int order = 0; order < exons->count; order++) {
    const Segment *e =
        &exons->items[strand == '+' ? order : exons->count - 1 - order];
    for (int64_t p = 0; p < e->end - e->start; p++)
      transcript[at++] = strand == '+'
                             ? sequence->bases[e->start + p]
                             : anno_complement(sequence->bases[e->end - 1 - p]);
  }
  int64_t orf_start = INT64_MAX, orf_end = 0;
  int found = longest_orf(transcript, length, &orf_start, &orf_end);
  free(transcript);
  if (!found)
    return 0;
  SegmentList utr5 = {0}, cds = {0}, utr3 = {0};
  map_range(exons, strand, 0, orf_start, &utr5);
  map_range(exons, strand, orf_start, orf_end, &cds);
  map_range(exons, strand, orf_end, length, &utr3);
  char gene[48], mrna[64], id[96];
  snprintf(gene, sizeof(gene), "gene%lu", gene_number);
  snprintf(mrna, sizeof(mrna), "%s.t1", gene);
  Record *records = NULL;
  int count = 0;
  int64_t start = exons->items[0].start,
          end = exons->items[exons->count - 1].end;
  record_add(&records, &count, "gene", start, end, -1, gene, "");
  record_add(&records, &count, "mRNA", start, end, -1, mrna, gene);
  for (int i = 0; i < exons->count; i++) {
    snprintf(id, sizeof(id), "%s.exon%d", mrna, i + 1);
    record_add(&records, &count, "exon", exons->items[i].start,
               exons->items[i].end, -1, id, mrna);
  }
  for (int i = 0; i < utr5.count; i++) {
    snprintf(id, sizeof(id), "%s.utr5p%d", mrna, i + 1);
    record_add(&records, &count, "five_prime_UTR", utr5.items[i].start,
               utr5.items[i].end, -1, id, mrna);
  }
  for (int i = 0; i < cds.count; i++) {
    int64_t before = 0; /* coding bases upstream in transcript orientation */
    for (int j = 0; j < cds.count; j++)
      if (strand == '+' ? j < i : j > i)
        before += cds.items[j].end - cds.items[j].start;
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
  if (!ok)
    anno_fail("internal error: inconsistent GFF3 records for %s", gene);
  for (int i = 0; ok && i < count; i++) {
    const Record *r = &records[i];
    char phase[12] = ".";
    if (r->phase >= 0)
      snprintf(phase, sizeof(phase), "%d", r->phase);
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
 * @param[in] cls Decoded annotation class.
 * @return '+' or '-', or 0 for intergenic.
 */
static char class_strand(uint8_t cls) {
  switch (cls) {
  case ANNO_PLUS_START_SINGLE_EXON:
  case ANNO_PLUS_END_EXON:
  case ANNO_PLUS_INTRON:
  case ANNO_PLUS_INTERNAL_EXON: return '+';
  case ANNO_MINUS_START_SINGLE_EXON:
  case ANNO_MINUS_END_EXON:
  case ANNO_MINUS_INTRON:
  case ANNO_MINUS_INTERNAL_EXON: return '-';
  default: return 0;
  }
}

/**
 * @brief Convert a decoded class path into GFF3 transcript records.
 * @param[in,out] output Destination GFF3 stream.
 * @param[in] sequence Sequence metadata and bases.
 * @param[in] path Per-base decoded classes.
 * @param[in,out] gene_counter Gene number counter, updated in place.
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
             (path[p] == ANNO_PLUS_INTRON || path[p] == ANNO_MINUS_INTRON) ==
                 intron)
        p++;
      if (!intron)
        anno_segments_push(&exons, start, p);
    }
    int result = exons.count ? write_transcript(output, sequence, &exons,
                                                strand, *gene_counter + 1)
                             : 0;
    free(exons.items);
    if (result < 0)
      return -1;
    if (result) {
      ++*gene_counter;
      written++;
    }
  }
  return written;
}

