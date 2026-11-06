/**
 * @file gff.c
 * @brief Logics for gene annotation decoding
 */

#include "include/surf.h"
#include "include/utils.h"
#include <stdlib.h>

/**
 * @brief Map a decoded class to its strand.
 * @param[in] cls Decoded annotation class.
 * @return '+' or '-', or 0 for intergenic.
 */
static char strand_of_class(uint8_t class) {
  switch (class) {
  case PLUS_START_SINGLE_EXON:
  case PLUS_END_EXON:
  case PLUS_INTRON:
  case PLUS_INTERNAL_EXON:
    return '+';
  case MINUS_START_SINGLE_EXON:
  case MINUS_END_EXON:
  case MINUS_INTRON:
  case MINUS_INTERNAL_EXON:
    return '-';
  default:
    return 0;
  }
}

/**
 * @brief Test whether three bases form a stop codon.
 * @param[in] c: Three-base codon.
 * @return 1 for a stop codon, otherwise 0.
 */
int is_stop(const char *c) {
  return c[0] == 'T' && ((c[1] == 'A' && (c[2] == 'A' || c[2] == 'G')) ||
                         (c[1] == 'G' && c[2] == 'A'));
}

/**
 * @brief Longest ATG..STOP ORF (stop included) in any frame of a spliced
 *        transcript; ties keep the 5'-most ORF.
 * @param[in] transcript: !!Spliced!! transcript bases.
 * @param[in] length: Transcript length.
 * @param[out] start: Receives the ORF start in transcript coordinates.
 * @param[out] end: Receives the exclusive ORF end in transcript coordinates.
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

static int write_record(FILE *output, const Sequence *sequence,
                        const char *type, int64_t start, int64_t end,
                        char strand, int phase, const char *id,
                        const char *parent) {
  char phase_text[2] = {phase < 0 ? '.' : (char)('0' + phase), '\0'};
  return fprintf(output, "%s\tanno\t%s\t%lld\t%lld\t.\t%c\t%s\tID=%s%s%s\n",
                 sequence->name, type, (long long)start + 1, (long long)end,
                 strand, phase_text, id, *parent ? ";Parent=" : "", parent) > 0;
}

static int write_region(FILE *output, const Sequence *sequence,
                        const IntervalList *exons, char strand, int64_t length,
                        int64_t from, int64_t to, const char *type,
                        const char *id_type, const char *parent, int coding) {
  int64_t offset = strand == '+' ? 0 : length;
  int count = 0;
  for (int i = 0; i < exons->count; i++) {
    const Interval *exon = &exons->items[i];
    int64_t width = exon->end - exon->start;
    if (strand == '-')
      offset -= width;
    int64_t a = from > offset ? from : offset;
    int64_t b = to < offset + width ? to : offset + width;
    if (a < b) {
      char id[96];
      snprintf(id, sizeof(id), "%s.%s%d", parent, id_type, ++count);
      if (!write_record(output, sequence, type,
                        strand == '+' ? exon->start + a - offset
                                      : exon->end - (b - offset),
                        strand == '+' ? exon->start + b - offset
                                      : exon->end - (a - offset),
                        strand, coding ? (int)((3 - (a - from) % 3) % 3) : -1,
                        id, parent))
        return 0;
    }
    if (strand == '+')
      offset += width;
  }
  return 1;
}

/**
 * @brief Emit one transcript (exons on one strand) with its longest ORF.
 * @param[in,out] output Destination GFF3 stream.
 * @param[in] sequence Source sequence record.
 * @param[in] exons Transcript exon intervals.
 * @param[in] strand Transcript strand.
 * @param[in] gene_number Gene number used in generated IDs.
 * @return 1 written, 0 skipped (no ORF), -1 on error.
 */
static int write_one_transcript(FILE *output, const Sequence *sequence,
                                const IntervalList *exons, char strand,
                                unsigned long gene_number) {
  int64_t length;
  char *transcript = splice_sequence(sequence, exons, strand, &length);
  int64_t orf_start = INT64_MAX, orf_end = 0;
  int found = longest_orf(transcript, length, &orf_start, &orf_end);
  free(transcript);
  if (!found)
    return 0;
  char gene[48], mrna[64], id[96];
  snprintf(gene, sizeof(gene), "gene%lu", gene_number);
  snprintf(mrna, sizeof(mrna), "%s.t1", gene);
  int64_t start = exons->items[0].start,
          end = exons->items[exons->count - 1].end;
  int ok = write_record(output, sequence, "gene", start, end, strand, -1, gene,
                        "") &&
           write_record(output, sequence, "mRNA", start, end, strand, -1, mrna,
                        gene);
  for (int i = 0; i < exons->count; i++) {
    snprintf(id, sizeof(id), "%s.exon%d", mrna, i + 1);
    ok = ok && write_record(output, sequence, "exon", exons->items[i].start,
                            exons->items[i].end, strand, -1, id, mrna);
  }
  ok = ok &&
       write_region(output, sequence, exons, strand, length, 0, orf_start,
                    "five_prime_UTR", "utr5p", mrna, 0) &&
       write_region(output, sequence, exons, strand, length, orf_start, orf_end,
                    "CDS", "cds", mrna, 1) &&
       write_region(output, sequence, exons, strand, length, orf_end, length,
                    "three_prime_UTR", "utr3p", mrna, 0);
  return ok ? 1 : -1;
}

/**
 * @brief Convert a decoded class path into GFF3 transcript records.
 * @param[in,out] output: Destination GFF3 stream.
 * @param[in] sequence: Sequence metadata and bases.
 * @param[in] path: Per-base decoded classes.
 * @param[in,out] gene_counter: Gene number counter, updated in place.
 * @return Number of transcripts written, or -1 on output error.
 */
long write_gff(FILE *output, const Sequence *sequence, const uint8_t *path,
               unsigned long *gene_counter) {
  long written = 0;
  for (int64_t p = 0; p < sequence->length;) {
    char strand = strand_of_class(path[p]);
    if (!strand) {
      p++;
      continue;
    }
    IntervalList exons = {0};
    while (p < sequence->length && strand_of_class(path[p]) == strand) {
      int intron = path[p] == PLUS_INTRON || path[p] == MINUS_INTRON;
      int64_t start = p;
      while (p < sequence->length && strand_of_class(path[p]) == strand &&
             (path[p] == PLUS_INTRON || path[p] == MINUS_INTRON) == intron)
        p++;
      if (!intron)
        push_interval(&exons, start, p);
    }
    int result = exons.count ? write_one_transcript(output, sequence, &exons,
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
