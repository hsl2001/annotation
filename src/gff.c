#include "include/surf.h"
#include "include/utils.h"

#include <stdint.h>
#include <stdio.h>

int is_stop(const char *c) {
  return c[0] == 'T' && ((c[1] == 'A' && (c[2] == 'A' || c[2] == 'G')) ||
                         (c[1] == 'G' && c[2] == 'A'));
}

static int write_record(FILE *output, const Sequence *sequence,
                        const char *type, int64_t start, int64_t end,
                        char strand, int phase, const char *id,
                        const char *parent) {
  int64_t a = strand == '+' ? start : sequence->length - end;
  int64_t b = strand == '+' ? end : sequence->length - start;
  char phase_text[2] = {phase < 0 ? '.' : (char)('0' + phase), '\0'};
  return fprintf(output, "%s\tsurf\t%s\t%lld\t%lld\t.\t%c\t%s\tID=%s%s%s\n",
                 sequence->name, type, (long long)a + 1, (long long)b, strand,
                 phase_text, id, *parent ? ";Parent=" : "", parent) > 0;
}

static char coding_base(const Sequence *sequence, int64_t p, char strand) {
  return strand == '+' ? sequence->bases[p]
                       : complement(sequence->bases[sequence->length - 1 - p]);
}

Interval longest_orf(const Sequence *sequence, const IntervalList *form,
                     char strand) {
  Interval best = {0, 0};
  int64_t at = 0, starts[3] = {-1, -1, -1};
  char codon[3] = {0};
  if (!sequence || !sequence->bases || !form || !form->items ||
      form->count <= 0)
    return best;
  for (int i = 0; i < form->count; ++i) {
    const Interval *e = &form->items[i];
    if (e->start < 0 || e->end > sequence->length || e->start >= e->end ||
        (i && e->start < form->items[i - 1].end))
      return (Interval){0, 0};
    for (int64_t p = e->start; p < e->end; ++p) {
      codon[0] = codon[1];
      codon[1] = codon[2];
      codon[2] = coding_base(sequence, p, strand);
      ++at;
      if (base_index(codon[2]) < 0)
        starts[0] = starts[1] = starts[2] = -1;
      if (at < 3)
        continue;
      int frame = (int)(at % 3);
      int64_t start = starts[frame];
      if (start < 0 && codon[0] == 'A' && codon[1] == 'T' && codon[2] == 'G')
        starts[frame] = at - 3;
      if (is_stop(codon)) {
        int64_t length = start < 0 ? 0 : at - start;
        if (start >= 0 &&
            (length > best.end - best.start ||
             (length == best.end - best.start && start < best.start)))
          best = (Interval){start, at};
        starts[frame] = -1;
      }
    }
  }
  return best;
}

int write_isoform(FILE *output, const Sequence *sequence,
                  const IntervalList *form, char strand,
                  unsigned long gene_number, size_t number) {
  if (!number || (strand != '+' && strand != '-'))
    return 0;
  Interval cds = longest_orf(sequence, form, strand);
  // Only the representative requires an ORF; alternatives use exon overlap.
  if (cds.start == cds.end && number == 1)
    return 0;
  if (!output)
    return 1;
  char gene[48], mrna[80], id[112];
  snprintf(gene, sizeof(gene), "gene%lu", gene_number);
  int64_t begin = form->items[0].start;
  int64_t end = form->items[form->count - 1].end;
  if (number == 1 &&
      !write_record(output, sequence, "gene", begin, end, strand, -1, gene, ""))
    return -1;
  snprintf(mrna, sizeof(mrna), "%s.t%zu", gene, number);
  if (!write_record(output, sequence, "mRNA", begin, end, strand, -1, mrna,
                    gene))
    return -1;
  const char *types[] = {"five_prime_UTR", "CDS", "three_prime_UTR"};
  int64_t offset = 0, coding = 0;
  for (int e = 0; e < form->count; ++e) {
    const Interval *exon = &form->items[e];
    snprintf(id, sizeof(id), "%s.exon%d", mrna, e + 1);
    if (!write_record(output, sequence, "exon", exon->start, exon->end, strand,
                      -1, id, mrna))
      return -1;
    if (cds.start == cds.end)
      continue;
    int64_t length = exon->end - exon->start;
    int64_t cuts[] = {0, cds.start - offset, cds.end - offset, length};
    for (int t = 0; t < 3; ++t) {
      int64_t lo = cuts[t] > 0 ? cuts[t] : 0;
      int64_t hi = cuts[t + 1] < length ? cuts[t + 1] : length;
      if (lo >= hi)
        continue;
      snprintf(id, sizeof(id), "%s.%s%d", mrna, t == 1 ? "cds" : types[t],
               e + 1);
      int phase = t == 1 ? (int)((3 - coding % 3) % 3) : -1;
      if (!write_record(output, sequence, types[t], exon->start + lo,
                        exon->start + hi, strand, phase, id, mrna))
        return -1;
      if (t == 1)
        coding += hi - lo;
    }
    offset += length;
  }
  return 1;
}
