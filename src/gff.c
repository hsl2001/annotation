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
  char phase_text[2] = {phase < 0 ? '.' : (char)('0' + phase), '\0'};
  return fprintf(output, "%s\tanno\t%s\t%lld\t%lld\t.\t%c\t%s\tID=%s%s%s\n",
                 sequence->name, type, (long long)start + 1, (long long)end,
                 strand, phase_text, id, *parent ? ";Parent=" : "", parent) > 0;
}

static char coding_base(const Sequence *sequence, int64_t p, char strand) {
  return strand == '+' ? sequence->bases[p]
                       : complement(sequence->bases[sequence->length - 1 - p]);
}

static int valid_form(const Sequence *sequence, const IntervalList *form,
                      char strand) {
  int64_t length = 0, at = 0;
  char codon[3];
  if (!form->count)
    return 0;
  for (int i = 0; i < form->count; ++i) {
    const Interval *e = &form->items[i];
    if (e->start < 0 || e->end > sequence->length || e->start >= e->end ||
        (i && e->start <= form->items[i - 1].end))
      return 0;
    length += e->end - e->start;
  }
  if (length < 6 || length % 3)
    return 0;
  for (int i = 0; i < form->count; ++i)
    for (int64_t p = form->items[i].start; p < form->items[i].end; ++p) {
      char c = coding_base(sequence, p, strand);
      if (c != 'A' && c != 'C' && c != 'G' && c != 'T')
        return 0;
      codon[at % 3] = c;
      ++at;
      if (at % 3)
        continue;
      if (at == 3 && (codon[0] != 'A' || codon[1] != 'T' || codon[2] != 'G'))
        return 0;
      if (is_stop(codon) != (at == length))
        return 0;
    }
  return 1;
}

int write_isoform(FILE *output, const Sequence *sequence,
                  const IntervalList *form, char strand,
                  unsigned long gene_number, size_t number) {
  if (!number || (strand != '+' && strand != '-') ||
      !valid_form(sequence, form, strand))
    return 0;
  char gene[48], mrna[80], id[112];
  snprintf(gene, sizeof(gene), "gene%lu", gene_number);
  int64_t begin = form->items[0].start;
  int64_t end = form->items[form->count - 1].end;
  int64_t a = strand == '+' ? begin : sequence->length - end;
  int64_t b = strand == '+' ? end : sequence->length - begin;
  if (number == 1 &&
      !write_record(output, sequence, "gene", a, b, strand, -1, gene, ""))
    return -1;
  snprintf(mrna, sizeof(mrna), "%s.t%zu", gene, number);
  if (!write_record(output, sequence, "mRNA", a, b, strand, -1, mrna, gene))
    return -1;
  int64_t length = 0;
  for (int e = 0; e < form->count; ++e) {
    const Interval *exon = &form->items[e];
    a = strand == '+' ? exon->start : sequence->length - exon->end;
    b = strand == '+' ? exon->end : sequence->length - exon->start;
    snprintf(id, sizeof(id), "%s.exon%d", mrna, e + 1);
    if (!write_record(output, sequence, "exon", a, b, strand, -1, id, mrna))
      return -1;
    snprintf(id, sizeof(id), "%s.cds%d", mrna, e + 1);
    if (!write_record(output, sequence, "CDS", a, b, strand,
                      (int)((3 - length % 3) % 3), id, mrna))
      return -1;
    length += exon->end - exon->start;
  }
  return 1;
}
