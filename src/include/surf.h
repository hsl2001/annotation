#ifndef SURF_H
#define SURF_H

#include <stdint.h>
#include <stdio.h>

#define NUCLEOTIDES 4
#define PATH_BUFFER 4096
#define DEFAULT_ITERATIONS 1000

typedef struct {
  int64_t start, end;
} Interval;

typedef struct {
  Interval *items;
  int count, capacity;
} IntervalList;

typedef struct {
  char *name, *bases;
  int64_t length;
  uint8_t *labels;
  int64_t max_intron, max_span;
  IntervalList introns[2], spans[2];
} Sequence;

void push_interval(IntervalList *list, int64_t start, int64_t end);
int compare_interval(const void *a, const void *b);

void read_fasta(const char *path, int *count, Sequence **out);
void free_sequences(Sequence *sequences, int count);
void read_gff3(const char *path, Sequence *sequences, int count);
int is_stop(const char *codon);
Interval longest_orf(const Sequence *sequence, const IntervalList *form,
                     char strand);
int write_isoform(FILE *output, const Sequence *sequence,
                  const IntervalList *form, char strand,
                  unsigned long gene_number, size_t number);

#endif
