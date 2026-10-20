#ifndef ANNO_H
#define ANNO_H

#include <complex.h>
#include <limits.h>
#include <stddef.h>
#include <stdint.h>

#define CWT_SCALES 4, 5, 6
#define WAVE_COUNT 3
#define CWT_WINDOW_SIZE 1024
#define CWT_CHANNELS (2 * WAVE_COUNT)

#define MIN_CDS 12
#define DEFAULT_ITERATIONS 100

typedef struct {
   double _Complex *kernel[WAVE_COUNT];
   int widths[WAVE_COUNT], max_width;
} Wavelets;
typedef struct {
   char *name, *seq;
   int length, genome_length, offset;
   float *cwt_power;
   uint8_t *label;
} Contig;

extern const int wave_sizes[WAVE_COUNT];
void anno_fail(const char *format, ...);
void *anno_alloc(size_t count, size_t size);

static inline int base_index(char base) {
  switch (base) {
  case 'A': case 'a': return 0;
  case 'C': case 'c': return 1;
  case 'G': case 'g': return 2;
  case 'T': case 't': return 3;
  default: return -1;
  }
}

Contig *read_fasta(const char *path, int *count);
int label_cds(Contig *contigs, int count, const char *gff, int *skipped);
void write_gff(const Contig *contig, const uint8_t *path, unsigned long *genes);

#endif
