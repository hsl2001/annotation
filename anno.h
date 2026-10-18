#ifndef ANNO_H
#define ANNO_H

#include <complex.h>
#include <limits.h>
#include <stddef.h>
#include <stdint.h>

#ifndef CWT_SCALES
#define CWT_SCALES 3, 4, 5
#endif
#define CWT_COUNT_IMPL(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, \
                       _17, _18, _19, _20, _21, _22, _23, _24, _25, _26, _27, _28, _29, _30, \
                       _31, _32, count, ...) count
#define CWT_COUNT(...) CWT_COUNT_IMPL(__VA_ARGS__, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, \
                                      21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, \
                                      6, 5, 4, 3, 2, 1)
#define WAVE_COUNT CWT_COUNT(CWT_SCALES)
#ifndef MAX_WAVE_SIZE
#define MAX_WAVE_SIZE 5
#endif
#define CWT_WINDOW_SIZE 1024
#define CWT_CHANNELS (2 * WAVE_COUNT)

#define MIN_CDS 12
#define DEFAULT_ITERATIONS 500

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
