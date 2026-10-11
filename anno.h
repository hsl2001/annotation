#ifndef ANNO_H
#define ANNO_H

#include <complex.h>
#include <stddef.h>
#include <stdint.h>

#define WAVE_COUNT 10
#ifndef MAX_WAVE_SIZE
#define MAX_WAVE_SIZE 9
#endif
#define CWT_CHANNELS (2 * WAVE_COUNT)

#define MIN_CDS 30
#define DEFAULT_ITERATIONS 12
#define DEFAULT_THREADS 1
#define DEFAULT_BLOCK 65536

typedef struct {
   double _Complex *kernel[WAVE_COUNT];
   int widths[WAVE_COUNT], max_width;
} Wavelets;
typedef struct {
   char *name, *seq;
   int length;
   float *cwt_power;
   uint8_t *label;
} Contig;

extern const int wave_sizes[WAVE_COUNT];
void anno_fail(const char *format, ...);
void *anno_alloc(size_t count, size_t size);
int base_index(char base);
void wavelets_init(Wavelets *wavelets);
void cwt_extract(const Wavelets *wavelets, const char *sequence, int length,
                 int start, int count, double *features);
Contig *read_fasta(const char *path, int *count);
int label_cds(Contig *contigs, int count, const char *gff, int *skipped);
void write_gff(const Contig *contig, const uint8_t *path, unsigned long *genes);

#endif
