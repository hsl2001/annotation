#ifndef ANNO_H
#define ANNO_H

#include <complex.h>
#include <stddef.h>
#include <stdint.h>

#define WAVE_COUNT 4
#define MAX_WAVE_SIZE 9
#define CWT_CHANNELS (2 * WAVE_COUNT)

/* Linear-chain CRF decoded jointly for both strands.
   0 intergenic; 1-3 CDS(+) codon position; 4-6 intron(+) carrying the codon position of the
   next exon base; 7-9 CDS(-) codon position; 10-12 intron(-). */
#define STATES 13
#define WINDOW 65536
#define MIN_CDS 90
#define LEARNING_RATE 0.1f
#define L2_PENALTY 1e-3f

typedef struct { double complex kernel[WAVE_COUNT][MAX_WAVE_SIZE]; } Wavelets;
typedef struct { char *name, *seq; int length; uint16_t *cwt; uint8_t *label; } Contig;

extern const int wave_sizes[WAVE_COUNT];
void anno_fail(const char *format, ...);
void *anno_alloc(size_t count, size_t size);
int base_index(char base);
void wavelets_init(Wavelets *wavelets);
void cwt_extract(const Wavelets *wavelets, const char *sequence, int length,
                 int start, int count, double *features);
uint16_t *cwt_features(const Wavelets *wavelets, const char *sequence, int length);
Contig *read_fasta(const char *path, int *count);
int label_cds(Contig *contigs, int count, const char *gff, int *skipped);
size_t crf_weights(void);
void crf_train(Contig *contigs, int count, float *weights, int epochs);
void crf_decode(const float *weights, const Contig *contig, uint8_t *path);
void write_gff(const Contig *contig, const uint8_t *path, unsigned long *genes);

#endif
