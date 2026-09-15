#ifndef ANNO_H
#define ANNO_H

#include <stddef.h>

#ifdef __cplusplus
#include <complex>
typedef std::complex<double> anno_complex;
#else
#include <complex.h>
typedef double complex anno_complex;
#endif

#define WAVE_COUNT 3
#define CWT_CHANNELS (2 * WAVE_COUNT)
#define MAX_WAVE_SIZE 8

typedef struct {
  anno_complex kernel[WAVE_COUNT][MAX_WAVE_SIZE];
} Wavelets;

extern const int wave_sizes[WAVE_COUNT];
void anno_fail(const char *format, ...);
void *anno_alloc(size_t count, size_t size);
int base_index(char base);
anno_complex base_signal(char base);
char *reverse_complement(const char *sequence, int length);
void wavelets_init(Wavelets *wavelets);
void cwt_extract(const Wavelets *wavelets, const char *sequence, int length,
                 int start, int count, double *features);

/* Boundary-peak exon caller tunables. */
#define POWER_THRESHOLD 1.5
#define EXON_WINDOW 200
#define MIN_EXON_PEAKS 3

/* Predicted exon interval, 0-based half-open [start, end). */
typedef struct {
  int start;
  int end;
  float score;
} Exon;

/* Detect threshold peaks independently in each scale, keep only positions that
   peak in every scale, group those common peaks that fall within EXON_WINDOW
   into exons, and keep clusters holding at least MIN_EXON_PEAKS peaks. `powers`
   is length*wave_count row-major per-scale energy; `exons` must have capacity
   for `length` entries. Returns exon count. */
int anno_call_exons(const float *powers, int length, int wave_count, Exon *exons);

#endif
