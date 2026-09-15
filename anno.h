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

#define WAVE_COUNT 4
#define CWT_CHANNELS (2 * WAVE_COUNT)
#define MAX_WAVE_SIZE 9

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

#define POWER_THRESHOLD 0.5
#define EXON_WINDOW 400
#define MIN_EXON_PEAKS 3

typedef struct {
  float power_threshold;
  int exon_window;
  int min_exon_peaks;
} AnnoParameters;

typedef struct {
  int start;
  int end;
  float score;
} Exon;

int anno_call_exons(const float *powers, int length, int wave_count, Exon *exons);
int anno_call_exons_with_params(const float *powers, int length, int wave_count,
                                const AnnoParameters *parameters, Exon *exons);

#endif
