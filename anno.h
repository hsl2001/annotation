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

#define GFF_BOUND_KINDS 4
#define GFF_CWT_RADIUS 26
#define GFF_CWT_WINDOW (2 * GFF_CWT_RADIUS + 1)
#define GFF_CWT_BASIS (GFF_BOUND_KINDS * GFF_CWT_WINDOW)
#define GFF_CWT_RIDGE 1.0

#define MIN_CDS 6
#define DEFAULT_ITERATIONS 10

typedef struct {
  double _Complex *kernel[WAVE_COUNT];
  int widths[WAVE_COUNT], max_width;
} Wavelets;

typedef struct {
  float bound[GFF_BOUND_KINDS][GFF_CWT_WINDOW][CWT_CHANNELS];
} GffCwtModel;

typedef struct {
  double *gram[2];
  double *rhs;
  uint64_t valid_bases;
} GffCwtAccumulator;

typedef struct {
  char *name, *seq;
  int length;
  float *cwt;
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

static inline int is_cds(int state) {
  return (state >= 1 && state <= 3) || (state >= 7 && state <= 9);
}

static inline int is_intron(int state) {
  return (state >= 4 && state <= 6) || (state >= 10 && state <= 12);
}

/* Frame-agnostic region class: 0=intergenic, 1=CDS, 2=intron (no codon periodicity). */
static inline int emit_region(int state) {
  return state <= 0 ? 0 : is_cds(state) ? 1 : 2;
}

static inline int boundary_kind(int previous, int current) {
  if ((previous >= 1 && previous <= 3 && current >= 4 && current <= 6) ||
      (previous >= 10 && current >= 7 && current <= 9)) return 0;
  if ((previous >= 4 && previous <= 6 && current >= 1 && current <= 3) ||
      (previous >= 7 && previous <= 9 && current >= 10)) return 1;
  if ((previous == 0 && current == 1) || (previous == 7 && current == 0)) return 2;
  if ((previous == 3 && current == 0) || (previous == 0 && current == 9)) return 3;
  return -1;
}

static inline int boundary_rev(int previous, int current) {
  return previous >= 7 || current >= 7;
}

static inline int gff_cwt_bound_pos(int boundary_pos, int rev, int slot) {
  int delta = slot - GFF_CWT_RADIUS;
  return rev ? boundary_pos - delta - 1 : boundary_pos + delta;
}

static inline double cwt_strand_sign(int rev, int ch) {
  return (rev && (ch & 1)) ? -1.0 : 1.0;
}

static inline double gff_cwt_bound_energy(const GffCwtModel *model, int kind) {
  double sum = 0.0;
  for (int slot = 0; slot < GFF_CWT_WINDOW; slot++)
    for (int ch = 0; ch < CWT_CHANNELS; ch++) {
      double v = (double)model->bound[kind][slot][ch];
      sum += v * v;
    }
  return sum;
}

Contig *read_fasta(const char *path, int *count);
int label_cds(Contig *contigs, int count, const char *gff, int *skipped);
void write_gff(const Contig *contig, const uint8_t *path, unsigned long *genes);

void gff_cwt_accumulator_init(GffCwtAccumulator *accumulator);
void gff_cwt_accumulate_range(GffCwtAccumulator *accumulator, const Contig *contig,
                              int start, int length, const float *cwt);
void gff_cwt_accumulator_finish(GffCwtAccumulator *accumulator, GffCwtModel *model);
void gff_cwt_predict_soft(const double *bound_marg, int length,
                          const GffCwtModel *model, double *out_cwt);
void gff_cwt_residual_adjoint(const double *residual, int length,
                              const GffCwtModel *model, double *grad_bound);

#endif
