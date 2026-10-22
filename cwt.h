#ifndef ANNO_CWT_H
#define ANNO_CWT_H

#include "anno.h"

#endif

#if defined(CWT_IMPLEMENTATION) && !defined(ANNO_CWT_IMPLEMENTATION)
#define ANNO_CWT_IMPLEMENTATION

#include <complex.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

_Static_assert(sizeof(float) == 4 && FLT_RADIX == 2 && FLT_MANT_DIG == 24,
               "CWT requires binary32 floats");

static double complex cwt_base_signal(char base) {
  const double complex mapping[4] = {1.0 + I, -1.0 + I, -1.0 - I, 1.0 - I};
  int index = base_index(base);
  return index < 0 ? 0.0 : mapping[index] / sqrt(2.0);
}

static double complex cwt_morlet_integral(double lower, double upper) {
  static const double nodes[8] = {
    0.09501250983763744, 0.2816035507792589, 0.4580167776572274, 0.6178762444026438,
    0.7554044083550030, 0.8656312023878318, 0.9445750230732326, 0.9894009349916499
  };
  static const double weights[8] = {
    0.1894506104550685, 0.1826034150449236, 0.1691565193950025, 0.1495959888165767,
    0.1246289712555339, 0.0951585116824928, 0.0622535239386479, 0.0271524594117541
  };
  lower = fmax(lower, -8.0);
  upper = fmin(upper, 8.0);
  if (lower >= upper) return 0.0;
  double normalization = 1.0 / sqrt(sqrt(acos(-1.0)) *
                                    (1.0 + exp(-36.0) - 2.0 * exp(-27.0)));
  int segments = (int)ceil((upper - lower) * 2.0);
  double step = (upper - lower) / segments;
  double complex integral = 0.0;
  for (int segment = 0; segment < segments; segment++) {
    double middle = lower + (segment + 0.5) * step;
    for (int node = 0; node < 8; node++) {
      double left = middle - nodes[node] * step / 2.0;
      double right = middle + nodes[node] * step / 2.0;
      integral += weights[node] * step / 2.0 * (
        exp(-0.5 * left * left) * (cexp(-I * 6.0 * left) - exp(-18.0)) +
        exp(-0.5 * right * right) * (cexp(-I * 6.0 * right) - exp(-18.0)));
    }
  }
  return normalization * integral;
}

static void wavelets_init(Wavelets *wavelets) {
  memset(wavelets, 0, sizeof(*wavelets));
  for (int scale = 0; scale < WAVE_COUNT; scale++) {
    double dilation = wave_sizes[scale];
    if (!isfinite(dilation) || dilation <= 0.0 || dilation != floor(dilation))
      anno_fail("CWT scales must be positive natural numbers");
    double radius_value = ceil(8.0 * dilation + 0.5);
    if (radius_value > (INT_MAX - 1) / 2.0) anno_fail("CWT scale exceeds supported kernel size");
    int radius = (int)radius_value;
    int width = 2 * radius + 1;
    wavelets->kernel[scale] = anno_alloc((size_t)width, sizeof(double complex));
    wavelets->widths[scale] = width;
    if (width > wavelets->max_width) wavelets->max_width = width;
    for (int tap = 0; tap < width; tap++) {
      double displacement = tap - radius;
      wavelets->kernel[scale][tap] = sqrt(dilation) * cwt_morlet_integral(
        (displacement - 0.5) / dilation, (displacement + 0.5) / dilation);
    }
  }
}

static void wavelets_free(Wavelets *wavelets) {
  for (int scale = 0; scale < WAVE_COUNT; scale++) free(wavelets->kernel[scale]);
}

static size_t convolution_fft_size(size_t signal_length, size_t kernel_length) {
  if (signal_length > SIZE_MAX - (kernel_length - 1)) anno_fail("CWT FFT size overflow");
  size_t required = signal_length + kernel_length - 1, size = 1;
  while (size < required) {
    if (size > SIZE_MAX / 2) anno_fail("CWT FFT size overflow");
    size <<= 1;
  }
  return size;
}

static void cwt_window_scale(
  const char *restrict sequence, size_t length, size_t window, int scale, int width,
  size_t max_fft, const double complex *restrict kernel_table,
  const double complex *restrict twiddle_table, double complex *restrict signal_fft,
  float *restrict complex_values
) {
  size_t start = window * (size_t)CWT_WINDOW_SIZE;
  size_t window_length = length - start;
  if (window_length > CWT_WINDOW_SIZE) window_length = CWT_WINDOW_SIZE;
  size_t end = start + window_length;
  size_t radius = (size_t)(width / 2);
  size_t input_start = start > radius ? start - radius : 0;
  size_t input_end = radius > length - end ? length : end + radius;
  size_t input_length = input_end - input_start;
  size_t required = input_length + (size_t)width - 1, fft_size = 1;
  int level = 0;
  while (fft_size < required) {
    fft_size <<= 1;
    level++;
  }

  const double complex *twiddle = twiddle_table + (fft_size >> 1);
  const double complex *kernel_fft =
    kernel_table + ((size_t)level * WAVE_COUNT + (size_t)scale) * max_fft;

  for (size_t i = 0; i < input_length; i++)
    signal_fft[i] = cwt_base_signal(sequence[input_start + i]);
  memset(signal_fft + input_length, 0, (fft_size - input_length) * sizeof(*signal_fft));

  fft_transform_radix2_twiddle(signal_fft, fft_size, false, twiddle);
  for (size_t i = 0; i < fft_size; i++) signal_fft[i] *= kernel_fft[i];
  fft_transform_radix2_twiddle(signal_fft, fft_size, true, twiddle);

  double inv_size = 1.0 / (double)fft_size;
  for (size_t position = start; position < end; position++) {
    size_t index = position - input_start + (size_t)(width - 1 - width / 2);
    double complex coefficient = signal_fft[index] * inv_size;
    size_t feature = position * CWT_CHANNELS + 2 * (size_t)scale;
    complex_values[feature] = (float)creal(coefficient);
    complex_values[feature + 1] = (float)cimag(coefficient);
  }
}

typedef struct {
  const Wavelets *wavelets;
  const char *sequence;
  size_t length;
  size_t window_count;
  size_t max_fft;
  const double complex *kernel_table;
  const double complex *twiddle_table;
  float *complex_values;
  atomic_size_t next_window;
} CwtTask;

static void *cwt_worker(void *arg) {
  CwtTask *task = arg;
  double complex *signal_fft = anno_alloc(task->max_fft, sizeof(*signal_fft));
  for (;;) {
    size_t window = atomic_fetch_add_explicit(&task->next_window, 1, memory_order_relaxed);
    if (window >= task->window_count) break;
    for (int scale = 0; scale < WAVE_COUNT; scale++)
      cwt_window_scale(task->sequence, task->length, window, scale,
                       task->wavelets->widths[scale], task->max_fft,
                       task->kernel_table, task->twiddle_table, signal_fft,
                       task->complex_values);
  }
  free(signal_fft);
  return NULL;
}

static void cwt_transform(const Wavelets *wavelets, const char *sequence, int length,
                          int threads, float *complex_values) {
  if (length <= 0) return;
  if (threads < 1) anno_fail("CWT requires at least one worker thread");
  if ((size_t)length > SIZE_MAX / (CWT_CHANNELS * sizeof(*complex_values)))
    anno_fail("CWT output size overflow");
  size_t sequence_length = (size_t)length;
  size_t window_count = sequence_length / CWT_WINDOW_SIZE +
                        (sequence_length % CWT_WINDOW_SIZE != 0);
  if (window_count > SIZE_MAX / WAVE_COUNT) anno_fail("CWT task count overflow");

  size_t max_input = CWT_WINDOW_SIZE + (size_t)wavelets->max_width;
  if (max_input > sequence_length) max_input = sequence_length;
  size_t max_fft = convolution_fft_size(max_input, (size_t)wavelets->max_width);
  int max_level = 0;
  for (size_t s = 1; s < max_fft; s <<= 1) max_level++;

  double complex *twiddle_table = anno_alloc(max_fft, sizeof(*twiddle_table));
  double complex *kernel_table =
    anno_alloc((size_t)(max_level + 1) * WAVE_COUNT * max_fft, sizeof(*kernel_table));
  for (int lvl = 1; lvl <= max_level; lvl++) {
    size_t fft_size = (size_t)1 << lvl;
    double complex *tw = twiddle_table + (fft_size >> 1);
    double base_angle = -2.0 * acos(-1.0) / (double)fft_size;
    for (size_t k = 0; k < fft_size / 2; k++) {
      double angle = base_angle * (double)k;
      tw[k] = cos(angle) + I * sin(angle);
    }
    for (int scale = 0; scale < WAVE_COUNT; scale++) {
      int width = wavelets->widths[scale];
      if ((size_t)width > fft_size) continue;
      double complex *kf = kernel_table + ((size_t)lvl * WAVE_COUNT + (size_t)scale) * max_fft;
      for (int tap = 0; tap < width; tap++)
        kf[width - 1 - tap] = wavelets->kernel[scale][tap];
      fft_transform_radix2_twiddle(kf, fft_size, false, tw);
    }
  }

  CwtTask task = {.wavelets = wavelets, .sequence = sequence, .length = sequence_length,
                  .window_count = window_count, .max_fft = max_fft,
                  .kernel_table = kernel_table, .twiddle_table = twiddle_table,
                  .complex_values = complex_values};
  atomic_init(&task.next_window, 0);
  int worker_count = window_count < (size_t)threads ? (int)window_count : threads;
  pthread_t *workers = worker_count > 1
      ? anno_alloc((size_t)(worker_count - 1), sizeof(*workers)) : NULL;
  for (int worker = 0; worker < worker_count - 1; worker++)
    if (pthread_create(workers + worker, NULL, cwt_worker, &task))
      anno_fail("Cannot spawn CWT thread");
  cwt_worker(&task);
  for (int worker = 0; worker < worker_count - 1; worker++)
    if (pthread_join(workers[worker], NULL)) anno_fail("Cannot join CWT thread");
  free(workers);
  free(kernel_table);
  free(twiddle_table);
}

static inline int gff_cwt_bound_basis(int kind, int slot) {
  return kind * GFF_CWT_WINDOW + slot;
}

void gff_cwt_predict_soft(const double *bound_marg, int length,
                          const GffCwtModel *model, double *out_cwt) {
  if (length <= 0) return;
  memset(out_cwt, 0, (size_t)length * CWT_CHANNELS * sizeof(*out_cwt));
  if (!bound_marg) return;
  for (int b = 1; b < length; b++) {
    const double *bm = bound_marg + (size_t)b * (2 * GFF_BOUND_KINDS);
    for (int rev = 0; rev < 2; rev++)
      for (int kind = 0; kind < GFF_BOUND_KINDS; kind++) {
        double p = bm[rev * GFF_BOUND_KINDS + kind];
        if (p == 0.0) continue;
        for (int slot = 0; slot < GFF_CWT_WINDOW; slot++) {
          int i = gff_cwt_bound_pos(b, rev, slot);
          if (i < 0 || i >= length) continue;
          double *dst = out_cwt + (size_t)i * CWT_CHANNELS;
          const float *w = model->bound[kind][slot];
          for (int ch = 0; ch < CWT_CHANNELS; ch++)
            dst[ch] += p * cwt_strand_sign(rev, ch) * (double)w[ch];
        }
      }
  }
}

void gff_cwt_residual_adjoint(const double *residual, int length,
                              const GffCwtModel *model, double *grad_bound) {
  if (length <= 0 || !grad_bound) return;
  memset(grad_bound, 0, (size_t)length * (2 * GFF_BOUND_KINDS) * sizeof(*grad_bound));
  for (int b = 1; b < length; b++) {
    double *gb = grad_bound + (size_t)b * (2 * GFF_BOUND_KINDS);
    for (int rev = 0; rev < 2; rev++)
      for (int kind = 0; kind < GFF_BOUND_KINDS; kind++) {
        double sum = 0.0;
        for (int slot = 0; slot < GFF_CWT_WINDOW; slot++) {
          int i = gff_cwt_bound_pos(b, rev, slot);
          if (i < 0 || i >= length) continue;
          const double *res = residual + (size_t)i * CWT_CHANNELS;
          const float *w = model->bound[kind][slot];
          for (int ch = 0; ch < CWT_CHANNELS; ch++)
            sum += res[ch] * cwt_strand_sign(rev, ch) * (double)w[ch];
        }
        gb[rev * GFF_BOUND_KINDS + kind] = sum;
      }
  }
}

void gff_cwt_accumulator_init(GffCwtAccumulator *accumulator) {
  memset(accumulator, 0, sizeof(*accumulator));
  for (int p = 0; p < 2; p++)
    accumulator->gram[p] = anno_alloc((size_t)GFF_CWT_BASIS * GFF_CWT_BASIS,
                                      sizeof(*accumulator->gram[p]));
  accumulator->rhs = anno_alloc((size_t)GFF_CWT_BASIS * CWT_CHANNELS,
                                sizeof(*accumulator->rhs));
}

void gff_cwt_accumulate_range(GffCwtAccumulator *accumulator, const Contig *contig,
                              int start, int length, const float *cwt) {
  if (!contig->label || !cwt || start < 0 || length < 0 ||
      start > contig->length - length)
    anno_fail("Invalid GFF->CWT training range");
  const uint8_t *labels = contig->label;
  const char *sequence = contig->seq + start;

  for (int i = 0; i < length; i++) {
    if (base_index(sequence[i]) < 0) continue;
    int pos = start + i;
    accumulator->valid_bases++;
    int active[2 * GFF_CWT_WINDOW + 4];
    int active_rev[2 * GFF_CWT_WINDOW + 4];
    int n_active = 0;

    int b_lo = pos - GFF_CWT_RADIUS - 1;
    if (b_lo < 1) b_lo = 1;
    int b_hi = pos + GFF_CWT_RADIUS + 1;
    if (b_hi >= contig->length) b_hi = contig->length - 1;
    for (int b = b_lo; b <= b_hi; b++) {
      int kind = boundary_kind(labels[b - 1], labels[b]);
      if (kind < 0) continue;
      int rev = boundary_rev(labels[b - 1], labels[b]);
      int delta = rev ? (b - 1 - pos) : (pos - b);
      if (delta < -GFF_CWT_RADIUS || delta > GFF_CWT_RADIUS) continue;
      active[n_active] = gff_cwt_bound_basis(kind, delta + GFF_CWT_RADIUS);
      active_rev[n_active++] = rev;
    }

    const float *obs = cwt + (size_t)i * CWT_CHANNELS;
    for (int a = 0; a < n_active; a++) {
      int ma = active[a], ra = active_rev[a];
      double *r_row = accumulator->rhs + (size_t)ma * CWT_CHANNELS;
      for (int ch = 0; ch < CWT_CHANNELS; ch++)
        r_row[ch] += cwt_strand_sign(ra, ch) * (double)obs[ch];
      double *g_re = accumulator->gram[0] + (size_t)ma * GFF_CWT_BASIS;
      double *g_im = accumulator->gram[1] + (size_t)ma * GFF_CWT_BASIS;
      for (int d = 0; d < n_active; d++) {
        g_re[active[d]] += 1.0;
        g_im[active[d]] += (ra == active_rev[d]) ? 1.0 : -1.0;
      }
    }
  }
}

void gff_cwt_accumulator_finish(GffCwtAccumulator *accumulator, GffCwtModel *model) {
  if (!accumulator->gram[0] || !accumulator->gram[1] ||
      !accumulator->rhs || !accumulator->valid_bases)
    anno_fail("GFF->CWT training requires at least one valid labeled base");
  double *rhs = accumulator->rhs;

  for (int p = 0; p < 2; p++) {
    double *gram = accumulator->gram[p];
    for (int m = 0; m < GFF_CWT_BASIS; m++)
      gram[(size_t)m * GFF_CWT_BASIS + m] += GFF_CWT_RIDGE;

    for (int i = 0; i < GFF_CWT_BASIS; i++)
      for (int j = 0; j <= i; j++) {
        double sum = gram[(size_t)i * GFF_CWT_BASIS + j];
        for (int k = 0; k < j; k++)
          sum -= gram[(size_t)i * GFF_CWT_BASIS + k] * gram[(size_t)j * GFF_CWT_BASIS + k];
        if (i == j) {
          if (!(sum > 0.0) || !isfinite(sum))
            anno_fail("Non-positive definite Gram matrix in GFF->CWT training");
          gram[(size_t)i * GFF_CWT_BASIS + j] = sqrt(sum);
        } else {
          gram[(size_t)i * GFF_CWT_BASIS + j] = sum / gram[(size_t)j * GFF_CWT_BASIS + j];
        }
      }

    for (int i = 0; i < GFF_CWT_BASIS; i++) {
      double inv_diag = 1.0 / gram[(size_t)i * GFF_CWT_BASIS + i];
      for (int ch = p; ch < CWT_CHANNELS; ch += 2) {
        double sum = rhs[(size_t)i * CWT_CHANNELS + ch];
        for (int k = 0; k < i; k++)
          sum -= gram[(size_t)i * GFF_CWT_BASIS + k] * rhs[(size_t)k * CWT_CHANNELS + ch];
        rhs[(size_t)i * CWT_CHANNELS + ch] = sum * inv_diag;
      }
    }

    for (int i = GFF_CWT_BASIS - 1; i >= 0; i--) {
      double inv_diag = 1.0 / gram[(size_t)i * GFF_CWT_BASIS + i];
      for (int ch = p; ch < CWT_CHANNELS; ch += 2) {
        double sum = rhs[(size_t)i * CWT_CHANNELS + ch];
        for (int k = i + 1; k < GFF_CWT_BASIS; k++)
          sum -= gram[(size_t)k * GFF_CWT_BASIS + i] * rhs[(size_t)k * CWT_CHANNELS + ch];
        rhs[(size_t)i * CWT_CHANNELS + ch] = sum * inv_diag;
      }
    }
    free(gram);
  }

  memset(model, 0, sizeof(*model));
  for (int kind = 0; kind < GFF_BOUND_KINDS; kind++)
    for (int slot = 0; slot < GFF_CWT_WINDOW; slot++) {
      int m = gff_cwt_bound_basis(kind, slot);
      for (int ch = 0; ch < CWT_CHANNELS; ch++)
        model->bound[kind][slot][ch] = (float)rhs[(size_t)m * CWT_CHANNELS + ch];
    }

  free(rhs);
  memset(accumulator, 0, sizeof(*accumulator));
}

#endif
