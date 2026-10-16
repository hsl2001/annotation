#ifndef ANNO_CWT_H
#define ANNO_CWT_H

#include "anno.h"

#endif

#ifdef CWT_IMPLEMENTATION

#include <errno.h>
#include <complex.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <unistd.h>

_Static_assert(sizeof(float) == 4 && FLT_RADIX == 2 && FLT_MANT_DIG == 24,
               "CWT export requires binary32 floats");

static double complex cwt_base_signal(char base) {
  const double complex mapping[4] = {1.0, I, -I, -1.0};
  int index = base_index(base);
  return index < 0 ? 0.0 : mapping[index];
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
    size_t capacity = MAX_WAVE_SIZE > 0 ? (size_t)MAX_WAVE_SIZE : 1;
    if (capacity > (size_t)width) capacity = (size_t)width;
    while (capacity < (size_t)width) capacity *= 2;
    wavelets->kernel[scale] = anno_alloc(capacity, sizeof(double complex));
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
  double *restrict features, float *restrict complex_values, float *restrict powers
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
    double real = creal(coefficient), imaginary = cimag(coefficient);
    size_t feature = position * CWT_CHANNELS + 2 * (size_t)scale;
    if (features)
      features[feature] = real, features[feature + 1] = imaginary;
    if (complex_values)
      complex_values[feature] = (float)real,
      complex_values[feature + 1] = (float)imaginary;
    if (powers)
      powers[position * WAVE_COUNT + (size_t)scale] =
        (float)(real * real + imaginary * imaginary);
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
  double *features;
  float *complex_values, *powers;
  atomic_size_t next_window;
} CwtTask;

static void *cwt_worker(void *arg) {
  CwtTask *task = arg;
  double complex *signal_fft = anno_alloc(task->max_fft, sizeof(*signal_fft));
  size_t window = atomic_fetch_add_explicit(&task->next_window, 1, memory_order_relaxed);
  while (window < task->window_count) {
    for (int scale = 0; scale < WAVE_COUNT; scale++) {
      cwt_window_scale(task->sequence, task->length, window, scale,
                       task->wavelets->widths[scale], task->max_fft,
                       task->kernel_table, task->twiddle_table, signal_fft,
                       task->features, task->complex_values, task->powers);
    }
    window = atomic_fetch_add_explicit(&task->next_window, 1, memory_order_relaxed);
  }
  free(signal_fft);
  return NULL;
}

static void cwt_transform(const Wavelets *wavelets, const char *sequence, int length,
                          int threads, double *features, float *complex_values, float *powers) {
  if (length <= 0) return;
  if (threads < 1) anno_fail("CWT requires at least one worker thread");
  if (features && (size_t)length > SIZE_MAX / (CWT_CHANNELS * sizeof(*features)))
    anno_fail("CWT feature size overflow");
  if (complex_values && (size_t)length > SIZE_MAX / (CWT_CHANNELS * sizeof(*complex_values)))
    anno_fail("CWT output size overflow");
  if (powers && (size_t)length > SIZE_MAX / (WAVE_COUNT * sizeof(*powers)))
    anno_fail("CWT power size overflow");
  if (features) memset(features, 0, (size_t)length * CWT_CHANNELS * sizeof(*features));
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
                  .features = features, .complex_values = complex_values, .powers = powers};
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

static void *cwt_map_file(FILE *stream, size_t bytes, off_t start, void **mapping,
                          size_t *mapping_size, const char *kind) {
  if (bytes > (size_t)INT64_MAX || start < 0 ||
      (uintmax_t)start > (uintmax_t)INT64_MAX - bytes)
    anno_fail("CWT %s file size overflow", kind);
  off_t end = start + (off_t)bytes;
  if ((uintmax_t)end != (uintmax_t)start + bytes)
    anno_fail("CWT %s file size overflow", kind);
  long page_size = sysconf(_SC_PAGE_SIZE);
  if (page_size < 1) anno_fail("Cannot determine memory page size");
  off_t aligned = start - start % page_size;
  size_t prefix = (size_t)(start - aligned);
  if (bytes > SIZE_MAX - prefix) anno_fail("CWT mapping size overflow");
  size_t length = prefix + bytes;
  if (ftruncate(fileno(stream), end)) anno_fail("Cannot size %s file: %s", kind, strerror(errno));
  void *base = mmap(NULL, length, PROT_READ | PROT_WRITE, MAP_SHARED, fileno(stream), aligned);
  if (base == MAP_FAILED) anno_fail("Cannot map %s file: %s", kind, strerror(errno));
  *mapping = base;
  *mapping_size = length;
  return (unsigned char *)base + prefix;
}

static float *cwt_map_output(FILE *stream, size_t count, void **mapping,
                             size_t *mapping_size, off_t *position) {
  if (count > SIZE_MAX / sizeof(float)) anno_fail("CWT output size overflow");
  *position = ftello(stream);
  return cwt_map_file(stream, count * sizeof(float), *position, mapping, mapping_size,
                      "CWT output");
}

static void cwt_unmap_output(FILE *stream, float *values, size_t count, void *mapping,
                             size_t mapping_size, off_t position) {
  uint16_t endian = 1;
  if (!*(unsigned char *)&endian) {
    unsigned char *bytes = (unsigned char *)values;
    for (size_t i = 0; i < count; i++) {
      unsigned char first = bytes[4 * i], second = bytes[4 * i + 1];
      bytes[4 * i] = bytes[4 * i + 3], bytes[4 * i + 1] = bytes[4 * i + 2];
      bytes[4 * i + 2] = second, bytes[4 * i + 3] = first;
    }
  }
  if (msync(mapping, mapping_size, MS_SYNC)) anno_fail("Cannot flush CWT output: %s", strerror(errno));
  if (munmap(mapping, mapping_size)) anno_fail("Cannot unmap CWT output: %s", strerror(errno));
  off_t end = position + (off_t)(count * sizeof(float));
  if (fseeko(stream, end, SEEK_SET)) anno_fail("Cannot seek CWT output: %s", strerror(errno));
}

#endif
