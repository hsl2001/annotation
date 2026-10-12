#ifndef ANNO_CWT_H
#define ANNO_CWT_H

#include "anno.h"

#include <errno.h>
#include <complex.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <unistd.h>

_Static_assert(sizeof(float) == 4 && FLT_RADIX == 2 && FLT_MANT_DIG == 24,
               "CWT export requires binary32 floats");

#define CWT_CHANNELS (2 * WAVE_COUNT)

static const int wave_sizes[WAVE_COUNT] = {CWT_SCALES};

typedef struct {
  double complex *kernel[WAVE_COUNT];
  int widths[WAVE_COUNT];
} Wavelets;

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
  const double omega_squared = MORLET_OMEGA0 * MORLET_OMEGA0;
  const double correction = exp(-0.5 * omega_squared);
  lower = fmax(lower, -MORLET_CUTOFF);
  upper = fmin(upper, MORLET_CUTOFF);
  if (lower >= upper) return 0.0;
  double normalization = 1.0 / sqrt(sqrt(acos(-1.0)) *
                                   (1.0 + exp(-omega_squared) - 2.0 * exp(-0.75 * omega_squared)));
  int segments = (int)ceil((upper - lower) * 2.0);
  double step = (upper - lower) / segments;
  double complex integral = 0.0;
  for (int segment = 0; segment < segments; segment++) {
    double middle = lower + (segment + 0.5) * step;
    for (int node = 0; node < 8; node++) {
      double left = middle - nodes[node] * step / 2.0;
      double right = middle + nodes[node] * step / 2.0;
      integral += weights[node] * step / 2.0 * (
        exp(-0.5 * left * left) * (cexp(-I * MORLET_OMEGA0 * left) - correction) +
        exp(-0.5 * right * right) * (cexp(-I * MORLET_OMEGA0 * right) - correction));
    }
  }
  return normalization * integral;
}

static void wavelets_init(Wavelets *wavelets) {
  for (int scale = 0; scale < WAVE_COUNT; scale++) {
    double dilation = wave_sizes[scale];
    if (!isfinite(dilation) || dilation <= 0.0 || dilation != floor(dilation))
      anno_fail("CWT scales must be positive natural numbers");
    double radius_value = ceil(MORLET_CUTOFF * dilation + 0.5);
    if (radius_value > (INT_MAX - 1) / 2.0) anno_fail("CWT scale exceeds supported kernel size");
    int radius = (int)radius_value;
    int width = 2 * radius + 1;
    wavelets->kernel[scale] = anno_alloc((size_t)width, sizeof(double complex));
    wavelets->widths[scale] = width;
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

/* Each CWT_WINDOW_SIZE block is convolved with its kernel halo by FFT; outputs may be NULL. */
static void cwt_transform(const Wavelets *wavelets, const char *sequence, int length,
                          float *complex_values, float *powers) {
  double complex *signal_fft = NULL, *kernel_fft = NULL;
  size_t capacity = 0, total = (size_t)length;
  for (size_t start = 0; start < total; start += CWT_WINDOW_SIZE) {
    size_t end = total - start > CWT_WINDOW_SIZE ? start + CWT_WINDOW_SIZE : total;
    for (int scale = 0; scale < WAVE_COUNT; scale++) {
      int width = wavelets->widths[scale];
      size_t radius = (size_t)(width / 2);
      size_t input_start = start > radius ? start - radius : 0;
      size_t input_end = radius > total - end ? total : end + radius;
      size_t input_length = input_end - input_start;
      size_t fft_size = convolution_fft_size(input_length, (size_t)width);
      if (fft_size > capacity) {
        free(signal_fft);
        free(kernel_fft);
        signal_fft = anno_alloc(fft_size, sizeof(*signal_fft));
        kernel_fft = anno_alloc(fft_size, sizeof(*kernel_fft));
        capacity = fft_size;
      }
      memset(signal_fft, 0, fft_size * sizeof(*signal_fft));
      memset(kernel_fft, 0, fft_size * sizeof(*kernel_fft));
      for (size_t position = input_start; position < input_end; position++)
        signal_fft[position - input_start] = cwt_base_signal(sequence[position]);
      for (int tap = 0; tap < width; tap++)
        kernel_fft[width - 1 - tap] = wavelets->kernel[scale][tap];
      fft_transform(signal_fft, fft_size, false);
      fft_transform(kernel_fft, fft_size, false);
      for (size_t index = 0; index < fft_size; index++)
        signal_fft[index] *= kernel_fft[index];
      fft_transform(signal_fft, fft_size, true);

      for (size_t position = start; position < end; position++) {
        size_t index = position - input_start + (size_t)(width - 1 - width / 2);
        double complex coefficient = signal_fft[index] / (double)fft_size;
        double real = creal(coefficient), imaginary = cimag(coefficient);
        size_t feature = position * CWT_CHANNELS + 2 * (size_t)scale;
        if (complex_values)
          complex_values[feature] = (float)real, complex_values[feature + 1] = (float)imaginary;
        if (powers)
          powers[position * WAVE_COUNT + (size_t)scale] =
            (float)(real * real + imaginary * imaginary);
      }
    }
  }
  free(signal_fft);
  free(kernel_fft);
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
