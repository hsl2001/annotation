#include "anno_signal.h"

#include <math.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const int wave_sizes[WAVE_COUNT] = {4, 5, 8};

double otsu_threshold(const float *values, int count) {
  double minimum = INFINITY, maximum = -INFINITY;
  for (int index = 0; index < count; index++) {
    if (!isfinite(values[index])) continue;
    minimum = fmin(minimum, values[index]);
    maximum = fmax(maximum, values[index]);
  }
  if (!(maximum > minimum)) return INFINITY;
  double histogram[256] = {0}, total = 0, sum = 0;
  double step = (maximum - minimum) / 256;
  for (int index = 0; index < count; index++) {
    if (!isfinite(values[index])) continue;
    int bin = (int)((values[index] - minimum) / step);
    if (bin > 255) bin = 255;
    histogram[bin]++;
    total++;
    sum += bin;
  }
  double weight = 0, partial = 0, best = -1, threshold = INFINITY;
  for (int bin = 0; bin < 255; bin++) {
    weight += histogram[bin];
    partial += bin * histogram[bin];
    if (weight == 0 || weight == total) continue;
    double difference = partial / weight - (sum - partial) / (total - weight);
    double variance = weight * (total - weight) * difference * difference;
    if (variance > best) {
      best = variance;
      threshold = minimum + (bin + 1) * step;
    }
  }
  return threshold;
}

void anno_fail(const char *format, ...) {
  va_list arguments;
  va_start(arguments, format);
  fputs("Error: ", stderr);
  vfprintf(stderr, format, arguments);
  fputc('\n', stderr);
  va_end(arguments);
  exit(EXIT_FAILURE);
}

void *anno_alloc(size_t count, size_t size) {
  if (size && count > SIZE_MAX / size)
    anno_fail("Allocation size overflow");
  void *memory = calloc(count ? count : 1, size);
  if (!memory)
    anno_fail("Out of memory");
  return memory;
}

int base_index(char base) {
  switch (base) {
  case 'A': case 'a': return 0;
  case 'C': case 'c': return 1;
  case 'G': case 'g': return 2;
  case 'T': case 't': return 3;
  default: return -1;
  }
}

double complex base_signal(char base) {
  const double complex mapping[4] = {1.0, I, -I, -1.0};
  int index = base_index(base);
  return index < 0 ? 0.0 : mapping[index];
}

char *reverse_complement(const char *sequence, int length) {
  char *result = anno_alloc((size_t)length + 1, sizeof(char));
  for (int position = 0; position < length; position++) {
    int base = base_index(sequence[length - position - 1]);
    result[position] = base < 0 ? 'N' : "TGCA"[base];
  }
  return result;
}

void wavelets_init(Wavelets *wavelets) {
  memset(wavelets, 0, sizeof(*wavelets));
  for (int scale = 0; scale < WAVE_COUNT; scale++) {
    int width = wave_sizes[scale];
    double complex mean = 0.0;
    for (int tap = 0; tap < width; tap++) {
      double coordinate = (tap - width / 2.0 + 0.5) / (width / 8.0);
      double complex value = exp(-0.5 * coordinate * coordinate) *
                             (cexp(I * 6.0 * coordinate) - exp(-18.0));
      wavelets->kernel[scale][tap] = value;
      mean += value / width;
    }
    double energy = 0.0;
    for (int tap = 0; tap < width; tap++) {
      wavelets->kernel[scale][tap] -= mean;
      double magnitude = cabs(wavelets->kernel[scale][tap]);
      energy += magnitude * magnitude;
    }
    for (int tap = 0; tap < width; tap++)
      wavelets->kernel[scale][tap] =
          conj(wavelets->kernel[scale][tap]) / sqrt(energy);
  }
}

void cwt_extract(const Wavelets *wavelets, const char *sequence, int length,
                 int start, int count, double *features) {
  memset(features, 0, (size_t)count * CWT_CHANNELS * sizeof(double));
  for (int offset = 0; offset < count; offset++) {
    long position = (long)start + offset;
    if (position < 0 || position >= length)
      continue;
    for (int scale = 0; scale < WAVE_COUNT; scale++) {
      int width = wave_sizes[scale];
      double complex coefficient = 0.0;
      for (int tap = 0; tap < width; tap++) {
        long context = position + tap - width / 2;
        if (context >= 0 && context < length)
          coefficient += base_signal(sequence[context]) *
                         wavelets->kernel[scale][tap];
      }
      features[(size_t)offset * CWT_CHANNELS + 2 * scale] = creal(coefficient);
      features[(size_t)offset * CWT_CHANNELS + 2 * scale + 1] = cimag(coefficient);
    }
  }
}