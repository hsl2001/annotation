#ifndef RFFT_H
#define RFFT_H

#include <complex.h>
#include <stdbool.h>
#include <stddef.h>

void fft_transform_radix2_twiddle(double _Complex *restrict vec, size_t n, bool inverse,
                                  const double _Complex *restrict twiddle);

#endif

#ifdef RFFT_IMPLEMENTATION

void fft_transform_radix2_twiddle(double complex *restrict vec, size_t n, bool inverse,
                                  const double complex *restrict twiddle) {
  int levels = 0;
  for (size_t k = 1; (k & n) == 0; k <<= 1) levels++;

  for (size_t i = 0; i < n; i++) {
    size_t j = 0, ii = i;
    for (int k = 0; k < levels; k++) {
      j = (j << 1) | (ii & 1);
      ii >>= 1;
    }
    if (j > i) {
      double complex tmp = vec[i];
      vec[i] = vec[j];
      vec[j] = tmp;
    }
  }

  for (size_t half = 1; half < n; half *= 2) {
    size_t size = 2 * half, step = n / size;
    for (size_t i = 0; i < n; i += size)
      for (size_t j = 0, k = 0; j < half; j++, k += step) {
        double complex w = inverse ? conj(twiddle[k]) : twiddle[k];
        double complex tmp = vec[i + j + half] * w;
        vec[i + j + half] = vec[i + j] - tmp;
        vec[i + j] += tmp;
      }
  }
}

#endif