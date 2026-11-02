#ifndef ANNO_RFFT_H
#define ANNO_RFFT_H

#include <complex.h>
#include <stdbool.h>
#include <stddef.h>

static void anno_fft_radix2(double complex *vector, size_t size, bool inverse) {
  for (size_t i = 1, j = 0; i < size; i++) {
    size_t bit = size >> 1;
    for (; j & bit; bit >>= 1)
      j ^= bit;
    j ^= bit;
    if (i < j) {
      double complex temporary = vector[i];
      vector[i] = vector[j];
      vector[j] = temporary;
    }
  }
  for (size_t half = 1; half < size; half <<= 1) {
    double angle = (inverse ? 1.0 : -1.0) * acos(-1.0) / (double)half;
    double complex step = cos(angle) + I * sin(angle);
    for (size_t start = 0; start < size; start += 2 * half) {
      double complex twiddle = 1.0;
      for (size_t offset = 0; offset < half; offset++) {
        double complex even = vector[start + offset];
        double complex odd = vector[start + offset + half] * twiddle;
        vector[start + offset] = even + odd;
        vector[start + offset + half] = even - odd;
        twiddle *= step;
      }
    }
  }
}

#endif
