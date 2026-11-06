/**
 * @file cwt.c
 * @brief CWT transforms and per-base model features.
 */

#include "include/cwt.h"
#include "include/surf.h"
#include "include/utils.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/** Morlet scales of the CWT channels. */
static const int cwt_scales[CWT_CHANNELS] = {3, 4, 5};

/**
 * @brief Transform a complex vector with an in-place radix-2 FFT.
 * @param[in,out] vector: Complex input vector; replaced by its transform.
 * @param[in] size: Power-of-two vector length.
 * @param[in] inverse Nonzero selects the inverse transform.
 * @return None.
 */
static void fft_radix2(double complex *vector, size_t size, int inverse) {
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
    double angle = (inverse == 1 ? 1.0 : -1.0) * acos(-1.0) / (double)half;
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

static void transform(const char *bases, int64_t length, int64_t start,
                      int64_t end, float *out, int power) {
  static const double complex signal[NUCLEOTIDES] = {1.0 + I, -1.0 + I,
                                                     -1.0 - I, 1.0 - I};
  static const double nodes[8] = {0.09501250983763744, 0.2816035507792589,
                                  0.4580167776572274,  0.6178762444026438,
                                  0.7554044083550030,  0.8656312023878318,
                                  0.9445750230732324,  0.9894009349916499};
  static const double weights[8] = {0.1894506104550685, 0.1826034150449236,
                                    0.1691565193950025, 0.1495959888165767,
                                    0.1246289712555339, 0.0951585116824928,
                                    0.0622535239386479, 0.0271524594117541};
  if (start >= end || start < 0 || end > length)
    return;
  int width[CWT_CHANNELS];
  const int max_width = 4 * CWT_RADIUS * cwt_scales[CWT_CHANNELS - 1] + 3;
  const double support = 2.0 * CWT_RADIUS;
  const double omega2 = MORLET_CENTRAL * MORLET_CENTRAL;
  const double correction = exp(-0.5 * omega2);
  const double normalization =
      1.0 /
      sqrt(sqrt(acos(-1.0)) * (1.0 + exp(-omega2) - 2.0 * exp(-0.75 * omega2)));
  double complex kernel[CWT_CHANNELS][max_width];
  for (int s = 0; s < CWT_CHANNELS; s++) {
    double scale = cwt_scales[s];
    width[s] = 4 * CWT_RADIUS * cwt_scales[s] + 3;
    int radius = width[s] / 2;
    for (int tap = 0; tap < width[s]; tap++) {
      double lower = fmax(((double)(tap - radius) - 0.5) / scale, -support);
      double upper = fmin(((double)(tap - radius) + 0.5) / scale, support);
      double complex integral = 0.0;
      /* Every scale is at least 3, so each bin needs only one segment. */
      if (lower < upper) {
        double half = (upper - lower) / 2.0, middle = (upper + lower) / 2.0;
        for (int node = 0; node < 8; node++) {
          double left = middle - nodes[node] * half;
          double right = middle + nodes[node] * half;
          integral += weights[node] * half *
                      (exp(-0.5 * left * left) *
                           (cexp(-I * MORLET_CENTRAL * left) - correction) +
                       exp(-0.5 * right * right) *
                           (cexp(-I * MORLET_CENTRAL * right) - correction));
        }
      }
      kernel[s][tap] = sqrt(scale) * normalization * integral;
    }
  }
  size_t input_length = (size_t)(end - start);
  if (input_length > CWT_WINDOW)
    input_length = CWT_WINDOW;
  input_length += max_width;
  size_t fft_size = 1;
  while (fft_size < input_length + max_width - 1)
    fft_size <<= 1;
  double complex *kernel_fft =
      scalloc(smul(CWT_CHANNELS, fft_size), sizeof(*kernel_fft));
  for (int s = 0; s < CWT_CHANNELS; s++) {
    for (int tap = 0; tap < width[s]; tap++)
      kernel_fft[(size_t)s * fft_size + width[s] - 1 - tap] = kernel[s][tap];
    fft_radix2(kernel_fft + (size_t)s * fft_size, fft_size, 0);
  }
  const int64_t window_count = (end - start + CWT_WINDOW - 1) / CWT_WINDOW;
#pragma omp parallel for schedule(static)
  for (int64_t window = 0; window < window_count; window++) {
    int64_t window_start = start + window * CWT_WINDOW;
    int64_t window_end = window_start + CWT_WINDOW;
    if (window_end > end)
      window_end = end;
    int radius = max_width / 2;
    int64_t input_start = window_start > radius ? window_start - radius : 0;
    int64_t input_end =
        window_end + radius < length ? window_end + radius : length;
    size_t input_count = (size_t)(input_end - input_start);
    double complex *transformed =
        scalloc(smul(fft_size, 2), sizeof(*transformed));
    double complex *spectrum = transformed + fft_size;
    for (size_t i = 0; i < input_count; i++) {
      int base = base_index(bases[input_start + (int64_t)i]);
      spectrum[i] = base < 0 ? 0.0 : signal[base] / sqrt(2.0);
    }
    fft_radix2(spectrum, fft_size, 0);
    for (int s = 0; s < CWT_CHANNELS; s++) {
      for (size_t i = 0; i < fft_size; i++)
        transformed[i] = spectrum[i] * kernel_fft[(size_t)s * fft_size + i];
      fft_radix2(transformed, fft_size, 1);
      for (int64_t position = window_start; position < window_end; position++) {
        size_t index = (size_t)(position - input_start) +
                       (size_t)(width[s] - 1 - width[s] / 2);
        double complex coefficient = transformed[index] / (double)fft_size;
        float real = (float)creal(coefficient),
              imag = (float)cimag(coefficient);
        size_t out_index =
            (size_t)(position - start) * CWT_CHANNELS + (size_t)s;
        if (power)
          out[out_index] = real * real + imag * imag;
        else {
          out[2 * out_index] = real;
          out[2 * out_index + 1] = imag;
        }
      }
    }
    free(transformed);
  }
  free(kernel_fft);
}

void cwt(const char *bases, int64_t length, int64_t start, int64_t end,
         float *out) {
  transform(bases, length, start, end, out, 0);
}

void cwt_power(const char *bases, int64_t length, int64_t start, int64_t end,
               float *out) {
  transform(bases, length, start, end, out, 1);
}

void compute_features(char base, const float *power, float *features) {
  memset(features, 0, NUCLEOTIDES * sizeof(*features));
  int index = base_index(base);
  if (index >= 0)
    features[index] = 1.0f;
  memcpy(features + NUCLEOTIDES, power, CWT_CHANNELS * sizeof(*power));
}
