/**
 * @file cwt.c
 * @brief CWT transforms and normalized per-base model features.
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

void cwt(const char *bases, int64_t length, int64_t start, int64_t end,
              float *out) {
  static const double complex signal[NUCLEOTIDES] = {
      1.0 + I, -1.0 + I, -1.0 - I, 1.0 - I};
  static const double nodes[8] = {
      0.09501250983763744, 0.2816035507792589, 0.4580167776572274,
      0.6178762444026438, 0.7554044083550030, 0.8656312023878318,
      0.9445750230732324, 0.9894009349916499};
  static const double weights[8] = {
      0.1894506104550685, 0.1826034150449236, 0.1691565193950025,
      0.1495959888165767, 0.1246289712555339, 0.0951585116824928,
      0.0622535239386479, 0.0271524594117541};
  if (start >= end || start < 0 || end > length)
    return;
  const int width[CWT_CHANNELS] = {67, 83, 99};
  const double support = 2.0 * CWT_RADIUS;
  double complex kernel[CWT_CHANNELS][99];
  for (int s = 0; s < CWT_CHANNELS; s++) {
    double scale = cwt_scales[s];
    int radius = width[s] / 2;
    for (int tap = 0; tap < width[s]; tap++) {
      double lower =
          fmax(((double)(tap - radius) - 0.5) / scale, -support);
      double upper =
          fmin(((double)(tap - radius) + 0.5) / scale, support);
      int segments = (int)ceil((upper - lower) * 2.0);
      double complex integral = 0.0;
      if (lower < upper) {
        double step = (upper - lower) / segments;
        for (int segment = 0; segment < segments; segment++) {
          double middle = lower + (segment + 0.5) * step;
          for (int node = 0; node < 8; node++) {
            double left = middle - nodes[node] * step / 2.0;
            double right = middle + nodes[node] * step / 2.0;
            integral += weights[node] * step / 2.0 *
                        (cexp(-0.5 * left * left) *
                             (cexp(-I * MORLET_CENTRAL * left) -
                              exp(-0.5 * MORLET_CENTRAL *
                                   MORLET_CENTRAL)) +
                         cexp(-0.5 * right * right) *
                             (cexp(-I * MORLET_CENTRAL * right) -
                              exp(-0.5 * MORLET_CENTRAL *
                                   MORLET_CENTRAL)));
          }
        }
      }
      double omega2 = MORLET_CENTRAL * MORLET_CENTRAL;
      double normalization =
          1.0 / sqrt(sqrt(acos(-1.0)) *
                     (1.0 + exp(-omega2) - 2.0 * exp(-0.75 * omega2)));
      kernel[s][tap] = sqrt(scale) * normalization * integral;
    }
  }
  const int max_width = 99;
  size_t input_length = (size_t)(end - start);
  if (input_length > CWT_WINDOW)
    input_length = CWT_WINDOW;
  input_length += max_width;
  size_t fft_size = 1;
  while (fft_size < input_length + max_width - 1)
    fft_size <<= 1;
  double complex *kernel_fft =
      calloc((size_t)CWT_CHANNELS * fft_size, sizeof(*kernel_fft));
  if (!kernel_fft)
    bye("out of memory");
  for (int s = 0; s < CWT_CHANNELS; s++) {
    for (int tap = 0; tap < width[s]; tap++)
      kernel_fft[(size_t)s * fft_size + width[s] - 1 - tap] = kernel[s][tap];
    fft_radix2(kernel_fft + (size_t)s * fft_size, fft_size, 0);
  }
  const int64_t window_count =
      (end - start + CWT_WINDOW - 1) / CWT_WINDOW;
#pragma omp parallel for schedule(static)
  for (int64_t window = 0; window < window_count; window++) {
    int64_t window_start = start + window * CWT_WINDOW;
    int64_t window_end = window_start + CWT_WINDOW;
    if (window_end > end)
      window_end = end;
    int radius = max_width / 2;
    int64_t input_start = window_start > radius ? window_start - radius : 0;
    int64_t input_end = window_end + radius < length ? window_end + radius
                                                      : length;
    size_t input_count = (size_t)(input_end - input_start);
    double complex *transformed = calloc(fft_size, sizeof(*transformed));
    if (!transformed)
      bye("out of memory");
    for (int s = 0; s < CWT_CHANNELS; s++) {
      for (size_t i = 0; i < input_count; i++) {
        int base = base_index(bases[input_start + (int64_t)i]);
        transformed[i] = base < 0 ? 0.0 : signal[base] / sqrt(2.0);
      }
      memset(transformed + input_count, 0,
             (fft_size - input_count) * sizeof(*transformed));
      fft_radix2(transformed, fft_size, 0);
      for (size_t i = 0; i < fft_size; i++)
        transformed[i] *= kernel_fft[(size_t)s * fft_size + i];
      fft_radix2(transformed, fft_size, 1);
      for (int64_t position = window_start; position < window_end; position++) {
        size_t index = (size_t)(position - input_start) +
                       (size_t)(width[s] - 1 - width[s] / 2);
        double complex coefficient = transformed[index] / (double)fft_size;
        float *value =
            out + (size_t)(position - start) * CWT_CHANNELS * 2 +
                  (size_t)s * 2;
        value[0] = (float)creal(coefficient);
        value[1] = (float)cimag(coefficient);
      }
    }
    free(transformed);
  }
  free(kernel_fft);
}

void cwt_power(const char *bases, int64_t length, int64_t start,
                    int64_t end, float *out) {
  for (int64_t block = start; block < end; block += CWT_WINDOW) {
    float values[CWT_WINDOW * CWT_CHANNELS * COMPLEX_COMPONENTS];
    int64_t stop =
        block + CWT_WINDOW < end ? block + CWT_WINDOW : end;
    cwt(bases, length, block, stop, values);
    for (int64_t i = 0; i < (stop - block) * CWT_CHANNELS; i++)
      out[(size_t)(block - start) * CWT_CHANNELS + (size_t)i] =
          values[COMPLEX_COMPONENTS * i] *
              values[COMPLEX_COMPONENTS * i] +
          values[COMPLEX_COMPONENTS * i + 1] *
              values[COMPLEX_COMPONENTS * i + 1];
  }
}

void compute_features(const Sequence *sequence, const float *power,
                   int64_t power_start, int64_t start, int64_t valid,
                   float *input) {
  float mean[CWT_CHANNELS] = {0}, variance[CWT_CHANNELS] = {0};
  for (int64_t q = 0; q < valid; q++)
    for (int s = 0; s < CWT_CHANNELS; s++)
      mean[s] +=
          power[(size_t)(start + q - power_start) * CWT_CHANNELS + s];
  for (int s = 0; s < CWT_CHANNELS; s++)
    mean[s] /= (float)valid;
  for (int64_t q = 0; q < valid; q++)
    for (int s = 0; s < CWT_CHANNELS; s++) {
      float d =
          power[(size_t)(start + q - power_start) * CWT_CHANNELS + s] -
          mean[s];
      variance[s] += d * d;
    }
  for (int s = 0; s < CWT_CHANNELS; s++)
    variance[s] = sqrtf(variance[s] / (float)valid + 1e-6f);
  for (int64_t q = 0; q < valid; q++) {
    int64_t p = start + q;
    int base = base_index(sequence->bases[p]);
    if (base >= 0)
      input[(size_t)base * (size_t)valid + q] = 1.0f;
    for (int s = 0; s < CWT_CHANNELS; s++)
      input[(size_t)(NUCLEOTIDES + s) * (size_t)valid + q] =
          (power[(size_t)(p - power_start) * CWT_CHANNELS + s] - mean[s]) /
          variance[s];
  }
}
