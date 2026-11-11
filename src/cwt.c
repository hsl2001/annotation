/**
 * @file cwt.c
 * @brief CWT transforms and per-base model features.
 */

#include "include/cwt.h"
#include "include/surf.h"
#include "include/utils.h"

#include <complex.h>
#include <math.h>

#include <string.h>

/** Morlet scales of the CWT channels. */
enum { CWT_MAX_SCALE = 5 };
static const int cwt_scales[CWT_CHANNELS] = {3, 4, CWT_MAX_SCALE};

void cwt_power(const char *bases, int64_t length, int64_t start, int64_t end,
               float *out) {
  static const double complex signal[NUCLEOTIDES] = {1.0 + I, -1.0 - I,
                                                     -1.0 + I, 1.0 - I};
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

  const double support = 2.0 * CWT_RADIUS;
  const double omega2 = MORLET_CENTRAL * MORLET_CENTRAL;
  const double correction = exp(-0.5 * omega2);
  const double normalization =
      1.0 /
      sqrt(sqrt(acos(-1.0)) * (1.0 + exp(-omega2) - 2.0 * exp(-0.75 * omega2)));
  double complex kernel[CWT_CHANNELS][4 * CWT_RADIUS * CWT_MAX_SCALE + 3];
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
#pragma omp parallel for schedule(static)
  for (int64_t position = start; position < end; position++) {
    for (int s = 0; s < CWT_CHANNELS; s++) {
      double complex coefficient = 0.0;
      for (int tap = 0; tap < width[s]; tap++) {
        int64_t p = position + tap - width[s] / 2;
        if (p < 0 || p >= length)
          continue;
        int base = base_index(bases[p]);
        if (base >= 0)
          coefficient += signal[base] / sqrt(2.0) * kernel[s][tap];
      }
      float real = (float)creal(coefficient), imag = (float)cimag(coefficient);
      size_t index = (size_t)(position - start) * CWT_CHANNELS + (size_t)s;
      out[index] = real * real + imag * imag;
    }
  }
}

void compute_features(char base, const float *power, float *features) {
  memset(features, 0, NUCLEOTIDES * sizeof(*features));
  int index = base_index(base);
  if (index >= 0)
    features[index] = 1.0f;
  memcpy(features + NUCLEOTIDES, power, CWT_CHANNELS * sizeof(*power));
}
