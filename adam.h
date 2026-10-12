#ifndef ADAM_H
#define ADAM_H

#include "anno.h"
#include <math.h>
#include <stdlib.h>

/* Adam ascent: weights move along the gradient of an objective to maximize. */
typedef struct {
  size_t size;
  long step;
  double *first, *second;
} Adam;

static Adam adam_init(size_t size) {
  return (Adam){size, 0, anno_alloc(size, sizeof(double)), anno_alloc(size, sizeof(double))};
}

static void adam_update(Adam *adam, double *weights, const double *gradient) {
  adam->step++;
  double first_correction = 1 - pow(ADAM_BETA1, (double)adam->step);
  double second_correction = 1 - pow(ADAM_BETA2, (double)adam->step);
  for (size_t i = 0; i < adam->size; i++) {
    adam->first[i] = ADAM_BETA1 * adam->first[i] + (1 - ADAM_BETA1) * gradient[i];
    adam->second[i] = ADAM_BETA2 * adam->second[i] + (1 - ADAM_BETA2) * gradient[i] * gradient[i];
    weights[i] += ADAM_LEARNING_RATE * (adam->first[i] / first_correction) /
                  (sqrt(adam->second[i] / second_correction) + ADAM_EPSILON);
    if (!isfinite(weights[i])) anno_fail("Non-finite training weight");
  }
}

static void adam_free(Adam *adam) {
  free(adam->first), free(adam->second);
}

#endif
