#ifndef LBFGS_H
#define LBFGS_H

/* Minimal L-BFGS: two-loop recursion with backtracking line search. */
#include <stddef.h>

typedef double (*lbfgs_evaluate)(const double *x, double *gradient, void *context);
typedef void (*lbfgs_progress)(void *context, int iteration, double objective, size_t n, const double *gradient);
int lbfgs_minimize(size_t n, double *x, lbfgs_evaluate evaluate, lbfgs_progress progress,
                   void *context, int iterations, int history);

#endif

#if defined(LBFGS_IMPLEMENTATION) && !defined(LBFGS_IMPLEMENTED)
#define LBFGS_IMPLEMENTED

#include <math.h>
#include <stdlib.h>
#include <string.h>

static double lbfgs_dot(size_t n, const double *a, const double *b) {
  double sum = 0;
  for (size_t i = 0; i < n; i++) sum += a[i] * b[i];
  return sum;
}

int lbfgs_minimize(size_t n, double *x, lbfgs_evaluate evaluate, lbfgs_progress progress,
                   void *context, int iterations, int history) {
  if (history < 1) history = 1;
  double *gradient = malloc(n * sizeof(*gradient)), *previous_gradient = malloc(n * sizeof(*previous_gradient));
  double *previous_x = malloc(n * sizeof(*previous_x)), *direction = malloc(n * sizeof(*direction));
  double *s = malloc((size_t)history * n * sizeof(*s)), *y = malloc((size_t)history * n * sizeof(*y));
  double *rho = malloc((size_t)history * sizeof(*rho)), *alpha = malloc((size_t)history * sizeof(*alpha));
  if (!gradient || !previous_gradient || !previous_x || !direction || !s || !y || !rho || !alpha) {
    free(gradient), free(previous_gradient), free(previous_x), free(direction);
    free(s), free(y), free(rho), free(alpha);
    return -1;
  }
  double objective = evaluate(x, gradient, context);
  int kept = 0, iteration = 0;
  while (iteration < iterations) {
    /* direction = -H * gradient via the two-loop recursion */
    memcpy(direction, gradient, n * sizeof(*direction));
    for (int i = kept - 1; i >= 0; i--) {
      alpha[i] = rho[i] * lbfgs_dot(n, s + (size_t)i * n, direction);
      for (size_t j = 0; j < n; j++) direction[j] -= alpha[i] * y[(size_t)i * n + j];
    }
    if (kept) {
      const double *sn = s + (size_t)(kept - 1) * n, *yn = y + (size_t)(kept - 1) * n;
      double gamma = lbfgs_dot(n, sn, yn) / lbfgs_dot(n, yn, yn);
      for (size_t j = 0; j < n; j++) direction[j] *= gamma;
    }
    for (int i = 0; i < kept; i++) {
      double beta = rho[i] * lbfgs_dot(n, y + (size_t)i * n, direction);
      for (size_t j = 0; j < n; j++) direction[j] += (alpha[i] - beta) * s[(size_t)i * n + j];
    }
    for (size_t j = 0; j < n; j++) direction[j] = -direction[j];
    double descent = lbfgs_dot(n, gradient, direction);
    if (!(descent < 0)) { /* lost positive definiteness: restart from steepest descent */
      kept = 0;
      for (size_t j = 0; j < n; j++) direction[j] = -gradient[j];
      descent = -lbfgs_dot(n, gradient, gradient);
      if (!descent) break;
    }
    memcpy(previous_x, x, n * sizeof(*x));
    memcpy(previous_gradient, gradient, n * sizeof(*previous_gradient));
    double previous_objective = objective, step = 1;
    int halving = 0;
    while (halving < 40) {
      for (size_t j = 0; j < n; j++) x[j] = previous_x[j] + step * direction[j];
      objective = evaluate(x, gradient, context);
      if (objective <= previous_objective + 1e-4 * step * descent) break;
      step *= 0.5;
      halving++;
    }
    if (halving == 40) {
      memcpy(x, previous_x, n * sizeof(*x));
      memcpy(gradient, previous_gradient, n * sizeof(*gradient));
      break;
    }
    iteration++;
    if (kept == history) { /* drop the oldest pair */
      memmove(s, s + n, (size_t)(history - 1) * n * sizeof(*s));
      memmove(y, y + n, (size_t)(history - 1) * n * sizeof(*y));
      memmove(rho, rho + 1, (size_t)(history - 1) * sizeof(*rho));
      kept--;
    }
    double *sk = s + (size_t)kept * n, *yk = y + (size_t)kept * n;
    for (size_t j = 0; j < n; j++) {
      sk[j] = x[j] - previous_x[j];
      yk[j] = gradient[j] - previous_gradient[j];
    }
    double ys = lbfgs_dot(n, yk, sk);
    if (ys > 1e-10) rho[kept++] = 1 / ys;
    if (progress) progress(context, iteration, objective, n, gradient);
    if (fabs(previous_objective - objective) <= 1e-7 * (1 + fabs(previous_objective))) break;
  }
  free(gradient), free(previous_gradient), free(previous_x), free(direction);
  free(s), free(y), free(rho), free(alpha);
  return iteration;
}

#endif
