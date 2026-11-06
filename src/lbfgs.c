/**
 * @file lbfgs.c
 * @brief L-BFGS optimizer implementation.
 */
#include "include/lbfgs.h"

#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

enum { LBFGS_HISTORY = 20, LBFGS_PAST = 10, LBFGS_MAX_LINE_SEARCH = 64 };
static const double lbfgs_epsilon = 1e-5;
static const double lbfgs_delta = 1e-5;
static const double lbfgs_ftol = 1e-4;
static const double lbfgs_gtol = 0.9;
static const double lbfgs_min_step = 1e-20;
static const double lbfgs_max_step = 1e20;

/**
 * @brief Compute the dot product of two vectors.
 * @param[in] n Number of elements.
 * @param[in] a First vector.
 * @param[in] b Second vector.
 * @return Dot product.
 */
static double lbfgs_dot(size_t n, const double *a, const double *b) {
  double sum = 0;
  for (size_t i = 0; i < n; i++)
    sum += a[i] * b[i];
  return sum;
}

/**
 * @brief Compute a stable Euclidean vector norm.
 * @param[in] n Number of elements.
 * @param[in] values Vector to measure.
 * @return Euclidean norm.
 */
static double lbfgs_norm(size_t n, const double *values) {
  double norm = 0;
  for (size_t i = 0; i < n; i++)
    norm = hypot(norm, values[i]);
  return norm;
}

/**
 * @brief Check whether all vector elements are finite.
 * @param[in] n Number of elements.
 * @param[in] values Vector to inspect.
 * @return Nonzero when every element is finite.
 */
static int lbfgs_finite_vector(size_t n, const double *values) {
  for (size_t i = 0; i < n; i++)
    if (!isfinite(values[i]))
      return 0;
  return 1;
}

/** Expand a descent step, then bisect the bracket using strong-Wolfe zoom. */
static LbfgsStatus lbfgs_line_search(size_t n, double *x, double *gradient,
                                     const double *base_x,
                                     const double *direction, double *objective,
                                     double *step, lbfgs_evaluate evaluate,
                                     void *context) {
  const double initial = *objective;
  const double dginit = lbfgs_dot(n, gradient, direction);
  if (!isfinite(dginit))
    return LBFGS_NONFINITE;
  if (dginit >= 0)
    return LBFGS_LINE_SEARCH_FAILURE;

  double lower = 0, upper = INFINITY, lower_objective = initial;
  *step = fmin(lbfgs_max_step, fmax(lbfgs_min_step, *step));
  for (int trial = 1; trial <= LBFGS_MAX_LINE_SEARCH; ++trial) {
    for (size_t j = 0; j < n; ++j)
      x[j] = base_x[j] + *step * direction[j];
    if (!lbfgs_finite_vector(n, x))
      return LBFGS_NONFINITE;
    double value = evaluate(x, gradient, context);
    if (!isfinite(value) || !lbfgs_finite_vector(n, gradient))
      return LBFGS_NONFINITE;
    double derivative = lbfgs_dot(n, gradient, direction);
    double sufficient = initial + *step * lbfgs_ftol * dginit;
    if (!isfinite(derivative) || !isfinite(sufficient))
      return LBFGS_NONFINITE;

    if (value > sufficient || value >= lower_objective) {
      upper = *step;
    } else {
      if (fabs(derivative) <= -lbfgs_gtol * dginit) {
        *objective = value;

        return LBFGS_CONVERGED;
      }
      if (derivative * (upper - lower) >= 0)
        upper = lower;
      lower = *step;
      lower_objective = value;
    }
    double next = isfinite(upper) ? lower + 0.5 * (upper - lower)
                                  : fmin(lbfgs_max_step, 2 * *step);
    if (next < lbfgs_min_step || next == lower || next == upper)
      return LBFGS_LINE_SEARCH_FAILURE;
    *step = next;
  }
  return LBFGS_LINE_SEARCH_FAILURE;
}

/**
 * @brief Apply the L-BFGS two-loop recursion to form a search direction.
 * @param[in] n Number of parameters.
 * @param[in] gradient Current gradient.
 * @param[in] s Stored parameter-difference vectors.
 * @param[in] y Stored gradient-difference vectors.
 * @param[in] rho Stored reciprocal curvature products.
 * @param[in] kept Number of valid history entries.
 * @param[out] alpha Scratch coefficients for the recursion.
 * @param[out] direction Receives the negative quasi-Newton direction.
 * @return Number of history entries used.
 */
static int lbfgs_direction(size_t n, const double *gradient, const double *s,
                           const double *y, const double *rho, int kept,
                           double *alpha, double *direction) {
  memcpy(direction, gradient, n * sizeof(*direction));
  for (int i = kept - 1; i >= 0; i--) {
    const double *si = s + (size_t)i * n, *yi = y + (size_t)i * n;
    alpha[i] = rho[i] * lbfgs_dot(n, si, direction);
    for (size_t j = 0; j < n; j++)
      direction[j] -= alpha[i] * yi[j];
  }
  if (kept) {
    const double *sn = s + (size_t)(kept - 1) * n;
    const double *yn = y + (size_t)(kept - 1) * n;
    double ys = lbfgs_dot(n, sn, yn), yy = lbfgs_dot(n, yn, yn);
    double gamma = ys > 0 && yy > 0 ? ys / yy : 0;
    if (isfinite(gamma) && gamma > 0) {
      for (size_t j = 0; j < n; j++)
        direction[j] *= gamma;
    } else {
      for (size_t j = 0; j < n; j++)
        direction[j] = gradient[j];
      kept = 0;
    }
  }
  for (int i = 0; i < kept; i++) {
    double beta = rho[i] * lbfgs_dot(n, y + (size_t)i * n, direction);
    for (size_t j = 0; j < n; j++)
      direction[j] += (alpha[i] - beta) * s[(size_t)i * n + j];
  }
  for (size_t j = 0; j < n; j++)
    direction[j] = -direction[j];
  return kept;
}

LbfgsStatus lbfgs_minimize(size_t n, double *x, lbfgs_evaluate evaluate,
                           void *context, int iterations, int *completed) {
  if (completed)
    *completed = 0;
  if (!n || n > INT_MAX || !x || !evaluate || !completed || iterations < 1 ||
      n > SIZE_MAX / sizeof(double) / (4 + 2 * LBFGS_HISTORY))
    return LBFGS_INVALID_ARGUMENT;

  double *work = malloc((4 + 2 * LBFGS_HISTORY) * n * sizeof(*work));
  if (!work)
    return LBFGS_ALLOCATION_FAILURE;
  double *gradient = work, *previous_gradient = gradient + n,
         *previous_x = previous_gradient + n, *direction = previous_x + n,
         *s = direction + n, *y = s + LBFGS_HISTORY * n;
  double rho[LBFGS_HISTORY], alpha[LBFGS_HISTORY];
  double past_objectives[LBFGS_PAST];

  LbfgsStatus status = LBFGS_MAX_ITERATIONS;
  double objective = evaluate(x, gradient, context);
  if (!isfinite(objective) || !lbfgs_finite_vector(n, gradient)) {
    status = LBFGS_NONFINITE;
    goto cleanup;
  }
  past_objectives[0] = objective;

  double xnorm = lbfgs_norm(n, x), gnorm = lbfgs_norm(n, gradient);
  if (!isfinite(xnorm) || !isfinite(gnorm)) {
    status = LBFGS_NONFINITE;
    goto cleanup;
  }
  if (xnorm < 1)
    xnorm = 1;
  if (gnorm / xnorm <= lbfgs_epsilon) {
    status = LBFGS_CONVERGED;
    goto cleanup;
  }

  for (size_t j = 0; j < n; j++)
    direction[j] = -gradient[j];
  double step = 1.0 / gnorm;
  int kept = 0;

  while (*completed < iterations) {
    memcpy(previous_x, x, n * sizeof(*x));
    memcpy(previous_gradient, gradient, n * sizeof(*gradient));

    status = lbfgs_line_search(n, x, gradient, previous_x, direction,
                               &objective, &step, evaluate, context);
    if (status != LBFGS_CONVERGED) {
      memcpy(x, previous_x, n * sizeof(*x));
      goto cleanup;
    }

    (*completed)++;
    xnorm = lbfgs_norm(n, x);
    gnorm = lbfgs_norm(n, gradient);
    if (!isfinite(xnorm) || !isfinite(gnorm)) {
      status = LBFGS_NONFINITE;
      goto cleanup;
    }

    if (xnorm < 1)
      xnorm = 1;
    if (gnorm / xnorm <= lbfgs_epsilon) {
      status = LBFGS_CONVERGED;
      break;
    }

    if (*completed >= LBFGS_PAST && objective > 0) {
      double rate =
          (past_objectives[*completed % LBFGS_PAST] - objective) / objective;
      if (fabs(rate) < lbfgs_delta) {
        status = LBFGS_STOPPING_CRITERION;
        break;
      }
    }
    past_objectives[*completed % LBFGS_PAST] = objective;

    if (*completed >= iterations) {
      status = LBFGS_MAX_ITERATIONS;
      break;
    }

    double ys = 0, yy = 0;
    for (size_t j = 0; j < n; j++) {
      double sj = x[j] - previous_x[j];
      double yj = gradient[j] - previous_gradient[j];
      ys += yj * sj;
      yy += yj * yj;
    }
    if (isfinite(ys) && isfinite(yy) && ys > 1e-10 * yy && isfinite(1 / ys)) {
      if (kept == LBFGS_HISTORY) {
        memmove(s, s + n, (LBFGS_HISTORY - 1) * n * sizeof(*s));
        memmove(y, y + n, (LBFGS_HISTORY - 1) * n * sizeof(*y));
        memmove(rho, rho + 1, (LBFGS_HISTORY - 1) * sizeof(*rho));
        kept--;
      }
      double *sk = s + (size_t)kept * n, *yk = y + (size_t)kept * n;
      for (size_t j = 0; j < n; j++) {
        sk[j] = x[j] - previous_x[j];
        yk[j] = gradient[j] - previous_gradient[j];
      }
      rho[kept++] = 1 / ys;
    }

    kept = lbfgs_direction(n, gradient, s, y, rho, kept, alpha, direction);
    double descent = lbfgs_dot(n, gradient, direction);
    if (!isfinite(descent)) {
      status = LBFGS_NONFINITE;
      goto cleanup;
    }
    if (!(descent < 0)) {
      kept = 0;
      for (size_t j = 0; j < n; j++)
        direction[j] = -gradient[j];
    }
    if (!lbfgs_finite_vector(n, direction)) {
      status = LBFGS_NONFINITE;
      goto cleanup;
    }
    step = kept ? 1.0 : 1.0 / gnorm;
  }

cleanup:
  free(work);
  return status;
}
