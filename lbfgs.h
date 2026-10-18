#ifndef LBFGS_H
#define LBFGS_H

#include <stddef.h>

typedef double (*lbfgs_evaluate)(const double *x, double *gradient, void *context);
typedef void (*lbfgs_progress)(void *context, int iteration, double objective, double step,
                               int line_search_trials, size_t n, const double *gradient);

typedef enum {
  LBFGS_MAX_ITERATIONS,
  LBFGS_CONVERGED,
  LBFGS_STOPPING_CRITERION,
  LBFGS_ALLOCATION_FAILURE = -1,
  LBFGS_NONFINITE = -2,
  LBFGS_LINE_SEARCH_FAILURE = -3,
  LBFGS_INVALID_ARGUMENT = -4
} LbfgsStatus;

LbfgsStatus lbfgs_minimize(size_t n, double *x, lbfgs_evaluate evaluate,
                           lbfgs_progress progress, void *context, int iterations,
                           int *completed);

#endif

#if defined(LBFGS_IMPLEMENTATION) && !defined(LBFGS_IMPLEMENTED)
#define LBFGS_IMPLEMENTED

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

enum { LBFGS_HISTORY = 6, LBFGS_PAST = 10, LBFGS_MAX_LINE_SEARCH = 20 };
static const double lbfgs_epsilon = 1e-5;
static const double lbfgs_delta = 1e-5;
static const double lbfgs_ftol = 1e-4;
static const double lbfgs_gtol = 0.9;
static const double lbfgs_xtol = 1e-16;
static const double lbfgs_min_step = 1e-20;
static const double lbfgs_max_step = 1e20;

static double lbfgs_dot(size_t n, const double *a, const double *b) {
  double sum = 0;
  for (size_t i = 0; i < n; i++) sum += a[i] * b[i];
  return sum;
}

static double lbfgs_norm(size_t n, const double *values) {
  double norm = 0;
  for (size_t i = 0; i < n; i++) norm = hypot(norm, values[i]);
  return norm;
}

static int lbfgs_finite_vector(size_t n, const double *values) {
  for (size_t i = 0; i < n; i++)
    if (!isfinite(values[i])) return 0;
  return 1;
}

static double lbfgs_cubic_step(double x, double fx, double dx,
                               double y, double fy, double dy) {
  double distance = y - x;
  if (distance == 0) return NAN;
  double theta = 3 * (fx - fy) / distance + dx + dy;
  double scale = fmax(fabs(theta), fmax(fabs(dx), fabs(dy)));
  if (!(scale > 0) || !isfinite(scale)) return NAN;
  double a = theta / scale, b = dx / scale, c = dy / scale;
  double gamma = scale * sqrt(fmax(0, a * a - b * c));
  if (y < x) gamma = -gamma;
  double numerator = gamma - dx + theta;
  double denominator = gamma - dx + gamma + dy;
  if (denominator == 0) return NAN;
  return x + (numerator / denominator) * distance;
}

static double lbfgs_quadratic_step(double x, double fx, double dx,
                                   double t, double ft) {
  double distance = t - x;
  if (distance == 0) return NAN;
  double denominator = (fx - ft) / distance + dx;
  if (denominator == 0) return NAN;
  return x + (0.5 * dx / denominator) * distance;
}

static double lbfgs_secant_step(double x, double dx, double t, double dt) {
  double denominator = dt - dx;
  if (denominator == 0) return NAN;
  return t + (dt / denominator) * (x - t);
}

static int lbfgs_update_interval(double *x, double *fx, double *dx,
                                  double *y, double *fy, double *dy,
                                  double *step, double fp, double dp,
                                  double lower, double upper, int *bracketed) {
  if (*bracketed &&
      (*step <= fmin(*x, *y) || *step >= fmax(*x, *y) ||
       *dx * (*step - *x) >= 0 || upper <= lower))
    return -1;

  int bound;
  double sgnd = dp * (*dx >= 0 ? 1 : -1);
  double cubic, secant, next;
  if (fp > *fx) {
    cubic = lbfgs_cubic_step(*x, *fx, *dx, *step, fp, dp);
    secant = lbfgs_quadratic_step(*x, *fx, *dx, *step, fp);
    if (!isfinite(cubic)) cubic = 0.5 * (*x + *step);
    if (!isfinite(secant)) secant = 0.5 * (*x + *step);
    next = fabs(cubic - *x) < fabs(secant - *x)
      ? cubic : cubic + 0.5 * (secant - cubic);
    *bracketed = 1;
    bound = 1;
  } else if (sgnd < 0) {
    cubic = lbfgs_cubic_step(*step, fp, dp, *x, *fx, *dx);
    secant = lbfgs_secant_step(*step, dp, *x, *dx);
    if (!isfinite(cubic)) cubic = 0.5 * (*step + *x);
    if (!isfinite(secant)) secant = 0.5 * (*step + *x);
    next = fabs(cubic - *step) > fabs(secant - *step) ? cubic : secant;
    *bracketed = 1;
    bound = 0;
  } else if (fabs(dp) < fabs(*dx)) {
    double step_limit = *step > *x
      ? fmin(upper, *step + 4 * (*step - *x))
      : fmax(lower, *step + 4 * (*step - *x));
    double distance = *step - *x;
    double theta = distance == 0 ? NAN : 3 * (*fx - fp) / distance + *dx + dp;
    double scale = fmax(fabs(theta), fmax(fabs(*dx), fabs(dp)));
    double gamma = NAN, ratio = NAN;
    if (scale > 0 && isfinite(scale)) {
      double a = theta / scale, b = *dx / scale, c = dp / scale;
      gamma = scale * sqrt(fmax(0, a * a - b * c));
      if (*step > *x) gamma = -gamma;
      double denominator = gamma + *dx - dp + gamma;
      if (denominator != 0) ratio = (gamma - dp + theta) / denominator;
    }
    cubic = ratio < 0 && gamma != 0
      ? *step + ratio * (*x - *step)
      : step_limit;
    secant = lbfgs_secant_step(*step, dp, *x, *dx);
    if (!isfinite(cubic)) cubic = step_limit;
    if (!isfinite(secant)) secant = step_limit;
    if (*bracketed) {
      next = fabs(*step - cubic) < fabs(*step - secant) ? cubic : secant;
    } else {
      next = fabs(*step - cubic) > fabs(*step - secant) ? cubic : secant;
    }
    next = *step > *x ? fmin(step_limit, next) : fmax(step_limit, next);
    bound = 1;
  } else {
    if (*bracketed) {
      next = lbfgs_cubic_step(*step, fp, dp, *y, *fy, *dy);
      if (!isfinite(next)) next = 0.5 * (*step + *y);
    } else {
      next = *step > *x ? upper : lower;
    }
    bound = 0;
  }

  if (fp > *fx) {
    *y = *step;
    *fy = fp;
    *dy = dp;
  } else {
    if (sgnd < 0) {
      *y = *x;
      *fy = *fx;
      *dy = *dx;
    }
    *x = *step;
    *fx = fp;
    *dx = dp;
  }

  if (*bracketed) {
    double new_lower = fmin(*x, *y);
    double new_upper = fmax(*x, *y);
    next = fmin(new_upper, fmax(new_lower, next));
    if (bound) {
      double interior = *x + 0.66 * (*y - *x);
      next = *y > *x ? fmin(interior, next) : fmax(interior, next);
    }
    if (!(next > new_lower && next < new_upper))
      next = *x + 0.5 * (*y - *x);
  } else {
    next = fmin(upper, fmax(lower, next));
  }
  if (!isfinite(next)) return -1;
  *step = next;
  return 0;
}

static LbfgsStatus lbfgs_morethuente(size_t n, double *x, double *gradient,
                                    const double *base_x, const double *direction,
                                    double *objective, double *step,
                                    lbfgs_evaluate evaluate, void *context,
                                    int *trials) {
  double dginit = lbfgs_dot(n, gradient, direction);
  if (!isfinite(dginit)) return LBFGS_NONFINITE;
  if (dginit >= 0) return LBFGS_LINE_SEARCH_FAILURE;

  int bracketed = 0, stage1 = 1, update_error = 0;
  double initial_objective = *objective;
  double dgtest = lbfgs_ftol * dginit;
  double width = lbfgs_max_step - lbfgs_min_step;
  double previous_width = 2 * width;
  double stx = 0, fx = initial_objective, dgx = dginit;
  double sty = 0, fy = initial_objective, dgy = dginit;
  double best_step = 0, best_objective = initial_objective, last_step = 0;
  double trial_objective = initial_objective;
  int count = 0;

  while (count < LBFGS_MAX_LINE_SEARCH) {
    double lower = bracketed ? fmin(stx, sty) : stx;
    double upper = bracketed ? fmax(stx, sty) : *step + 4 * (*step - stx);
    *step = fmin(lbfgs_max_step, fmax(lbfgs_min_step, *step));
    if (bracketed &&
        (*step <= lower || *step >= upper || update_error ||
         upper - lower <= lbfgs_xtol * fmax(1.0, upper)))
      break;

    for (size_t j = 0; j < n; j++) x[j] = base_x[j] + *step * direction[j];
    trial_objective = evaluate(x, gradient, context);
    if (!isfinite(trial_objective) || !lbfgs_finite_vector(n, gradient))
      return LBFGS_NONFINITE;
    double dg = lbfgs_dot(n, gradient, direction);
    if (!isfinite(dg)) return LBFGS_NONFINITE;
    double ftest = initial_objective + *step * dgtest;
    if (!isfinite(ftest)) return LBFGS_NONFINITE;
    last_step = *step;
    count++;

    if (trial_objective <= ftest && trial_objective < best_objective) {
      best_step = *step;
      best_objective = trial_objective;
    }

    if (trial_objective <= ftest && fabs(dg) <= lbfgs_gtol * (-dginit)) {
      *objective = trial_objective;
      *trials = count;
      return LBFGS_CONVERGED;
    }
    if (*step == lbfgs_max_step && trial_objective <= ftest && dg <= dgtest)
      break;
    if (*step == lbfgs_min_step && (ftest < trial_objective || dgtest <= dg))
      break;
    if (bracketed && upper - lower <= lbfgs_xtol * fmax(1.0, upper))
      break;
    if (count >= LBFGS_MAX_LINE_SEARCH)
      break;

    if (stage1 && trial_objective <= ftest &&
        fmin(lbfgs_ftol, lbfgs_gtol) * dginit <= dg)
      stage1 = 0;

    if (stage1 && ftest < trial_objective && trial_objective <= fx) {
      double fm = trial_objective - *step * dgtest;
      double fxm = fx - stx * dgtest;
      double fym = fy - sty * dgtest;
      double dgm = dg - dgtest;
      double dgxm = dgx - dgtest;
      double dgym = dgy - dgtest;
      update_error = lbfgs_update_interval(&stx, &fxm, &dgxm,
                                           &sty, &fym, &dgym, step, fm, dgm,
                                           lower, upper, &bracketed);
      fx = fxm + stx * dgtest;
      fy = fym + sty * dgtest;
      dgx = dgxm + dgtest;
      dgy = dgym + dgtest;
    } else {
      update_error = lbfgs_update_interval(&stx, &fx, &dgx, &sty, &fy, &dgy,
                                           step, trial_objective, dg,
                                           lower, upper, &bracketed);
    }

    if (bracketed) {
      if (0.66 * previous_width <= fabs(sty - stx))
        *step = stx + 0.5 * (sty - stx);
      previous_width = width;
      width = fabs(sty - stx);
    }
  }

  if (best_step > 0) {
    *step = best_step;
    if (last_step != best_step) {
      for (size_t j = 0; j < n; j++) x[j] = base_x[j] + *step * direction[j];
      trial_objective = evaluate(x, gradient, context);
      if (!isfinite(trial_objective) || !lbfgs_finite_vector(n, gradient))
        return LBFGS_NONFINITE;
      count++;
    }
    if (trial_objective <= initial_objective + *step * dgtest) {
      *objective = trial_objective;
      *trials = count;
      return LBFGS_CONVERGED;
    }
  }
  return LBFGS_LINE_SEARCH_FAILURE;
}

static int lbfgs_direction(size_t n, const double *gradient,
                           const double *s, const double *y, const double *rho,
                           int kept, double *alpha, double *direction) {
  memcpy(direction, gradient, n * sizeof(*direction));
  for (int i = kept - 1; i >= 0; i--) {
    const double *si = s + (size_t)i * n, *yi = y + (size_t)i * n;
    alpha[i] = rho[i] * lbfgs_dot(n, si, direction);
    for (size_t j = 0; j < n; j++) direction[j] -= alpha[i] * yi[j];
  }
  if (kept) {
    const double *sn = s + (size_t)(kept - 1) * n;
    const double *yn = y + (size_t)(kept - 1) * n;
    double ys = lbfgs_dot(n, sn, yn), yy = lbfgs_dot(n, yn, yn);
    double gamma = ys > 0 && yy > 0 ? ys / yy : 0;
    if (isfinite(gamma) && gamma > 0) {
      for (size_t j = 0; j < n; j++) direction[j] *= gamma;
    } else {
      for (size_t j = 0; j < n; j++) direction[j] = gradient[j];
      kept = 0;
    }
  }
  for (int i = 0; i < kept; i++) {
    double beta = rho[i] * lbfgs_dot(n, y + (size_t)i * n, direction);
    for (size_t j = 0; j < n; j++)
      direction[j] += (alpha[i] - beta) * s[(size_t)i * n + j];
  }
  for (size_t j = 0; j < n; j++) direction[j] = -direction[j];
  return kept;
}

LbfgsStatus lbfgs_minimize(size_t n, double *x, lbfgs_evaluate evaluate,
                           lbfgs_progress progress, void *context, int iterations,
                           int *completed) {
  if (completed) *completed = 0;
  if (!n || n > INT_MAX || !x || !evaluate || !completed || iterations < 1 ||
      n > SIZE_MAX / sizeof(double) ||
      (size_t)LBFGS_HISTORY > SIZE_MAX / n / sizeof(double))
    return LBFGS_INVALID_ARGUMENT;

  size_t history_size = (size_t)LBFGS_HISTORY * n;
  double *gradient = malloc(n * sizeof(*gradient));
  double *previous_gradient = malloc(n * sizeof(*previous_gradient));
  double *previous_x = malloc(n * sizeof(*previous_x));
  double *direction = malloc(n * sizeof(*direction));
  double *s = malloc(history_size * sizeof(*s));
  double *y = malloc(history_size * sizeof(*y));
  double *rho = malloc(LBFGS_HISTORY * sizeof(*rho));
  double *alpha = malloc(LBFGS_HISTORY * sizeof(*alpha));
  double *past_objectives = malloc(LBFGS_PAST * sizeof(*past_objectives));
  if (!gradient || !previous_gradient || !previous_x || !direction || !s || !y ||
      !rho || !alpha || !past_objectives) {
    free(gradient), free(previous_gradient), free(previous_x), free(direction);
    free(s), free(y), free(rho), free(alpha), free(past_objectives);
    return LBFGS_ALLOCATION_FAILURE;
  }

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
  if (xnorm < 1) xnorm = 1;
  if (gnorm / xnorm <= lbfgs_epsilon) {
    status = LBFGS_CONVERGED;
    goto cleanup;
  }

  for (size_t j = 0; j < n; j++) direction[j] = -gradient[j];
  double step = 1.0 / gnorm;
  int kept = 0;

  while (*completed < iterations) {
    memcpy(previous_x, x, n * sizeof(*x));
    memcpy(previous_gradient, gradient, n * sizeof(*gradient));
    double previous_objective = objective;
    int trials = 0;
    status = lbfgs_morethuente(n, x, gradient, previous_x, direction,
                               &objective, &step, evaluate, context, &trials);
    if (status == LBFGS_LINE_SEARCH_FAILURE && kept > 0) {
      kept = 0;
      memcpy(x, previous_x, n * sizeof(*x));
      memcpy(gradient, previous_gradient, n * sizeof(*gradient));
      objective = previous_objective;
      for (size_t j = 0; j < n; j++) direction[j] = -gradient[j];
      step = 1.0 / gnorm;
      status = lbfgs_morethuente(n, x, gradient, previous_x, direction,
                                 &objective, &step, evaluate, context, &trials);
    }
    if (status != LBFGS_CONVERGED) {
      memcpy(x, previous_x, n * sizeof(*x));
      memcpy(gradient, previous_gradient, n * sizeof(*gradient));
      objective = previous_objective;
      goto cleanup;
    }

    (*completed)++;
    xnorm = lbfgs_norm(n, x);
    gnorm = lbfgs_norm(n, gradient);
    if (!isfinite(xnorm) || !isfinite(gnorm)) {
      status = LBFGS_NONFINITE;
      goto cleanup;
    }
    if (progress) progress(context, *completed, objective, step, trials, n, gradient);

    if (xnorm < 1) xnorm = 1;
    if (gnorm / xnorm <= lbfgs_epsilon) {
      status = LBFGS_CONVERGED;
      break;
    }

    if (*completed >= LBFGS_PAST && objective > 0) {
      double rate = (past_objectives[*completed % LBFGS_PAST] - objective) / objective;
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
      for (size_t j = 0; j < n; j++) direction[j] = -gradient[j];
    }
    if (!lbfgs_finite_vector(n, direction)) {
      status = LBFGS_NONFINITE;
      goto cleanup;
    }
    step = kept ? 1.0 : 1.0 / gnorm;
  }

cleanup:
  free(gradient), free(previous_gradient), free(previous_x), free(direction);
  free(s), free(y), free(rho), free(alpha), free(past_objectives);
  return status;
}

#endif
