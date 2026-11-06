/**
 * @file lbfgs.h
 * @brief L-BFGS optimizer.
 */
#ifndef LBFGS_H
#define LBFGS_H

#include <stddef.h>

/**
 * @brief Evaluate the objective and its gradient at a parameter vector.
 * @param[in] x: Current parameter vector.
 * @param[out] gradient: Receives the objective gradient.
 * @param[in,out] context: Caller-owned evaluation context.
 * @return Objective value, or a non-finite value on failure.
 */
typedef double (*lbfgs_evaluate)(const double *x, double *gradient,
                                 void *context);

/** Status returned by the L-BFGS optimizer. */
typedef enum {
  LBFGS_MAX_ITERATIONS,     // The configured iteration limit was reached.
  LBFGS_CONVERGED,          // The gradient-based convergence criterion was met.
  LBFGS_STOPPING_CRITERION, // Progress fell below the stopping threshold.
  LBFGS_ALLOCATION_FAILURE = -1, // A required workspace allocation failed.
  LBFGS_NONFINITE =
      -2, // An objective, gradient, or search value was non-finite.
  LBFGS_LINE_SEARCH_FAILURE =
      -3, // The line search could not find an acceptable step.
  LBFGS_INVALID_ARGUMENT = -4 // An input argument was invalid.
} LbfgsStatus;

/**
 * @brief Minimize an objective with limited-memory BFGS.
 * @param[in] n: Number of parameters in the vector.
 * @param[in,out] x: Initial values; receives the best accepted values.
 * @param[in] evaluate: Objective and gradient callback.

 * @param[in,out] context: Caller context passed to evaluate.
 * @param[in] iterations: Maximum accepted iterations.
 * @param[out] completed: Receives the number of accepted iterations.
 * @return Optimizer status, including convergence and failure cases.
 */
LbfgsStatus lbfgs_minimize(size_t n, double *x, lbfgs_evaluate evaluate,
                           void *context, int iterations, int *completed);

#endif /* LBFGS_H */
