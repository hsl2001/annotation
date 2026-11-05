/**
 * @file crf_internal.h
 * @brief Internal CRF objective and model-validation interfaces.
 */
#ifndef CRF_INTERNAL_H
#define CRF_INTERNAL_H

#include "crf.h"

/**
 * @brief Check that every CRF parameter is finite.
 * @param[in] model Model to inspect.
 * @return Nonzero when all parameters are finite.
 */
int crf_model_is_finite(const CrfModel *model);

/** Context shared with the CRF objective callback. */
typedef struct {
  const CrfGrammar *grammar;
  const CrfExample *examples;
  size_t example_count;
  uint64_t total_length;
  int threads;
  int failure_status;
} Crf__ObjectiveContext;

/**
 * @brief Compute normalized conditional NLL, regularization, and gradient.
 * @param[in] x Current parameter vector.
 * @param[out] gradient Receives the objective gradient.
 * @param[in,out] opaque CRF objective context; records any worker failure.
 * @return Objective value, or NaN if evaluation fails.
 */
double crf__objective(const double *x, double *gradient, void *opaque);

#endif /* CRF_INTERNAL_H */
