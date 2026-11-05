/**
 * @file crf.h
 * @brief Linear-chain CRF
 */
#ifndef CRF_H
#define CRF_H

#include "surf.h"

#include <stddef.h>
#include <stdint.h>

/** Parameters of CRF */
#define EMISSIONS (CLASSES * (INPUT_CHANNELS + 1))
#define TRANSITIONS (CLASSES * CLASSES)
#define PARAMETERS (EMMISIONS * TRANSITIONS)
#define L2 1e-4

/** Learned emission and transition parameters for the observed classes. */
typedef struct {
  double emission_weights[EMISSIONS]; //Class/channel weights and biases.
  double transition_weights[EMISSIONS]; // Class-pair transition weights.
} Model;

/** Just one sampled, compact, and position-major feature window */
typedef struct {
    int64_t length; // Number of labeled positions
    const uint8_t *base; // A, C, G, T as 0..3, (4 for ambiguous); CAUTION: other implemations are 0, 1, 2, 3 and -1
    const float *cwt; // Position-major normalized CWT power array
    const uint8_t *classes; // Observed class labels per position
} One;

/** Expanede grammar for decoding (includes inner exons) */
typedef struct {
    int states;// Number of expanded states
    const uint8_t *state_class; // Observed class emitted by each state
    int max_backpointers; // Row stride for backpointer array
    const uint8_t *backpointer_count; // Incoming-edge per state
    const int *backpointers; // Flattened state-major backpointer indcies
} LongGrammar;

/**
 * @brief Evaluate all class emissions for one position.
 * @param[in] model: CRF model containing emission weights.
 * @param[in] features: One input-feature row.
 * @param[out] scores: Receives one score for each observed class.
 * @return None; invalid pointers leave the output unchanged.
 */
void evaluate_emissions(const Model *model, const float *features,
                   double scores[CLASSES]);

/**
 * @brief Optimize a CRF model with L-BFGS over labeled one examples.
 * @param[in,out] model: Initial model parameters; replaced with the trained parameters on success.
 * @param[in] grammar: Duration-expanded state grammar.
 * @param[in] ones: Sampled labeled training windows.
 * @param[in] one_count Number of training windows.
 * @param[in] iterations: Maximum L-BFGS iterations.
 * @param[in] threads: Number of objective workers.
 * @return An LbfgsStatus value describing convergence or failure.
 */
int train_model(Model *model, const LongGrammar *grammar,
              const One *ones, size_t one_count, int iterations,
              int threads);

/**
 * @brief Check that every CRF parameter is finite.
 * @param[in] model: Model to inspect.
 * @return Nonzero when all parameters are finite.
 */
int model_is_finite(const Model *model);

/** Context shared with the CRF objective callback. */
typedef struct {
  const LongGrammar *grammar;
  const One *examples;
  size_t example_count;
  uint64_t total_length;
  int threads;
  int failure_status;
} ModelObjectiveContext;

/**
 * @brief Compute normalized conditional NLL, regularization, and gradient.
 * @param[in] x: Current parameter vector.
 * @param[out] gradient: Receives the objective gradient.
 * @param[in,out] opaque: CRF objective context; records any worker failure.
 * @return Objective value, or NaN if evaluation fails.
 */
double model_objective(const double *x, double *gradient, void *opaque);

#endif /* CRF_H */
