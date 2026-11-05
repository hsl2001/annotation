/**
 * @file crf.h
 * @brief Linear-chain CRF types, configuration, and application-facing API.
 */
#ifndef CRF_H
#define CRF_H

#include "anno.h"
#include "lbfgs.h"

#include <stddef.h>
#include <stdint.h>

#define CRF_EMISSION_WEIGHT_COUNT (ANNO_CLASSES * (ANNO_INPUT_CHANNELS + 1))
#define CRF_TRANSITION_WEIGHT_COUNT (ANNO_CLASSES * ANNO_CLASSES)
#define CRF_PARAMETER_COUNT (CRF_EMISSION_WEIGHT_COUNT + CRF_TRANSITION_WEIGHT_COUNT)
#ifndef CRF_L2_REGULARIZATION
#define CRF_L2_REGULARIZATION 1e-4
#endif
#ifndef CRF_MAX_WORK_BYTES
#define CRF_MAX_WORK_BYTES ((size_t)512 * 1024 * 1024)
#endif

/** Learned emission and transition parameters for the seven observed classes. */
typedef struct {
  double emission_weights[CRF_EMISSION_WEIGHT_COUNT]; /**< Class/channel weights and biases. */
  double transition_weights[CRF_TRANSITION_WEIGHT_COUNT]; /**< Class-pair transition weights. */
} CrfModel;

/** One compact, position-major labeled training window. */
typedef struct {
  int64_t length; /**< Number of labeled positions. */
  const uint8_t *bases; /**< A/C/G/T as 0..3, ambiguous as 4. */
  const float *cwt; /**< Position-major normalized CWT power. */
  const uint8_t *labels; /**< Per-position observed class labels. */
} CrfExample;

/** Duration-expanded state grammar and its incoming-edge lists. */
typedef struct {
  int states; /**< Number of duration-expanded states. */
  const uint8_t *state_class; /**< Observed class emitted by each state. */
  int max_predecessors; /**< Row stride of the predecessor array. */
  const uint8_t *predecessor_count; /**< Incoming-edge count per state. */
  const int *predecessors; /**< Flattened state-major predecessor indices. */
} CrfGrammar;

/**
 * @brief Evaluate all class emissions for one position.
 * @param[in] model CRF model containing emission weights.
 * @param[in] features One input-feature row.
 * @param[out] scores Receives one score for each observed class.
 * @return None; invalid pointers leave the output unchanged.
 */
void crf_emissions(const CrfModel *model, const float *features,
                   double scores[ANNO_CLASSES]);

/**
 * @brief Optimize a CRF model with L-BFGS over labeled examples.
 * @param[in,out] model Initial model parameters; replaced with the trained parameters on success.
 * @param[in] grammar Duration-expanded state grammar.
 * @param[in] examples Compact labeled training windows.
 * @param[in] example_count Number of training windows.
 * @param[in] iterations Maximum L-BFGS iterations.
 * @param[in] threads Number of objective workers.
 * @return An LbfgsStatus value describing convergence or failure.
 */
int crf_train(CrfModel *model, const CrfGrammar *grammar,
              const CrfExample *examples, size_t example_count, int iterations,
              int threads);



#endif /* CRF_H */
