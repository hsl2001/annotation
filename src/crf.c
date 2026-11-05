/**
 * @file crf.c
 * @brief Linear-chain CRF objective, emissions, and L-BFGS training.
 */
#include "include/crf_internal.h"

#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/**
 * @brief Compute log(exp(a) + exp(b)) without avoidable overflow.
 * @param[in] a First log-space value.
 * @param[in] b Second log-space value.
 * @return Logarithm of the sum in log space.
 */
static inline double crf__logadd(double a, double b) {
  if (a == -INFINITY) return b;
  if (b == -INFINITY) return a;
  if (a < b) {
    double t = a;
    a = b;
    b = t;
  }
  return a + log1p(exp(b - a));
}

/**
 * @brief Return the flattened emission-parameter offset.
 * @param[in] klass Observed class index.
 * @param[in] channel Input channel or bias index.
 * @return Offset of the class/channel parameter.
 */
static inline size_t crf__emission_index(int klass, int channel) {
  return (size_t)klass * (ANNO_INPUT_CHANNELS + 1) + (size_t)channel;
}

/**
 * @brief Score one class from a dense feature row.
 * @param[in] weights Class-major emission weights and bias.
 * @param[in] features One row of ANNO_INPUT_CHANNELS values.
 * @param[in] klass Observed class index.
 * @return Linear emission score.
 */
static inline double crf__state_emission(const double *weights,
                                         const float *features, int klass) {
  const size_t base = crf__emission_index(klass, 0);
  double score = weights[base + ANNO_INPUT_CHANNELS];
  for (int c = 0; c < ANNO_INPUT_CHANNELS; ++c)
    score += weights[base + (size_t)c] * (double)features[c];
  return score;
}

/**
 * @brief Return the learned score for a class transition.
 * @param[in] weights Class-pair transition weights.
 * @param[in] from_class Previous observed class.
 * @param[in] to_class Next observed class.
 * @return Transition score; same-class transitions have score zero.
 */
static inline double crf__transition(const double *weights, int from_class,
                                     int to_class) {
  if (from_class == to_class) return 0.0;
  return weights[(size_t)from_class * ANNO_CLASSES + (size_t)to_class];
}

/**
 * @brief Evaluate all class emissions for one position.
 * @param[in] model CRF model containing emission weights.
 * @param[in] features One input-feature row.
 * @param[out] scores Receives one score for each observed class.
 * @return None; invalid pointers leave the output unchanged.
 */
void crf_emissions(const CrfModel *model, const float *features,
                   double scores[ANNO_CLASSES]) {
  if (!model || !features || !scores) return;
  for (int c = 0; c < ANNO_CLASSES; ++c)
    scores[c] = crf__state_emission(model->emission_weights, features, c);
}

/**
 * @brief Validate grammar dimensions, classes, and predecessor indices.
 * @param[in] grammar Grammar to inspect.
 * @return Nonzero when the grammar is structurally valid.
 */
static inline int crf__grammar_valid(const CrfGrammar *grammar) {
  if (!grammar || grammar->states <= 0 || grammar->max_predecessors < 0 ||
      !grammar->state_class || !grammar->predecessor_count ||
      (grammar->max_predecessors > 0 && !grammar->predecessors))
    return 0;
  if ((size_t)grammar->states > SIZE_MAX / (size_t)(grammar->max_predecessors
                                                       ? grammar->max_predecessors
                                                       : 1))
    return 0;
  for (int s = 0; s < grammar->states; ++s) {
    if (grammar->state_class[s] >= ANNO_CLASSES ||
        grammar->predecessor_count[s] > grammar->max_predecessors)
      return 0;
    for (int k = 0; k < grammar->predecessor_count[s]; ++k) {
      int p = grammar->predecessors[(size_t)s * grammar->max_predecessors + k];
      if (p < 0 || p >= grammar->states) return 0;
    }
  }
  return 1;
}

/**
 * @brief Score one class directly from a compact training example.
 * @param[in] weights Class-major emission weights and bias.
 * @param[in] example Compact sequence features.
 * @param[in] t Position within the example.
 * @param[in] klass Observed class index.
 * @return Linear emission score.
 */
static inline double crf__example_emission(const double *weights,
                                           const CrfExample *example, size_t t,
                                           int klass) {
  const size_t base = crf__emission_index(klass, 0);
  double score = weights[base + ANNO_INPUT_CHANNELS];
  uint8_t nucleotide = example->bases[t];
  if (nucleotide < ANNO_NUCLEOTIDES)
    score += weights[base + nucleotide];
  const float *cwt = example->cwt + t * ANNO_CWT_CHANNELS;
  for (int c = 0; c < ANNO_CWT_CHANNELS; ++c)
    score += weights[base + ANNO_NUCLEOTIDES + (size_t)c] * cwt[c];
  return score;
}

/** Plan O(states * (block + length / block)) capped partition workspace. */
int crf__workspace_plan(size_t length, size_t states, size_t requested_block,
                        Crf__WorkspacePlan *plan) {
  if (!length || !states || !plan) return 0;
  size_t block = requested_block;
  if (!block) {
    size_t root = (size_t)sqrt((double)length);
    if (!root) root = 1;
    if (root < length / root + (length % root != 0)) ++root;
    block = CRF_CHECKPOINT_BLOCK > 0 ? (size_t)CRF_CHECKPOINT_BLOCK : 1;
    if (block < root) block = root;
  }
  if (block > length) block = length;
  const size_t blocks = (length - 1) / block + 1;
  const size_t max_cells = CRF_MAX_WORK_BYTES / sizeof(double);
  /* (blocks - 1) checkpoints, (block + 1) alpha rows, two beta rows,
   * and block scalar normalizers. Subtraction keeps every size check bounded. */
  if (block > max_cells || block > SIZE_MAX - 2 ||
      blocks > SIZE_MAX - block - 2)
    return 0;
  const size_t rows = blocks + block + 2;
  if (states > (max_cells - block) / rows) return 0;
  plan->block_length = block;
  plan->blocks = blocks;
  plan->bytes = (rows * states + block) * sizeof(double);
  return 1;
}

/**
 * @brief Check example pointers, labels, bases, and workspace bounds.
 * @param[in] example Training example to validate.
 * @param[in] grammar Grammar used to decode the example.
 * @return Nonzero when the example is valid for CRF training.
 */
static inline int crf__example_valid(const CrfExample *example,
                                     const CrfGrammar *grammar) {
  if (!example || example->length < 0 ||
      (example->boundary != CRF_BOUNDARY_FREE &&
       example->boundary != CRF_BOUNDARY_INTERGENIC) ||
      (example->boundary == CRF_BOUNDARY_INTERGENIC &&
       grammar->state_class[0] != ANNO_INTERGENIC))
    return 0;
  if (example->length == 0) return 1;
  if (!example->bases || !example->cwt || !example->labels) return 0;
  if ((uint64_t)example->length > SIZE_MAX / sizeof(float) /
                                      ANNO_CWT_CHANNELS ||
      (uint64_t)example->length > SIZE_MAX / sizeof(uint8_t))
    return 0;
  Crf__WorkspacePlan plan;
  if (!crf__workspace_plan((size_t)example->length, (size_t)grammar->states, 0,
                           &plan))
    return 0;
  for (int64_t t = 0; t < example->length; ++t)
    if (example->labels[t] >= ANNO_CLASSES || example->bases[t] > ANNO_NUCLEOTIDES)
      return 0;
  return 1;
}

/** Compute a max-normalized forward row, retaining its scalar normalizer. */
static int crf__forward_row(const double *parameters, const CrfGrammar *grammar,
                             const CrfExample *example, int gold_only, size_t t,
                             const double *previous, double *row, double *scale) {
  const double *transition = parameters + CRF_EMISSION_WEIGHT_COUNT;
  double emission[ANNO_CLASSES];
  for (int c = 0; c < ANNO_CLASSES; ++c) {
    emission[c] = gold_only && c != example->labels[t]
                      ? 0.0 : crf__example_emission(parameters, example, t, c);
    if (!isfinite(emission[c])) return 0;
  }
  double top = -INFINITY;
  for (int s = 0; s < grammar->states; ++s) {
    const int klass = grammar->state_class[s];
    double value = -INFINITY;
    if ((!gold_only || klass == example->labels[t]) &&
        (t != 0 || example->boundary == CRF_BOUNDARY_FREE || s == 0)) {
      if (t == 0) {
        value = emission[klass];
      } else {
        for (int k = 0; k < grammar->predecessor_count[s]; ++k) {
          const int p = grammar->predecessors[(size_t)s *
                                              grammar->max_predecessors + k];
          value = crf__logadd(value, previous[p] +
              crf__transition(transition, grammar->state_class[p], klass));
        }
        value += emission[klass];
      }
    }
    row[s] = value;
    if (value > top) top = value;
  }
  if (!isfinite(top)) return 0;
  for (int s = 0; s < grammar->states; ++s) row[s] -= top;
  *scale = top;
  return 1;
}

/** Evaluate exact marginals, recomputing alpha blocks from saved checkpoints. */
int crf__partition_gradient(const double *parameters, const CrfGrammar *grammar,
                            const CrfExample *example, int gold_only, double sign,
                            double *gradient, double *log_z_out,
                            size_t checkpoint_block) {
  const size_t states = (size_t)grammar->states;
  const size_t length = (size_t)example->length;
  Crf__WorkspacePlan plan;
  if (!crf__workspace_plan(length, states, checkpoint_block, &plan)) return 0;
  double *work = (double *)malloc(plan.bytes);
  if (!work) return -1;
  const size_t block = plan.block_length;
  double *checkpoints = work;
  double *alpha = checkpoints + (plan.blocks - 1) * states;
  double *beta = alpha + (block + 1) * states;
  double *previous_beta = beta + states;
  double *scales = previous_beta + states;
  const double *transition = parameters + CRF_EMISSION_WEIGHT_COUNT;

  /* The first pass needs only two alpha rows. Checkpoints are the normalized
   * row immediately before each block, not fresh initial distributions. */
  double *previous = alpha, *current = alpha + states;
  double log_scale = 0.0;
  for (size_t t = 0; t < length; ++t) {
    double scale;
    if (!crf__forward_row(parameters, grammar, example, gold_only, t,
                           previous, current, &scale))
      goto nonfinite;
    log_scale += scale;
    if ((t + 1) % block == 0 && t + 1 < length)
      memcpy(checkpoints + ((t + 1) / block - 1) * states, current,
             states * sizeof(*current));
    double *swap = previous;
    previous = current;
    current = swap;
  }
  double terminal = -INFINITY;
  for (size_t s = 0; s < states; ++s)
    if (example->boundary == CRF_BOUNDARY_FREE || s == 0)
      terminal = crf__logadd(terminal, previous[s]);
  if (!isfinite(terminal) || !isfinite(log_scale + terminal)) goto nonfinite;
  *log_z_out = log_scale + terminal;
  for (size_t s = 0; s < states; ++s)
    beta[s] = example->boundary == CRF_BOUNDARY_FREE || s == 0
                  ? -terminal : -INFINITY;

  for (size_t b = plan.blocks; b-- > 0;) {
    const size_t start = b * block;
    const size_t width = length - start < block ? length - start : block;
    if (b)
      memcpy(alpha, checkpoints + (b - 1) * states, states * sizeof(*alpha));
    for (size_t offset = 0; offset < width; ++offset)
      if (!crf__forward_row(parameters, grammar, example, gold_only,
                             start + offset, alpha + offset * states,
                             alpha + (offset + 1) * states, &scales[offset]))
        goto nonfinite;

    for (size_t offset = width; offset-- > 0;) {
      const size_t t = start + offset;
      const double *row = alpha + (offset + 1) * states;
      double emission[ANNO_CLASSES];
      for (int c = 0; c < ANNO_CLASSES; ++c)
        emission[c] = gold_only && c != example->labels[t]
                          ? 0.0 : crf__example_emission(parameters, example, t, c);
      for (size_t p = 0; p < states; ++p) previous_beta[p] = -INFINITY;
      for (size_t s = 0; s < states; ++s) {
        const int klass = grammar->state_class[s];
        if (gold_only && klass != example->labels[t]) continue;
        const double probability = exp(row[s] + beta[s]);
        const size_t base = crf__emission_index(klass, 0);
        const uint8_t nucleotide = example->bases[t];
        if (nucleotide < ANNO_NUCLEOTIDES)
          gradient[base + nucleotide] += sign * probability;
        const float *cwt = example->cwt + t * ANNO_CWT_CHANNELS;
        for (int c = 0; c < ANNO_CWT_CHANNELS; ++c)
          gradient[base + ANNO_NUCLEOTIDES + (size_t)c] +=
              sign * probability * (double)cwt[c];
        gradient[base + ANNO_INPUT_CHANNELS] += sign * probability;

        if (!t) continue;
        const double suffix = emission[klass] + beta[s] - scales[offset];
        for (int k = 0; k < grammar->predecessor_count[s]; ++k) {
          const int p = grammar->predecessors[s *
                              (size_t)grammar->max_predecessors + (size_t)k];
          const int from_class = grammar->state_class[p];
          const double edge = crf__transition(transition, from_class, klass);
          previous_beta[p] = crf__logadd(previous_beta[p], edge + suffix);
          if (from_class != klass) {
            const size_t index = CRF_EMISSION_WEIGHT_COUNT +
                                 (size_t)from_class * ANNO_CLASSES + (size_t)klass;
            gradient[index] += sign * exp(alpha[offset * states + (size_t)p] +
                                           edge + suffix);
          }
        }
      }
      if (t) {
        double *swap = beta;
        beta = previous_beta;
        previous_beta = swap;
      }
    }
  }
  free(work);
  return 1;

nonfinite:
  free(work);
  return 0;
}


/**
 * @brief Compute normalized conditional NLL, regularization, and gradient.
 * @param[in] x Current parameter vector.
 * @param[out] gradient Receives the objective gradient.
 * @param[in,out] opaque CRF objective context; records any worker failure.
 * @return Objective value, or NaN if evaluation fails.
 */
double crf__objective(const double *x, double *gradient, void *opaque) {
  Crf__ObjectiveContext *context = (Crf__ObjectiveContext *)opaque;
  memset(gradient, 0, CRF_PARAMETER_COUNT * sizeof(*gradient));
  context->failure_status = 0;
  const size_t workers = (size_t)context->threads;
  if (!workers || workers > SIZE_MAX / CRF_PARAMETER_COUNT ||
      workers * CRF_PARAMETER_COUNT > SIZE_MAX / sizeof(double)) {
    context->failure_status = LBFGS_INVALID_ARGUMENT;
    return NAN;
  }
  double *thread_gradient =
      (double *)calloc(workers * CRF_PARAMETER_COUNT, sizeof(*thread_gradient));
  double *thread_objective = (double *)calloc(workers, sizeof(*thread_objective));
  int *thread_status = (int *)calloc(workers, sizeof(*thread_status));
  if (!thread_gradient || !thread_objective || !thread_status) {
    free(thread_gradient);
    free(thread_objective);
    free(thread_status);
    context->failure_status = LBFGS_ALLOCATION_FAILURE;
    return NAN;
  }
#ifdef _OPENMP
#pragma omp parallel for num_threads(context->threads) schedule(static)
#endif
  for (size_t i = 0; i < context->example_count; ++i) {
    int thread = 0;
#ifdef _OPENMP
    thread = omp_get_thread_num();
#endif
    double *local_gradient =
        thread_gradient + (size_t)thread * CRF_PARAMETER_COUNT;
    const CrfExample *example = &context->examples[i];
    if (!example->length) continue;
    double log_z_all, log_z_gold;
    int all_status = crf__partition_gradient(
        x, context->grammar, example, 0, 1.0, local_gradient, &log_z_all, 0);
    int gold_status = all_status > 0
        ? crf__partition_gradient(x, context->grammar, example, 1, -1.0,
                                  local_gradient, &log_z_gold, 0)
        : all_status;
    if (all_status <= 0 || gold_status <= 0) {
      int status = all_status < 0 || gold_status < 0
                       ? LBFGS_ALLOCATION_FAILURE : LBFGS_NONFINITE;
      thread_status[thread] = status;
      continue;
    }
    thread_objective[thread] += log_z_all - log_z_gold;
  }
  double objective = 0.0;
  for (size_t thread = 0; thread < workers; ++thread) {
    if (thread_status[thread])
      context->failure_status = thread_status[thread];
    objective += thread_objective[thread];
    for (size_t i = 0; i < CRF_PARAMETER_COUNT; ++i)
      gradient[i] += thread_gradient[thread * CRF_PARAMETER_COUNT + i];
  }
  free(thread_gradient);
  free(thread_objective);
  free(thread_status);
  if (context->failure_status)
    return NAN;
  const double inverse_total_length = 1.0 / (double)context->total_length;
  objective *= inverse_total_length;
  for (size_t i = 0; i < CRF_PARAMETER_COUNT; ++i) {
    gradient[i] *= inverse_total_length;
    objective += 0.5 * CRF_L2_REGULARIZATION * x[i] * x[i];
    gradient[i] += CRF_L2_REGULARIZATION * x[i];
  }
  if (!isfinite(objective)) {
    context->failure_status = LBFGS_NONFINITE;
    return NAN;
  }
  return objective;
}

/**
 * @brief Check that every model parameter is finite.
 * @param[in] model Model to inspect.
 * @return Nonzero when all parameters are finite.
 */
int crf_model_is_finite(const CrfModel *model) {
  for (size_t i = 0; i < CRF_PARAMETER_COUNT; ++i) {
    const double value = i < CRF_EMISSION_WEIGHT_COUNT
                             ? model->emission_weights[i]
                             : model->transition_weights[i - CRF_EMISSION_WEIGHT_COUNT];
    if (!isfinite(value)) return 0;
  }
  return 1;
}

/**
 * @brief Optimize a CRF model with L-BFGS over labeled examples.
 * @param[in,out] model Initial model parameters; replaced with the trained parameters on success.
 * @param[in] grammar Duration-expanded state grammar.
 * @param[in] examples Compact labeled sequences and endpoint policies.
 * @param[in] example_count Number of sequences.
 * @param[in] iterations Maximum L-BFGS iterations.
 * @param[in] threads Number of objective workers.
 * @return An LbfgsStatus value describing convergence or failure.
 */
int crf_train(CrfModel *model, const CrfGrammar *grammar,
              const CrfExample *examples, size_t example_count, int iterations,
              int threads) {
  if (!model || !crf__grammar_valid(grammar) || !examples || !example_count ||
      iterations < 1 || threads < 1 || !crf_model_is_finite(model) ||
      CRF_PARAMETER_COUNT > SIZE_MAX / sizeof(double) ||
      CRF_PARAMETER_COUNT > INT_MAX)
    return LBFGS_INVALID_ARGUMENT;
  uint64_t total_length = 0;
  for (size_t i = 0; i < example_count; ++i) {
    if (!crf__example_valid(&examples[i], grammar)) return LBFGS_INVALID_ARGUMENT;
    const uint64_t length = (uint64_t)examples[i].length;
    if (length > UINT64_MAX - total_length) return LBFGS_INVALID_ARGUMENT;
    total_length += length;
  }
  if (total_length == 0) return LBFGS_INVALID_ARGUMENT;

  double parameters[CRF_PARAMETER_COUNT];
  memcpy(parameters, model->emission_weights, sizeof(model->emission_weights));
  memcpy(parameters + CRF_EMISSION_WEIGHT_COUNT, model->transition_weights,
         sizeof(model->transition_weights));
  if ((size_t)threads > example_count)
    threads = (int)example_count;
  Crf__ObjectiveContext context = {grammar, examples, example_count,
                                   total_length, threads, 0};
  int completed = 0;
  LbfgsStatus status = lbfgs_minimize(CRF_PARAMETER_COUNT, parameters,
                                     crf__objective, NULL, &context, iterations,
                                     &completed);
  if (status == LBFGS_NONFINITE && context.failure_status)
    status = (LbfgsStatus)context.failure_status;
  if (status == LBFGS_MAX_ITERATIONS || status == LBFGS_CONVERGED ||
      status == LBFGS_STOPPING_CRITERION) {
    memcpy(model->emission_weights, parameters, sizeof(model->emission_weights));
    memcpy(model->transition_weights, parameters + CRF_EMISSION_WEIGHT_COUNT,
           sizeof(model->transition_weights));
  }
  return (int)status;
}
