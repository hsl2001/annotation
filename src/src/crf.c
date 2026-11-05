/**
 * @file crf.c
 * @brief Linear-chain CRF objective, emissions, and L-BFGS training.
 */
#include "crf_internal.h"

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

/**
 * @brief Check example pointers, labels, bases, and workspace bounds.
 * @param[in] example Training example to validate.
 * @param[in] grammar Grammar used to decode the example.
 * @return Nonzero when the example is valid for CRF training.
 */
static inline int crf__example_valid(const CrfExample *example,
                                     const CrfGrammar *grammar) {
  if (!example || example->length < 0) return 0;
  if (example->length == 0) return 1;
  if (!example->bases || !example->cwt || !example->labels) return 0;
  if ((uint64_t)example->length > SIZE_MAX / sizeof(float) /
                                      ANNO_CWT_CHANNELS ||
      (uint64_t)example->length > SIZE_MAX / sizeof(uint8_t))
    return 0;
  size_t cells;
  const size_t max_cells = CRF_MAX_WORK_BYTES / sizeof(double);
  const size_t states = (size_t)grammar->states;
  if ((uint64_t)example->length > SIZE_MAX / states ||
      states > max_cells / 2)
    return 0;
  cells = (size_t)example->length * states;
  if (cells > max_cells - 2 * states) return 0;
  for (int64_t t = 0; t < example->length; ++t)
    if (example->labels[t] >= ANNO_CLASSES || example->bases[t] > ANNO_NUCLEOTIDES)
      return 0;
  return 1;
}

/**
 * @brief Evaluate a grammar partition and add its signed marginals.
 * @param[in] parameters Emission and transition parameter vector.
 * @param[in] grammar Duration-expanded CRF grammar.
 * @param[in] example Sequence features and observed labels.
 * @param[in] gold_only Restrict paths to the observed class at each position.
 * @param[in] sign Sign applied to this partition's gradient contribution.
 * @param[in,out] gradient Accumulates exact state and edge marginals.
 * @param[out] log_z_out Receives the log partition value on success.
 * @return 1 on success, 0 for an invalid/non-finite partition, or -1 on allocation failure.
 */
static inline int crf__partition_gradient(const double *parameters,
                                         const CrfGrammar *grammar,
                                         const CrfExample *example, int gold_only,
                                         double sign, double *gradient,
                                         double *log_z_out) {
  const size_t states = (size_t)grammar->states;
  const size_t length = (size_t)example->length;
  const size_t cells = length * states;
  double *alpha = (double *)malloc(cells * sizeof(*alpha));
  double *beta_a = (double *)malloc(states * sizeof(*beta_a));
  double *beta_b = (double *)malloc(states * sizeof(*beta_b));
  if (!alpha || !beta_a || !beta_b) {
    free(alpha);
    free(beta_a);
    free(beta_b);
    return -1;
  }
  const double *emission_weights = parameters;
  const double *transition_weights = parameters + CRF_EMISSION_WEIGHT_COUNT;
  double log_z = -INFINITY;
  double *beta = beta_a;
  double *previous_beta = beta_b;

  for (size_t t = 0; t < length; ++t) {
    for (size_t s = 0; s < states; ++s) {
      const int klass = grammar->state_class[s];
      double value = -INFINITY;
      if (!gold_only || klass == example->labels[t]) {
        double emission =
            crf__example_emission(emission_weights, example, t, klass);
        if (!isfinite(emission)) goto nonfinite;
        if (t == 0) {
          value = emission;
        } else {
          for (int k = 0; k < grammar->predecessor_count[s]; ++k) {
            const int p = grammar->predecessors[s * (size_t)grammar->max_predecessors +
                                                 (size_t)k];
            const int from_class = grammar->state_class[p];
            double candidate = alpha[(t - 1) * states + (size_t)p] +
                               crf__transition(transition_weights, from_class, klass);
            value = crf__logadd(value, candidate);
          }
          value += emission;
        }
      }
      alpha[t * states + s] = value;
    }
  }

  for (size_t s = 0; s < states; ++s)
    log_z = crf__logadd(log_z, alpha[(length - 1) * states + s]);
  if (!isfinite(log_z)) goto nonfinite;
  *log_z_out = log_z;

  for (size_t s = 0; s < states; ++s) beta_a[s] = 0.0;

  for (size_t t = length; t-- > 0;) {
    for (size_t s = 0; s < states; ++s) {
      const int klass = grammar->state_class[s];
      const double a = alpha[t * states + s];
      if (a != -INFINITY && (!gold_only || klass == example->labels[t])) {
        const double probability = exp(a + beta[s] - log_z);
        const size_t base = crf__emission_index(klass, 0);
        uint8_t nucleotide = example->bases[t];
        if (nucleotide < ANNO_NUCLEOTIDES)
          gradient[base + nucleotide] += sign * probability;
        const float *cwt = example->cwt + t * ANNO_CWT_CHANNELS;
        for (int c = 0; c < ANNO_CWT_CHANNELS; ++c)
          gradient[base + ANNO_NUCLEOTIDES + (size_t)c] +=
              sign * probability * (double)cwt[c];
        gradient[base + ANNO_INPUT_CHANNELS] += sign * probability;
      }
    }

    if (t > 0) {
      for (size_t s = 0; s < states; ++s) {
        const int to_class = grammar->state_class[s];
        if (gold_only && to_class != example->labels[t]) continue;
        const double emission =
            crf__example_emission(emission_weights, example, t, to_class);
        const double suffix = emission + beta[s];
        for (int k = 0; k < grammar->predecessor_count[s]; ++k) {
          const int p = grammar->predecessors[s * (size_t)grammar->max_predecessors +
                                               (size_t)k];
          const double before = alpha[(t - 1) * states + (size_t)p];
          if (before == -INFINITY) continue;
          const int from_class = grammar->state_class[p];
          const double edge = crf__transition(transition_weights, from_class, to_class);
          const double probability = exp(before + edge + suffix - log_z);
          if (from_class != to_class) {
            const size_t index = CRF_EMISSION_WEIGHT_COUNT +
                                 (size_t)from_class * ANNO_CLASSES +
                                 (size_t)to_class;
            gradient[index] += sign * probability;
          }
        }
      }

      for (size_t p = 0; p < states; ++p) previous_beta[p] = -INFINITY;
      for (size_t s = 0; s < states; ++s) {
        const int to_class = grammar->state_class[s];
        if (gold_only && to_class != example->labels[t]) continue;
        const double emission =
            crf__example_emission(emission_weights, example, t, to_class);
        const double suffix = emission + beta[s];
        for (int k = 0; k < grammar->predecessor_count[s]; ++k) {
          const int p = grammar->predecessors[s * (size_t)grammar->max_predecessors +
                                               (size_t)k];
          const int from_class = grammar->state_class[p];
          previous_beta[p] = crf__logadd(
              previous_beta[p],
              crf__transition(transition_weights, from_class, to_class) + suffix);
        }
      }
      double *swap = beta;
      beta = previous_beta;
      previous_beta = swap;
    }
  }

  free(alpha);
  free(beta_a);
  free(beta_b);
  return 1;

nonfinite:
  free(alpha);
  free(beta_a);
  free(beta_b);
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
        x, context->grammar, example, 0, 1.0, local_gradient, &log_z_all);
    int gold_status = all_status > 0
        ? crf__partition_gradient(x, context->grammar, example, 1, -1.0,
                                  local_gradient, &log_z_gold)
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
 * @param[in] examples Compact labeled training windows.
 * @param[in] example_count Number of training windows.
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

