/**
* @file crf.c
* @brief Training linear-chain CRF
*/

#ifndef CRF_C
#define CRF_C

#include "include/surf.h"
#include "include/crf.h"
#include "include/lbfgs.h"
#include "utils.c"

#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/**
 * @brief Prepare samples with window-local feature normalization.
 * @param[in] sequences: Input sequences with per-base labels.
 * @param[in] sequence_count: Number of input sequences.
 * @param[out] sample_count: Receives the number of prepared sequences.
 * @return Allocated samples and their owned feature arrays.
 */
static Sample *prepare_samples(const Sequence *sequences,
                                         int sequence_count,
                                         size_t *sample_count) {
  if (sequence_count < 1)
    bye("no sequences to train CRF");
  size_t count = (size_t)sequence_count;
  Sample *samples = scalloc(count, sizeof(*samples));
  float *input = scalloc((size_t)INPUT_CHANNELS * FEATURE_WINDOW,
                              sizeof(*input));
  float *power = scalloc((size_t)FEATURE_WINDOW * CWT_CHANNELS,
                              sizeof(*power));
  for (int s = 0; s < sequence_count; ++s) {
    const Sequence *sequence = &sequences[s];
    if (!sequence->labels || sequence->length < 1)
      bye("sequence %s has no labels", sequence->name);
    if ((uint64_t)sequence->length > SIZE_MAX)
      bye("sequence %s is too long", sequence->name);
    size_t positions = (size_t)sequence->length;
    size_t cwt_count = smul(positions, CWT_CHANNELS);
    uint8_t *bases = scalloc(positions, sizeof(*bases));
    float *cwt = scalloc(cwt_count, sizeof(*cwt));
    uint8_t *labels = scalloc(positions, sizeof(*labels));
    memcpy(labels, sequence->labels, positions * sizeof(*labels));
    samples[s] = (Sample){sequence->length, bases, cwt, labels};
    for (int64_t start = 0; start < sequence->length;) {
      int64_t length = sequence->length - start;
      if (length > FEATURE_WINDOW)
        length = FEATURE_WINDOW;
      cwt_power(sequence->bases, sequence->length, start, start + length, power);
      memset(input, 0, (size_t)INPUT_CHANNELS * FEATURE_WINDOW *
                          sizeof(*input));
      compute_features(sequence, power, start, start, length, input);
      for (size_t p = 0; p < (size_t)length; ++p) {
        size_t position = (size_t)start + p;
        int base = base_index(sequence->bases[position]);
        bases[position] = (uint8_t)(base < 0 ? NUCLEOTIDES : base);
        for (int c = 0; c < CWT_CHANNELS; ++c)
          cwt[position * CWT_CHANNELS + (size_t)c] =
              input[(size_t)(NUCLEOTIDES + c) * FEATURE_WINDOW + p];
      }
      start += length;
    }
  }
  free(power);
  free(input);
  *sample_count = count;
  return samples;
}

/**
 * @brief Release a set of training samples and all owned feature buffers.
 * @param[in,out] samples: samples to destroy; NULL is accepted.
 * @param[in] count: Number of samples.
 * @return None.
 */
static void prepare_samples_destroy(Sample *samples, size_t count) {
  if (!samples)
    return;
  for (size_t i = 0; i < count; ++i) {
    free((void *)samples[i].bases);
    free((void *)samples[i].cwt);
    free((void *)samples[i].classes);
  }
  free(samples);
}

/**
 * @brief Expose the decoder grammar in the CRF training representation.
 * @param[in] grammar: Expanded decoder grammar (longgrammar).
 * @param[out] classes: Receives the observed class for each grammar state.
 * @param[out] counts: Receives each state's backpointer count.
 * @param[out] backpointers: Receives flattened backpointer indices.
 * @return Grammar view referencing the caller-owned output arrays.
 */
LongGrammar long2short(const ShortGrammar *grammar, uint8_t *classes,
                                   uint8_t *counts, int *backpointers) {
  for (int s = 0; s < grammar->count; ++s) {
    classes[s] = grammar->cls[s];
    counts[s] = (uint8_t)grammar->backpointers[s];
    for (int k = 0; k < grammar->backpointers[s]; ++k)
      backpointers[(size_t)s * SHORTGRAMMAR_MAX_BACKPOINTERS + (size_t)k] =
          grammar->from[s][k];
  }
  return (LongGrammar){grammar->count, classes, SHORTGRAMMAR_MAX_BACKPOINTERS,
                      counts, backpointers};
}

int train(Model *model, const LongGrammar *grammar,
              const Sample *samples, size_t sample_count, int iterations,
              int threads) {
  if (!model || !samples || !sample_count ||
      iterations < 1 || threads < 1 || !is_model_finite(model) ||
      PARAMETERS > SIZE_MAX / sizeof(double))
    return LBFGS_INVALID_ARGUMENT;
  uint64_t total_length = 0;
  for (size_t i = 0; i < sample_count; ++i) {
    const uint64_t length = (uint64_t)samples[i].length;
    if (length > UINT64_MAX - total_length) return LBFGS_INVALID_ARGUMENT;
    total_length += length;
  }
  if (total_length == 0) return LBFGS_INVALID_ARGUMENT;

  double parameters[PARAMETERS];
  memcpy(parameters, model->emission_weights, sizeof(model->emission_weights));
  memcpy(parameters + EMISSIONS, model->transition_weights,
         sizeof(model->transition_weights));
  if ((size_t)threads > sample_count)
    threads = (int)sample_count;
  ModelObjectiveContext context = {grammar, samples, sample_count,
                                   total_length, threads, 0};
  int completed = 0;
  LbfgsStatus status = lbfgs_minimize(PARAMETERS, parameters,
                                     model_objective, NULL, &context, iterations,
                                     &completed);
  if (status == LBFGS_NONFINITE && context.failure_status)
    status = (LbfgsStatus)context.failure_status;
  if (status == LBFGS_MAX_ITERATIONS || status == LBFGS_CONVERGED ||
      status == LBFGS_STOPPING_CRITERION) {
    memcpy(model->emission_weights, parameters, sizeof(model->emission_weights));
    memcpy(model->transition_weights, parameters + EMISSIONS,
           sizeof(model->transition_weights));
  }
  return (int)status;
}

/**
 * @brief Train a CRF from labeled sequences using L-BFGS.
 * @param[in] sequences: Sequences with per-base training labels.
 * @param[in] count: Number of sequences.
 * @param[in] iterations: Maximum optimizer iterations.
 * @param[in] threads: Number of parallel objective workers.
 * @return Newly allocated trained CRF model.
 */
Model *train_pipeline(const Sequence *sequences, int count,
                           int iterations, int threads) {
  ShortGrammar grammar;
  build_grammar(&grammar);
  uint8_t classes[SHORTGRAMMER_MAX_STATES], backpointer_count[SHORTGRAMMER_MAX_STATES];
  int backpointers[SHORTGRAMMER_MAX_STATES * SHORTGRAMMAR_MAX_BACKPOINTERS];
  LongGrammar view = long2short(&grammar, classes, backpointer_count,
                                    backpointers);
  size_t sample_count = 0;
  Sample *sample = prepare_samples(sequences, count, &sample_count);
  Model *model = scalloc(1, sizeof(*model));
  double started = whattime();
  cry("Training: %zu sequences, %d iterations", sample_count,
      iterations);
  int status = train(model, &view, sample, sample_count, iterations,
                         threads);
  prepare_samples_destroy(sample, sample_count);
  if (status != LBFGS_MAX_ITERATIONS && status != LBFGS_CONVERGED &&
      status != LBFGS_STOPPING_CRITERION)
    bye("CRF L-BFGS training failed (status %d)", status);
  cry("CRF L-BFGS status %d after %.1fs", status, whattime() - started);
  return model;
}

/**
 * @brief Compute log(exp(a) + exp(b)) without avoidable overflow.
 * @param[in] a: First log-space value.
 * @param[in] b: Second log-space value.
 * @return Logarithm of the sum in log space.
 */
static inline double logadd(double a, double b) {
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
 * @param[in] klass: Observed class index.
 * @param[in] channel: Input channel or bias index.
 * @return Offset of the class/channel parameter.
 */
static inline size_t emission_index(int klass, int channel) {
  return (size_t)klass * (INPUT_CHANNELS + 1) + (size_t)channel;
}



#endif /** CRF_C */
