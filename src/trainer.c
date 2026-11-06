/**
* @file crf.c
* @brief Training linear-chain CRF
*/

#ifndef CRF_C
#define CRF_C

#include "include/surf.h"
#include "include/crf.h"
#include "include/lbfgs.h"
#include "include/utils.h"

#include <math.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>

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
 * @param[in] class: Observed class index.
 * @param[in] channel: Input channel or bias index.
 * @return Offset of the class/channel parameter.
 */
static inline size_t emission_index(int class, int channel) {
  return (size_t)class * (INPUT_CHANNELS + 1) + (size_t)channel;
}

/**
 * @brief Score one class from a dense feature row.
 * @param[in] weights: Class-major emission weights and bias.
 * @param[in] features: One row of INPUT_CHANNELS values.
 * @param[in] class: Observed class index.
 * @return Linear emission score.
 */
static inline double emission_scores(const double *weights,
                                         const float *features, int class) {
  const size_t base = emission_index(class, 0);
  double score = weights[base + INPUT_CHANNELS];
  for (int c = 0; c < INPUT_CHANNELS; ++c)
    score += weights[base + (size_t)c] * (double)features[c];
  return score;
}

/**
 * @brief Score classes of samples directly from a training sample.
 * @param[in] weights: Class-major emission weights and bias.
 * @param[in] sample: Compact sequence features.
 * @param[in] t: Position within the sample.
 * @param[in] class: Observed class index.
 * @return Linear emission score.
 */
static inline double sample_emission_scores(const double *weights,
                                           const Sample *sample, size_t t,
                                           int class) {
  const size_t base = emission_index(class, 0);
  double score = weights[base + INPUT_CHANNELS];
  uint8_t nucleotide = sample->bases[t];
  if (nucleotide < NUCLEOTIDES)
    score += weights[base + nucleotide];
  const float *cwt = sample->cwt + t * CWT_CHANNELS;
  for (int c = 0; c < CWT_CHANNELS; ++c)
    score += weights[base + NUCLEOTIDES + (size_t)c] * cwt[c];
  return score;
}

/**
 * @brief Evaluate all class emissions for one position.
 * @param[in] model CRF model containing emission weights.
 * @param[in] features One input-feature row.
 * @param[out] scores Receives one score for each observed class.
 * @return None; invalid pointers leave the output unchanged.
 */
void compute_emissions(const Model *model, const float *features,
                       double scores[CLASSES]) {
  if (!model || !features || !scores) return;
  for (int c = 0; c < CLASSES; ++c)
    scores[c] = emission_scores(model->emission_weights, features, c);
}

/**
 * @brief Return the learned score for a class transition.
 * @param[in] weights: Class-pair transition weights.
 * @param[in] from_class: Previous observed class.
 * @param[in] to_class: Next observed class.
 * @return Transition score; same-class transitions have score zero.
 */
static inline double transition_score(const double *weights, int from_class,
                                     int to_class) {
  if (from_class == to_class) return 0.0;
  return weights[(size_t)from_class * CLASSES + (size_t)to_class];
}

static inline int model_class(int class) {
  return class < CLASSES ? class : INTERGENIC;
}

static int forward_row(const double *parameters, const LongGrammar *grammar,
                             const Sample *sample, int gold_only, size_t t,
                             const double *previous, double *row, double *scale) {
  const double *transition = parameters + EMISSIONS;
  double emission[CLASSES];
  for (int c = 0; c < CLASSES; ++c) {
    emission[c] = gold_only && c != sample->classes[t]
                      ? 0.0 : sample_emission_scores(parameters, sample, t, c);
    if (!isfinite(emission[c])) return 0;
  }
  double top = -INFINITY;
  for (int s = 0; s < grammar->states; ++s) {
    const int state_class = grammar->state_class[s];
    const int class = model_class(state_class);
    double value = -INFINITY;
    if ((!gold_only || state_class == sample->classes[t])) {
      if (t == 0) {
        value = emission[class];
      } else {
        for (int k = 0; k < grammar->backpointer_count[s]; ++k) {
          const int p = grammar->backpointers[(size_t)s *
                                              grammar->max_backpointers + k];
          value = logadd(value, previous[p] +
              transition_score(transition,
                              model_class(grammar->state_class[p]), class));
        }
        value += emission[class];
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

/** Plan O(states * (block + length / block)) capped partition workspace. */
int workspace_plan(size_t length, size_t states, size_t requested_block,
                        WorkspacePlan *plan) {
  if (!length || !states || !plan) return 0;
  size_t block = requested_block;
  if (!block) {
    size_t root = (size_t)sqrt((double)length);
    if (!root) root = 1;
    if (root < length / root + (length % root != 0)) ++root;
    block = CHECKPOINT_BLOCK > 0 ? (size_t)CHECKPOINT_BLOCK : 1;
    if (block < root) block = root;
  }
  if (block > length) block = length;
  const size_t blocks = (length - 1) / block + 1;
  const size_t max_cells = MAX_WORK_BYTES / sizeof(double);
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

/** Evaluate exact marginals, recomputing alpha blocks from saved checkpoints. */
int compute_gradient(const double *parameters, const LongGrammar *grammar,
                            const Sample *sample, int gold_only, double sign,
                            double *gradient, double *log_z_out,
                            size_t checkpoint_block) {
  const size_t states = (size_t)grammar->states;
  const size_t length = (size_t)sample->length;
  WorkspacePlan plan;
  if (!workspace_plan(length, states, checkpoint_block, &plan)) return 0;
  double *work = (double *)malloc(plan.bytes);
  if (!work) return -1;
  const size_t block = plan.block_length;
  double *checkpoints = work;
  double *alpha = checkpoints + (plan.blocks - 1) * states;
  double *beta = alpha + (block + 1) * states;
  double *previous_beta = beta + states;
  double *scales = previous_beta + states;
  const double *transition = parameters + EMISSIONS;

  /* The first pass needs only two alpha rows. Checkpoints are the normalized
   * row immediately before each block, not fresh initial distributions. */
  double *previous = alpha, *current = alpha + states;
  double log_scale = 0.0;
  for (size_t t = 0; t < length; ++t) {
    double scale;
    if (!forward_row(parameters, grammar, sample, gold_only, t,
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
    if (s==0)
      terminal = logadd(terminal, previous[s]);
  if (!isfinite(terminal) || !isfinite(log_scale + terminal)) goto nonfinite;
  *log_z_out = log_scale + terminal;
  for (size_t s = 0; s < states; ++s)
    beta[s] = s == 0 ? -terminal : -INFINITY;

  for (size_t b = plan.blocks; b-- > 0;) {
    const size_t start = b * block;
    const size_t width = length - start < block ? length - start : block;
    if (b)
      memcpy(alpha, checkpoints + (b - 1) * states, states * sizeof(*alpha));
    for (size_t offset = 0; offset < width; ++offset)
      if (!forward_row(parameters, grammar, sample, gold_only,
                             start + offset, alpha + offset * states,
                             alpha + (offset + 1) * states, &scales[offset]))
        goto nonfinite;

    for (size_t offset = width; offset-- > 0;) {
      const size_t t = start + offset;
      const double *row = alpha + (offset + 1) * states;
      double emission[CLASSES];
      for (int c = 0; c < CLASSES; ++c)
        emission[c] = gold_only && c != sample->classes[t]
                          ? 0.0 : sample_emission_scores(parameters, sample, t, c);
      for (size_t p = 0; p < states; ++p) previous_beta[p] = -INFINITY;
      for (size_t s = 0; s < states; ++s) {
        const int state_class = grammar->state_class[s];
        const int class = model_class(state_class);
        if (gold_only && state_class != sample->classes[t]) continue;
        const double probability = exp(row[s] + beta[s]);
        const size_t base = emission_index(class, 0);
        const uint8_t nucleotide = sample->bases[t];
        if (nucleotide < NUCLEOTIDES)
          gradient[base + nucleotide] += sign * probability;
        const float *cwt = sample->cwt + t * CWT_CHANNELS;
        for (int c = 0; c < CWT_CHANNELS; ++c)
          gradient[base + NUCLEOTIDES + (size_t)c] +=
              sign * probability * (double)cwt[c];
        gradient[base + INPUT_CHANNELS] += sign * probability;

        if (!t) continue;
        const double suffix = emission[class] + beta[s] - scales[offset];
        for (int k = 0; k < grammar->backpointer_count[s]; ++k) {
          const int p = grammar->backpointers[s *
                              (size_t)grammar->max_backpointers + (size_t)k];
          const int from_class =
              model_class(grammar->state_class[p]);
          const double edge = transition_score(transition, from_class, class);
          previous_beta[p] = logadd(previous_beta[p], edge + suffix);
          if (from_class != class) {
            const size_t index = EMISSIONS +
                                 (size_t)from_class * CLASSES + (size_t)class;
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

double model_objective(const double *x, double *gradient, void *opaque) {
  ModelObjectiveContext *context = (ModelObjectiveContext *)opaque;
  memset(gradient, 0, PARAMETERS * sizeof(*gradient));
  context->failure_status = 0;
  const size_t workers = (size_t)context->threads;
  if (!workers || workers > SIZE_MAX / PARAMETERS ||
      workers * PARAMETERS > SIZE_MAX / sizeof(double)) {
    context->failure_status = LBFGS_INVALID_ARGUMENT;
    return NAN;
  }
  double *thread_gradient =
      (double *)calloc(workers * PARAMETERS, sizeof(*thread_gradient));
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
  for (size_t i = 0; i < context->sample_count; ++i) {
    int thread = 0;
#ifdef _OPENMP
    thread = omp_get_thread_num();
#endif
    double *local_gradient =
        thread_gradient + (size_t)thread * PARAMETERS;
    const Sample *sample = &context->samples[i];
    if (!sample->length) continue;
    double log_z_all, log_z_gold;
    int all_status = compute_gradient(
        x, context->grammar, sample, 0, 1.0, local_gradient, &log_z_all, 0);
    int gold_status = all_status > 0
        ? compute_gradient(x, context->grammar, sample, 1, -1.0,
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
    for (size_t i = 0; i < PARAMETERS; ++i)
      gradient[i] += thread_gradient[thread * PARAMETERS + i];
  }
  free(thread_gradient);
  free(thread_objective);
  free(thread_status);
  if (context->failure_status)
    return NAN;
  const double inverse_total_length = 1.0 / (double)context->total_length;
  objective *= inverse_total_length;
  for (size_t i = 0; i < PARAMETERS; ++i) {
    gradient[i] *= inverse_total_length;
    objective += 0.5 * L2 * x[i] * x[i];
    gradient[i] += L2 * x[i];
  }
  if (!isfinite(objective)) {
    context->failure_status = LBFGS_NONFINITE;
    return NAN;
  }
  return objective;
}

/**
 * @brief Check that every model parameter is finite.
 * @param[in] model: Model to inspect.
 * @return Nonzero when all parameters are finite.
 */
int is_model_finite(const Model *model) {
  for (size_t i = 0; i < PARAMETERS; ++i) {
    const double value = i < EMISSIONS
                             ? model->emission_weights[i]
                             : model->transition_weights[i - EMISSIONS];
    if (!isfinite(value)) return 0;
  }
  return 1;
}

int save_model(const char *path, const Model *model) {
  static const char magic[4] = {'S', 'U', 'R', 'F'};
  if (!path || !model || !is_model_finite(model))
    return EINVAL;
  FILE *output = fopen(path, "wb");
  if (!output)
    return errno;
  uint32_t count = PARAMETERS;
  int ok = fwrite(magic, sizeof(magic), 1, output) == 1 &&
           fwrite(&count, sizeof(count), 1, output) == 1 &&
           fwrite(model, sizeof(*model), 1, output) == 1;
  int closed = fclose(output) == 0;
  return ok && closed ? 0 : EIO;
}

int load_model(const char *path, Model *model) {
  static const char magic[4] = {'S', 'U', 'R', 'F'};
  if (!path || !model)
    return EINVAL;
  FILE *input = fopen(path, "rb");
  if (!input)
    return errno;
  char found[sizeof(magic)];
  uint32_t count = 0;
  Model candidate;
  int ok = fread(found, sizeof(found), 1, input) == 1 &&
           fread(&count, sizeof(count), 1, input) == 1 &&
           !memcmp(found, magic, sizeof(magic)) && count == PARAMETERS &&
           fread(&candidate, sizeof(candidate), 1, input) == 1;
  if (ok) {
    int extra = fgetc(input);
    if (extra != EOF || ferror(input) || !is_model_finite(&candidate))
      ok = 0;
  }
  int closed = fclose(input) == 0;
  if (!ok || !closed)
    return EINVAL;
  *model = candidate;
  return 0;
}

#endif /** CRF_C */
