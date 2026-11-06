/**
 * @file trainer.c
 * @brief Linear-chain CRF training and model serialization.
 */
#include "include/crf.h"
#include "include/lbfgs.h"
#include "include/utils.h"

#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static Sample *prepare_samples(const Sequence *sequences, int sequence_count) {
  if (sequence_count < 1)
    bye("no sequences to train CRF");
  Sample *samples = scalloc((size_t)sequence_count, sizeof(*samples));
  for (int s = 0; s < sequence_count; ++s) {
    const Sequence *sequence = &sequences[s];
    if (!sequence->labels || sequence->length < 1)
      bye("sequence %s has no labels", sequence->name);
    size_t length = (size_t)sequence->length;
    uint8_t *bases = scalloc(length, sizeof(*bases));
    float *cwt = scalloc(smul(length, CWT_CHANNELS), sizeof(*cwt));
    samples[s] = (Sample){sequence->length, bases, cwt, sequence->labels};
    cwt_power(sequence->bases, sequence->length, 0, sequence->length, cwt);
    for (size_t p = 0; p < length; p++) {
      int base = base_index(sequence->bases[p]);
      bases[p] = (uint8_t)(base < 0 ? NUCLEOTIDES : base);
    }
  }

  return samples;
}

static int is_grammar_valid(const Grammar *grammar) {
  if (!grammar || grammar->count <= 0 || grammar->count > GRAMMAR_MAX_STATES)
    return 0;
  for (int state = 0; state < grammar->count; ++state) {
    if (grammar->cls[state] >= CLASSES || grammar->backpointers[state] < 0 ||
        grammar->backpointers[state] > GRAMMAR_MAX_BACKPOINTERS)
      return 0;
    for (int edge = 0; edge < grammar->backpointers[state]; ++edge) {
      int predecessor = grammar->from[state][edge];
      if (predecessor < 0 || predecessor >= grammar->count)
        return 0;
    }
  }
  return 1;
}

static void prepare_samples_destroy(Sample *samples, size_t count) {
  if (!samples)
    return;
  for (size_t i = 0; i < count; ++i) {
    free((void *)samples[i].bases);
    free((void *)samples[i].cwt);
  }
  free(samples);
}

static int is_sample_valid(const Sample *sample, const Grammar *grammar) {
  if (sample->length < 0)
    return 0;
  if (sample->length == 0)
    return 1;
  if (!sample->bases || !sample->cwt || !sample->classes ||
      (uint64_t)sample->length > SIZE_MAX / sizeof(float) / CWT_CHANNELS)
    return 0;
  WorkspacePlan plan;
  if (!workspace_plan((size_t)sample->length, (size_t)grammar->count, 0, &plan))
    return 0;
  for (int64_t position = 0; position < sample->length; ++position)
    if (sample->classes[position] >= CLASSES ||
        sample->bases[position] > NUCLEOTIDES)
      return 0;
  return 1;
}

int train(Model *model, const Grammar *grammar, const Sample *samples,
          size_t sample_count, int iterations, int threads) {
  if (!model || !is_grammar_valid(grammar) || !samples || !sample_count ||
      iterations < 1 || threads < 1 || !is_model_finite(model))
    return LBFGS_INVALID_ARGUMENT;
  uint64_t total_length = 0;
  for (size_t i = 0; i < sample_count; ++i) {
    if (!is_sample_valid(&samples[i], grammar))
      return LBFGS_INVALID_ARGUMENT;
    const uint64_t length = (uint64_t)samples[i].length;
    if (length > UINT64_MAX - total_length)
      return LBFGS_INVALID_ARGUMENT;
    total_length += length;
  }
  if (total_length == 0)
    return LBFGS_INVALID_ARGUMENT;

  double parameters[PARAMETERS];
  memcpy(parameters, model->emission_weights, sizeof(model->emission_weights));
  memcpy(parameters + EMISSIONS, model->transition_weights,
         sizeof(model->transition_weights));
  memcpy(parameters + EMISSIONS + TRANSITIONS, model->boundary_weights,
         sizeof(model->boundary_weights));
  if ((size_t)threads > sample_count)
    threads = (int)sample_count;
  ModelObjectiveContext context = {grammar,      samples, sample_count,
                                   total_length, threads, 0};
  int completed = 0;
  LbfgsStatus status = lbfgs_minimize(PARAMETERS, parameters, model_objective,
                                      NULL, &context, iterations, &completed);
  if (status == LBFGS_NONFINITE && context.failure_status)
    status = (LbfgsStatus)context.failure_status;
  if (status == LBFGS_MAX_ITERATIONS || status == LBFGS_CONVERGED ||
      status == LBFGS_STOPPING_CRITERION) {
    memcpy(model->emission_weights, parameters,
           sizeof(model->emission_weights));
    memcpy(model->transition_weights, parameters + EMISSIONS,
           sizeof(model->transition_weights));
    memcpy(model->boundary_weights, parameters + EMISSIONS + TRANSITIONS,
           sizeof(model->boundary_weights));
  }
  return (int)status;
}

static int boundary_features(const Sample *sample, size_t position,
                             size_t indices[2 * BOUNDARY_RADIUS + 1]) {
  size_t begin = position > BOUNDARY_RADIUS ? position - BOUNDARY_RADIUS : 0;
  size_t end = (size_t)sample->length - position > BOUNDARY_RADIUS
                   ? position + BOUNDARY_RADIUS + 1
                   : (size_t)sample->length;
  int count = 0;
  for (size_t p = begin; p < end; ++p)
    if (sample->bases[p] < NUCLEOTIDES)
      indices[count++] =
          (p + BOUNDARY_RADIUS - position) * NUCLEOTIDES + sample->bases[p];
  return count;
}

Model *train_pipeline(const Sequence *sequences, int count, int iterations,
                      int threads) {
  Grammar grammar;
  build_grammar(&grammar);
  Sample *sample = prepare_samples(sequences, count);
  Model *model = scalloc(1, sizeof(*model));
  double started = whattime();
  /* Accumulate sufficient statistics in the zero-initialized model. */
  for (int s = 0; s < count; ++s)
    for (size_t t = 0; t < (size_t)sample[s].length; ++t) {
      int to = sample[s].classes[t];
      if (to >= CLASSES)
        bye("invalid CRF training label");
      if (sample[s].bases[t] < NUCLEOTIDES)
        model->emission_weights[to * (INPUT_CHANNELS + 1) +
                                sample[s].bases[t]]++;
      if (!t)
        continue;
      int from = sample[s].classes[t - 1];
      model->transition_weights[from * CLASSES + to]++;
      if (from != to) {
        size_t indices[2 * BOUNDARY_RADIUS + 1];
        int n = boundary_features(&sample[s], t, indices);
        for (int i = 0; i < n; ++i) {
          model->boundary_weights[from * BOUNDARY_FEATURES + indices[i]]++;
          model->boundary_weights[(CLASSES + to) * BOUNDARY_FEATURES +
                                  indices[i]]++;
        }
      }
    }
  for (int c = 0; c < CLASSES; ++c) {
    double *emission = model->emission_weights + c * (INPUT_CHANNELS + 1);
    double *transition = model->transition_weights + c * CLASSES;
    double bases = NUCLEOTIDES, outgoing = CLASSES;
    for (int b = 0; b < NUCLEOTIDES; ++b)
      bases += emission[b];
    for (int to = 0; to < CLASSES; ++to)
      outgoing += transition[to];
    for (int b = 0; b < NUCLEOTIDES; ++b)
      emission[b] = log((emission[b] + 1) / bases);
    for (int to = 0; to < CLASSES; ++to)
      transition[to] = log((transition[to] + 1) / outgoing);
    emission[INPUT_CHANNELS] = transition[c];
  }
  /* Each change cancels the destination's stay score already in its bias. */
  for (int from = 0; from < CLASSES; ++from)
    for (int to = 0; to < CLASSES; ++to)
      model->transition_weights[from * CLASSES + to] -=
          model->emission_weights[to * (INPUT_CHANNELS + 1) + INPUT_CHANNELS];
  for (size_t i = 0; i < BOUNDARY_WEIGHTS; i += NUCLEOTIDES) {
    double total = NUCLEOTIDES;
    for (int b = 0; b < NUCLEOTIDES; ++b)
      total += model->boundary_weights[i + b];
    for (int b = 0; b < NUCLEOTIDES; ++b)
      model->boundary_weights[i + b] =
          log(NUCLEOTIDES * (model->boundary_weights[i + b] + 1) / total);
  }
  cry("Training: %d sequences, %d iterations", count, iterations);
  int status =
      train(model, &grammar, sample, (size_t)count, iterations, threads);
  prepare_samples_destroy(sample, (size_t)count);
  if (status != LBFGS_MAX_ITERATIONS && status != LBFGS_CONVERGED &&
      status != LBFGS_STOPPING_CRITERION)
    bye("CRF L-BFGS training failed (status %d)", status);
  cry("CRF L-BFGS status %d after %.1fs", status, whattime() - started);
  return model;
}

static double logadd(double a, double b) {
  if (a == -INFINITY)
    return b;
  if (b == -INFINITY)
    return a;
  if (a < b) {
    double t = a;
    a = b;
    b = t;
  }
  return a + log1p(exp(b - a));
}

static size_t emission_index(int class, int channel) {
  return (size_t)class * (INPUT_CHANNELS + 1) + (size_t)channel;
}

static double sample_emission_scores(const double *weights,
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

void compute_emissions(const Model *model, const float *features,
                       double scores[CLASSES]) {
  if (!model || !features || !scores)
    return;
  for (int c = 0; c < CLASSES; ++c) {
    const double *weights = model->emission_weights + emission_index(c, 0);
    scores[c] = weights[INPUT_CHANNELS];
    for (int k = 0; k < INPUT_CHANNELS; ++k)
      scores[c] += weights[k] * features[k];
  }
}

int compute_transitions(const double *bias, const double *pwm,
                        const Sample *sample, size_t position,
                        double scores[TRANSITIONS]) {
  size_t indices[2 * BOUNDARY_RADIUS + 1];
  int count = boundary_features(sample, position, indices);
  double context[2 * CLASSES] = {0};
  for (int c = 0; c < 2 * CLASSES; ++c)
    for (int i = 0; i < count; ++i)
      context[c] += pwm[c * BOUNDARY_FEATURES + indices[i]];
  for (int from = 0; from < CLASSES; ++from)
    for (int to = 0; to < CLASSES; ++to) {
      int pair = from * CLASSES + to;
      scores[pair] =
          from == to ? 0 : bias[pair] + context[from] + context[CLASSES + to];
      if (!isfinite(scores[pair]))
        return 0;
    }
  return 1;
}

int workspace_plan(size_t length, size_t states, size_t requested_block,
                   WorkspacePlan *plan) {
  if (!length || !states || !plan)
    return 0;
  const size_t max_cells = SIZE_MAX / sizeof(double);
  if (states > max_cells / 2 ||
      length > (max_cells - 2 * states) / (states + CLASSES + TRANSITIONS))
    return 0;
  size_t block = requested_block ? requested_block : length;
  if (block > length)
    block = length;
  *plan =
      (WorkspacePlan){block, (length - 1) / block + 1,
                      (length * (states + CLASSES + TRANSITIONS) + 2 * states) *
                          sizeof(double)};
  return 1;
}

static int forward_row(const Grammar *grammar, const double *emission,
                       const double *transition, const double *previous,
                       double *row, int observed, int initial) {
  double top = -INFINITY;
  for (int s = 0; s < grammar->count; ++s) {
    int class = grammar->cls[s];
    double value = -INFINITY;
    if (observed < 0 || class == observed) {
      if (initial) {
        if (s == 0)
          value = 0;
      } else {
        for (int k = 0; k < grammar->backpointers[s]; ++k) {
          int p = grammar->from[s][k];
          value =
              logadd(value, previous[p] +
                                transition[grammar->cls[p] * CLASSES + class]);
        }
      }
      value += emission[class];
    }
    if (value != -INFINITY && !isfinite(value))
      return 0;
    row[s] = value;
    if (value > top)
      top = value;
  }
  return isfinite(top);
}

static int sample_objective(const double *parameters, const Grammar *grammar,
                            const Sample *sample, double *gradient,
                            double *objective) {
  size_t states = (size_t)grammar->count, length = (size_t)sample->length;
  WorkspacePlan plan;
  if (!workspace_plan(length, states, 0, &plan))
    return 0;
  double *work = malloc(plan.bytes);
  if (!work)
    return -1;
  double *alpha = work, *emissions = alpha + length * states;
  double *transitions = emissions + length * CLASSES;
  double *gold_previous = transitions + length * TRANSITIONS;
  double *gold_current = gold_previous + states;
  for (size_t t = 0; t < length; ++t) {
    double *emission = emissions + t * CLASSES;
    double *transition = transitions + t * TRANSITIONS;
    for (int c = 0; c < CLASSES; ++c) {
      emission[c] = sample_emission_scores(parameters, sample, t, c);
      if (!isfinite(emission[c]))
        goto nonfinite;
    }
    if (!compute_transitions(parameters + EMISSIONS,
                             parameters + EMISSIONS + TRANSITIONS, sample, t,
                             transition) ||
        !forward_row(grammar, emission, transition,
                     t ? alpha + (t - 1) * states : NULL, alpha + t * states,
                     -1, t == 0) ||
        !forward_row(grammar, emission, transition, gold_previous, gold_current,
                     sample->classes[t], t == 0))
      goto nonfinite;
    double *swap = gold_previous;
    gold_previous = gold_current;
    gold_current = swap;
  }
  double log_z = alpha[(length - 1) * states];
  if (!isfinite(log_z) || !isfinite(gold_previous[0]))
    goto nonfinite;
  *objective = log_z - gold_previous[0];
  double *beta = gold_previous, *previous_beta = gold_current;
  for (size_t s = 0; s < states; ++s)
    beta[s] = s == 0 ? -log_z : -INFINITY;

  /* All gold-compatible state paths have identical class features. Their
   * forward sum retains multiplicity; its derivative is the observed counts. */
  for (size_t t = length; t-- > 0;) {
    const double *emission = emissions + t * CLASSES;
    const double *transition = transitions + t * TRANSITIONS;
    double mass[CLASSES] = {0}, edges[TRANSITIONS] = {0};
    for (size_t s = 0; s < states; ++s)
      previous_beta[s] = -INFINITY;
    for (size_t s = 0; s < states; ++s) {
      int class = grammar->cls[s];
      mass[class] += exp(alpha[t * states + s] + beta[s]);
      if (!t)
        continue;
      double suffix = emission[class] + beta[s];
      for (int k = 0; k < grammar->backpointers[s]; ++k) {
        int p = grammar->from[s][k], from = grammar->cls[p];
        size_t pair = (size_t)from * CLASSES + (size_t)class;
        double value = transition[pair] + suffix;
        previous_beta[p] = logadd(previous_beta[p], value);
        if (from != class)
          edges[pair] += exp(alpha[(t - 1) * states + (size_t)p] + value);
      }
    }
    for (int c = 0; c < CLASSES; ++c) {
      double residual = mass[c] - (sample->classes[t] == c);
      size_t base = emission_index(c, 0);
      if (sample->bases[t] < NUCLEOTIDES)
        gradient[base + sample->bases[t]] += residual;
      for (int k = 0; k < CWT_CHANNELS; ++k)
        gradient[base + NUCLEOTIDES + (size_t)k] +=
            residual * sample->cwt[t * CWT_CHANNELS + (size_t)k];
      gradient[base + INPUT_CHANNELS] += residual;
    }
    if (t) {
      int from = sample->classes[t - 1], to = sample->classes[t];
      if (from != to)
        edges[from * CLASSES + to] -= 1;
      size_t indices[2 * BOUNDARY_RADIUS + 1];
      int count = boundary_features(sample, t, indices);
      double changes[2 * CLASSES] = {0};
      for (int pair = 0; pair < TRANSITIONS; ++pair) {
        gradient[EMISSIONS + pair] += edges[pair];
        changes[pair / CLASSES] += edges[pair];
        changes[CLASSES + pair % CLASSES] += edges[pair];
      }
      for (int c = 0; c < 2 * CLASSES; ++c)
        if (changes[c] != 0)
          for (int i = 0; i < count; ++i)
            gradient[EMISSIONS + TRANSITIONS + c * BOUNDARY_FEATURES +
                     indices[i]] += changes[c];
    }
    double *swap = beta;
    beta = previous_beta;
    previous_beta = swap;
  }
  free(work);
  return 1;

nonfinite:
  free(work);
  return 0;
}

double model_objective(const double *x, double *gradient, void *opaque) {
  typedef struct {
    double gradient[PARAMETERS], objective;
    int status;
  } Worker;
  ModelObjectiveContext *context = opaque;
  memset(gradient, 0, PARAMETERS * sizeof(*gradient));
  context->failure_status = 0;
  size_t count = (size_t)context->threads;
  if (!count || count > SIZE_MAX / sizeof(Worker) || !context->total_length) {
    context->failure_status = LBFGS_INVALID_ARGUMENT;
    return NAN;
  }
  for (size_t i = 0; i < PARAMETERS; ++i)
    if (!isfinite(x[i])) {
      context->failure_status = LBFGS_NONFINITE;
      return NAN;
    }
  Worker *workers = calloc(count, sizeof(*workers));
  if (!workers) {
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
    Worker *worker = &workers[thread];
    const Sample *sample = &context->samples[i];
    if (!sample->length)
      continue;
    double value;
    int status =
        sample_objective(x, context->grammar, sample, worker->gradient, &value);
    if (status > 0)
      worker->objective += value;
    else
      worker->status = status < 0 ? LBFGS_ALLOCATION_FAILURE : LBFGS_NONFINITE;
  }
  double objective = 0, inverse_length = 1.0 / (double)context->total_length;
  for (size_t thread = 0; thread < count; ++thread) {
    if (workers[thread].status)
      context->failure_status = workers[thread].status;
    objective += workers[thread].objective * inverse_length;
    for (size_t i = 0; i < PARAMETERS; ++i) {
      gradient[i] += workers[thread].gradient[i] * inverse_length;
      if (!isfinite(gradient[i]))
        context->failure_status = LBFGS_NONFINITE;
    }
  }
  free(workers);
  if (!isfinite(objective))
    context->failure_status = LBFGS_NONFINITE;
  return context->failure_status ? NAN : objective;
}

int is_model_finite(const Model *model) {
  if (!model)
    return 0;
  const double *weights[] = {model->emission_weights, model->transition_weights,
                             model->boundary_weights};
  const size_t counts[] = {EMISSIONS, TRANSITIONS, BOUNDARY_WEIGHTS};
  for (size_t group = 0; group < 3; ++group)
    for (size_t i = 0; i < counts[group]; ++i)
      if (!isfinite(weights[group][i]))
        return 0;
  return 1;
}

static const char model_magic[4] = {'S', 'U', 'R', '4'};

int save_model(const char *path, const Model *model) {
  if (!path || !model || !is_model_finite(model))
    return EINVAL;
  FILE *output = fopen(path, "wb");
  if (!output)
    return errno;
  uint32_t count = PARAMETERS;
  int ok = fwrite(model_magic, sizeof(model_magic), 1, output) == 1 &&
           fwrite(&count, sizeof(count), 1, output) == 1 &&
           fwrite(model, sizeof(*model), 1, output) == 1;
  int closed = fclose(output) == 0;
  return ok && closed ? 0 : EIO;
}

int load_model(const char *path, Model *model) {
  if (!path || !model)
    return EINVAL;
  FILE *input = fopen(path, "rb");
  if (!input)
    return errno;
  char found[sizeof(model_magic)];
  uint32_t count = 0;
  Model candidate;
  int ok = fread(found, sizeof(found), 1, input) == 1 &&
           fread(&count, sizeof(count), 1, input) == 1 &&
           !memcmp(found, model_magic, sizeof(model_magic)) &&
           count == PARAMETERS &&
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
