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

static Sample *prepare_samples(const Sequence *sequences, int sequence_count,
                               size_t *count) {
  Sample *samples = NULL;
  size_t capacity = 0;
  *count = 0;
  for (int s = 0; s < sequence_count; ++s) {
    const Sequence *sequence = &sequences[s];
    if (!sequence->labels || sequence->length < 1)
      bye("sequence %s has no labels", sequence->name);
    size_t length = (size_t)sequence->length;
    uint8_t *bases = scalloc(length, sizeof(*bases));
    float *cwt = scalloc(smul(length, CWT_CHANNELS), sizeof(*cwt));

    cwt_power(sequence->bases, sequence->length, 0, sequence->length, cwt);
    for (size_t p = 0; p < length; p++) {
      int base = base_index(sequence->bases[p]);
      bases[p] = (uint8_t)(base < 0 ? NUCLEOTIDES : base);
    }
    /* Condition windows on intergenic endpoints, not internal-exon background.
       Every base is retained, genes remain intact, and DNA context is preserved. */
    for (size_t begin = 0, end; begin < length; begin = end) {
      end = begin + 16384 < length ? begin + 16384 : length;
      size_t previous = end;
      while (previous > begin && sequence->labels[previous - 1] == INTERGENIC)
        --previous;
      uint8_t last = previous > begin ? sequence->labels[previous - 1] : INTERGENIC;
      while (end < length && (sequence->labels[end - 1] != INTERGENIC ||
                              sequence->labels[end] != INTERGENIC ||
                              last == PLUS_INTRON || last == MINUS_INTRON)) {
        if (sequence->labels[end] != INTERGENIC)
          last = sequence->labels[end];
        ++end;
      }
      if (*count == capacity) {
        capacity = capacity ? 2 * capacity : 64;
        samples = srealloc(samples, capacity, sizeof(*samples));
      }
      samples[(*count)++] = (Sample){
          (int64_t)(end - begin),   bases + begin,  cwt + begin * CWT_CHANNELS,
          sequence->labels + begin, (int64_t)begin, sequence->length, NULL, 0};
    }
  }
  return samples;
}

static size_t workspace_size(size_t length) {
  const size_t max_cells = SIZE_MAX / sizeof(double);
  if (!length ||
      length > (max_cells - 2 * CLASSES) / (2 * CLASSES + TRANSITIONS))
    return 0;
  return (length * (2 * CLASSES + TRANSITIONS) + 2 * CLASSES) * sizeof(double);
}

static int is_sample_valid(const Sample *sample) {
  if (sample->length < 0)
    return 0;
  if (sample->length == 0)
    return 1;
  if (!sample->bases || !sample->cwt || !sample->classes ||
      sample->offset < 0 ||
      (sample->sequence_length &&
       (sample->offset > sample->sequence_length ||
        sample->length > sample->sequence_length - sample->offset)) ||
      (!sample->sequence_length && sample->offset) ||
      (uint64_t)sample->length > SIZE_MAX / sizeof(float) / CWT_CHANNELS)
    return 0;
  if (!workspace_size((size_t)sample->length))
    return 0;
  for (int64_t position = 0; position < sample->length; ++position)
    if (sample->classes[position] >= CLASSES ||
        sample->bases[position] > NUCLEOTIDES)
      return 0;
  return 1;
}

int train(Model *model, const Sample *samples, size_t sample_count,
          int iterations, int threads) {
  if (!model || !samples || !sample_count || iterations < 1 || threads < 1 ||
      !is_model_finite(model))
    return LBFGS_INVALID_ARGUMENT;
  uint64_t total_length = 0;
  for (size_t i = 0; i < sample_count; ++i) {
    if (!is_sample_valid(&samples[i]))
      return LBFGS_INVALID_ARGUMENT;
    const uint64_t length = (uint64_t)samples[i].length;
    if (length > UINT64_MAX - total_length)
      return LBFGS_INVALID_ARGUMENT;
    total_length += length;
  }
  if (total_length == 0)
    return LBFGS_INVALID_ARGUMENT;

  Model candidate = *model;
  if ((size_t)threads > sample_count)
    threads = (int)sample_count;
  ModelObjectiveContext context = {samples, sample_count, total_length, threads, 0};
  int completed = 0;
  LbfgsStatus status =
      lbfgs_minimize(PARAMETERS, candidate.weights, model_objective, &context,
                     iterations, &completed);
  if (status == LBFGS_NONFINITE && context.failure_status)
    status = (LbfgsStatus)context.failure_status;
  if (status == LBFGS_MAX_ITERATIONS || status == LBFGS_CONVERGED ||
      status == LBFGS_STOPPING_CRITERION)
    *model = candidate;
  return (int)status;
}

static int boundary_features(const Sample *sample, size_t position,
                             size_t indices[2 * BOUNDARY_RADIUS + 1]) {
  position += sample->offset;
  size_t length = sample->sequence_length ? (size_t)sample->sequence_length
                                          : (size_t)sample->length;
  const uint8_t *bases = sample->bases - sample->offset;
  size_t begin = position > BOUNDARY_RADIUS ? position - BOUNDARY_RADIUS : 0;
  size_t end = position + BOUNDARY_RADIUS + 1;
  if (end > length)
    end = length;
  int count = 0;
  for (size_t p = begin; p < end; ++p)
    if (bases[p] < NUCLEOTIDES)
      indices[count++] =
          (p + BOUNDARY_RADIUS - position) * INPUT_CHANNELS + bases[p];
  return count;
}

static int cwt_boundary_slot(int slot) {
  int incoming = slot >= CLASSES, c = slot % CLASSES;
  return c == PLUS_INTRON || c == MINUS_INTRON ||
         (incoming ? c == PLUS_START_SINGLE_EXON || c == MINUS_END_EXON
                   : c == PLUS_END_EXON || c == MINUS_START_SINGLE_EXON);
}

static void boundary_cwt(const Sample *sample, size_t position,
                         float profile[BOUNDARY_FEATURES]) {
  memset(profile, 0, BOUNDARY_FEATURES * sizeof(*profile));
  position += sample->offset;
  size_t length = sample->sequence_length ? (size_t)sample->sequence_length
                                          : (size_t)sample->length;
  size_t begin = position > BOUNDARY_RADIUS ? position - BOUNDARY_RADIUS : 0;
  size_t end = position + BOUNDARY_RADIUS + 1;
  if (end > length)
    end = length;
  const float *power = sample->context_cwt
                           ? sample->context_cwt
                           : sample->cwt - sample->offset * CWT_CHANNELS;
  int64_t origin = sample->context_cwt ? sample->cwt_origin : 0;
  for (size_t p = begin; p < end; ++p)
    for (int k = 0; k < CWT_CHANNELS; ++k)
      profile[(p + BOUNDARY_RADIUS - position) * INPUT_CHANNELS +
              NUCLEOTIDES + k] = power[(p - origin) * CWT_CHANNELS + k];
}

Model *train_pipeline(const Sequence *sequences, int count, int iterations,
                      int threads) {
  double started = whattime();

  size_t sample_count;
  Sample *sample = prepare_samples(sequences, count, &sample_count);
  Model *model = scalloc(1, sizeof(*model));
  double *weights = model->weights;
#pragma omp parallel for num_threads(threads)                                  \
    reduction(+ : weights[ : PARAMETERS])
  for (size_t s = 0; s < sample_count; ++s)
    for (size_t t = 0; t < (size_t)sample[s].length; ++t) {
      int to = sample[s].classes[t];
      if (to >= CLASSES)
        bye("invalid CRF training label");
      if (sample[s].bases[t] < NUCLEOTIDES)
        weights[to * (INPUT_CHANNELS + 1) + sample[s].bases[t]]++;
      if (!t)
        continue;
      int from = sample[s].classes[t - 1];
      weights[EMISSIONS + from * CLASSES + to]++;
      if (from != to) {
        size_t indices[2 * BOUNDARY_RADIUS + 1];
        int n = boundary_features(&sample[s], t, indices);
        for (int i = 0; i < n; ++i) {
          weights[EMISSIONS + TRANSITIONS + from * BOUNDARY_FEATURES +
                  indices[i]]++;
          weights[EMISSIONS + TRANSITIONS + (CLASSES + to) * BOUNDARY_FEATURES +
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
  for (int from = 0; from < CLASSES; ++from)
    for (int to = 0; to < CLASSES; ++to)
      model->transition_weights[from * CLASSES + to] -=
          model->emission_weights[to * (INPUT_CHANNELS + 1) + INPUT_CHANNELS];
  for (size_t i = 0; i < BOUNDARY_WEIGHTS; i += INPUT_CHANNELS) {
    double total = NUCLEOTIDES;
    for (int b = 0; b < NUCLEOTIDES; ++b)
      total += model->boundary_weights[i + b];
    for (int b = 0; b < NUCLEOTIDES; ++b)
      model->boundary_weights[i + b] =
          log(NUCLEOTIDES * (model->boundary_weights[i + b] + 1) / total);
  }
  cry("Training: %d sequences, %zu windows, %d threads, %d iterations", count,
      sample_count, threads, iterations);
  int status = train(model, sample, sample_count, iterations, threads);
  for (size_t i = 0; i < sample_count; ++i)
    if (!sample[i].offset) {
      free((void *)sample[i].bases);
      free((void *)sample[i].cwt);
    }
  free(sample);
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


static double sample_emission_scores(const double *weights,
                                     const Sample *sample, size_t t,
                                     int class) {
  const size_t base = (size_t)class * (INPUT_CHANNELS + 1);
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
    const double *weights = model->emission_weights + c * (INPUT_CHANNELS + 1);
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
  float profile[BOUNDARY_FEATURES];
  boundary_cwt(sample, position, profile);
  double context[2 * CLASSES] = {0};
  for (int c = 0; c < 2 * CLASSES; ++c) {
    for (int i = 0; i < count; ++i)
      context[c] += pwm[c * BOUNDARY_FEATURES + indices[i]];
    if (cwt_boundary_slot(c))
      for (int p = 0; p < 2 * BOUNDARY_RADIUS + 1; ++p)
        for (int k = 0; k < CWT_CHANNELS; ++k) {
          int f = p * INPUT_CHANNELS + NUCLEOTIDES + k;
          context[c] += pwm[c * BOUNDARY_FEATURES + f] * profile[f];
        }
  }
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

static int forward_row(const double *emission, const double *transition,
                       const double *previous, double *row) {
  double top = -INFINITY;
  for (int s = 0; s < CLASSES; ++s) {
    double value = !previous && s == INTERGENIC ? 0 : -INFINITY;
    if (previous)
      for (int k = 0; k < gene_grammar.backpointers[s]; ++k) {
        int p = gene_grammar.from[s][k];
        value = logadd(value, previous[p] + transition[p * CLASSES + s]);
      }
    value += emission[s];
    if (value != -INFINITY && !isfinite(value))
      return 0;
    row[s] = value;
    if (value > top)
      top = value;
  }
  return isfinite(top);
}

static int sample_objective(const double *parameters, const Sample *sample,
                            double *gradient, double *objective) {
  size_t length = (size_t)sample->length;
  size_t bytes = workspace_size(length);
  if (!bytes || sample->classes[0] != INTERGENIC ||
      sample->classes[length - 1] != INTERGENIC)
    return 0;
  double *work = malloc(bytes);
  if (!work)
    return -1;
  double *alpha = work, *emissions = alpha + length * CLASSES;
  double *transitions = emissions + length * CLASSES;
  double *beta = transitions + length * TRANSITIONS;
  double *previous_beta = beta + CLASSES;
  double gold = 0;
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
        !forward_row(emission, transition,
                     t ? alpha + (t - 1) * CLASSES : NULL, alpha + t * CLASSES))
      goto nonfinite;
    int to = sample->classes[t];
    if (to >= CLASSES)
      goto nonfinite;
    gold += emission[to];
    if (t) {
      int from = sample->classes[t - 1], allowed = 0;
      for (int k = 0; k < gene_grammar.backpointers[to]; ++k)
        allowed |= gene_grammar.from[to][k] == from;
      if (!allowed)
        goto nonfinite;
      gold += transition[from * CLASSES + to];
    }
  }
  double log_z = alpha[(length - 1) * CLASSES + INTERGENIC];
  if (!isfinite(log_z) || !isfinite(gold))
    goto nonfinite;
  *objective = log_z - gold;
  for (int s = 0; s < CLASSES; ++s)
    beta[s] = s == INTERGENIC ? -log_z : -INFINITY;

  for (size_t t = length; t-- > 0;) {
    const double *emission = emissions + t * CLASSES;
    const double *transition = transitions + t * TRANSITIONS;
    double mass[CLASSES], edges[TRANSITIONS] = {0};
    for (int s = 0; s < CLASSES; ++s)
      previous_beta[s] = -INFINITY;
    for (int s = 0; s < CLASSES; ++s) {
      mass[s] = exp(alpha[t * CLASSES + s] + beta[s]);
      if (!t)
        continue;
      double suffix = emission[s] + beta[s];
      for (int k = 0; k < gene_grammar.backpointers[s]; ++k) {
        int from = gene_grammar.from[s][k];
        int pair = from * CLASSES + s;
        double value = transition[pair] + suffix;
        previous_beta[from] = logadd(previous_beta[from], value);
        if (from != s)
          edges[pair] += exp(alpha[(t - 1) * CLASSES + from] + value);
      }
    }
    for (int c = 0; c < CLASSES; ++c) {
      double residual = mass[c] - (sample->classes[t] == c);
      size_t base = (size_t)c * (INPUT_CHANNELS + 1);
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
      float profile[BOUNDARY_FEATURES];
      boundary_cwt(sample, t, profile);
      for (int c = 0; c < 2 * CLASSES; ++c)
        if (changes[c] != 0) {
          double *g = gradient + EMISSIONS + TRANSITIONS + c * BOUNDARY_FEATURES;
          for (int i = 0; i < count; ++i)
            g[indices[i]] += changes[c];
          if (cwt_boundary_slot(c))
            for (int p = 0; p < 2 * BOUNDARY_RADIUS + 1; ++p)
              for (int k = 0; k < CWT_CHANNELS; ++k) {
                int f = p * INPUT_CHANNELS + NUCLEOTIDES + k;
                g[f] += changes[c] * profile[f];
              }
        }
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
#pragma omp parallel for num_threads(context->threads) schedule(dynamic)
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
    int status = sample_objective(x, sample, worker->gradient, &value);
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
    for (size_t i = 0; i < PARAMETERS; ++i)
      gradient[i] += workers[thread].gradient[i] * inverse_length;
  }
  free(workers);
  for (size_t i = 0; i < PARAMETERS; ++i) {
    objective += 0.5 * L2 * x[i] * x[i];
    gradient[i] += L2 * x[i];
    if (!isfinite(gradient[i]))
      context->failure_status = LBFGS_NONFINITE;
  }
  if (!isfinite(objective))
    context->failure_status = LBFGS_NONFINITE;
  return context->failure_status ? NAN : objective;
}

int is_model_finite(const Model *model) {
  if (!model)
    return 0;
  for (size_t i = 0; i < PARAMETERS; ++i)
    if (!isfinite(model->weights[i]))
      return 0;
  return 1;
}

static const char model_magic[4] = {'S', 'U', 'R', 'F'};

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
