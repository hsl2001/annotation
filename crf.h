#ifndef CRF_H
#define CRF_H

#include "anno.h"
#include <stdint.h>

/* Linear-chain CRF decoded jointly for both strands.
   0 intergenic; 1-3 CDS(+) codon position; 4-6 intron(+) carrying the codon position of the
   next exon base; 7-9 CDS(-) codon position; 10-12 intron(-). */
enum { STATES = 13, HIDDEN_STATES = 47 };
enum { K2 = 17, K4 = 257, KNUC = 5, SPLICE_RADIUS = 12,
       FLANKS = SPLICE_RADIUS, EMISSION_FEATURES = 2, TRANSITION_FEATURES = 2 + FLANKS };
enum { OFF_ENUC = 0, OFF_E4 = OFF_ENUC + STATES * KNUC,
       OFF_EP = OFF_E4 + STATES * K4, OFF_ED = OFF_EP + STATES * WAVE_COUNT,
       OFF_TNUC = OFF_ED + STATES * WAVE_COUNT,
       OFF_T4 = OFF_TNUC + STATES * STATES * KNUC,
       OFF_TP = OFF_T4 + STATES * STATES * K4, OFF_TD = OFF_TP + STATES * STATES * WAVE_COUNT,
       OFF_SPLICE = OFF_TD + STATES * STATES * WAVE_COUNT,
       WEIGHTS = OFF_SPLICE + 4 * FLANKS * K2 };

void crf_init(void);
void crf_train(Contig *contigs, int count, double *weights, int iterations, int block, int threads);
void crf_decode(const double *weights, const Contig *contig, uint8_t *path, int block);

#endif

#ifdef CRF_IMPLEMENTATION

#include "lbfgs.h"
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  uint16_t k4[2];
  float cwt_power[WAVE_COUNT];
  int16_t derivative[WAVE_COUNT];
  uint8_t nuc, flank[2][FLANKS];
} Site;

static const int hidden_offsets[STATES] = {0, 1, 4, 8, 10, 11, 14, 18, 20, 28, 33, 41, 46};
static const int hidden_counts[STATES] = {1, 3, 4, 2, 1, 3, 4, 2, 8, 5, 8, 5, 1};
static int trans_a[STATES * STATES], trans_b[STATES * STATES], ntrans;
static uint8_t hidden_state[HIDDEN_STATES], hidden_context[HIDDEN_STATES];
typedef struct { uint8_t previous, current; uint16_t transition; } HiddenEdge;
static HiddenEdge hidden_edges[KNUC][HIDDEN_STATES * 4];
static int hidden_edge_count[KNUC];
static uint8_t initial_states[KNUC][HIDDEN_STATES];

static int kmer(const char *seq, int length, int from, int size, int reverse) {
  int code = 0;
  for (int offset = 0; offset < size; offset++) {
    int position = from + (reverse ? size - 1 - offset : offset);
    int base = position < 0 || position >= length ? -1 : base_index(seq[position]);
    if (base < 0) return 1 << (2 * size);
    code = code * 4 + (reverse ? 3 - base : base);
  }
  return code;
}

static Site site_at(const Contig *c, int i) {
  int base = base_index(c->seq[i]);
  Site s = {.k4 = {kmer(c->seq, c->length, i - 2, 4, 0),
                   kmer(c->seq, c->length, i - 2, 4, 1)}, .nuc = base < 0 ? 4 : base};
  for (int offset = 0; offset < FLANKS; offset++) {
    int displacement = -SPLICE_RADIUS + 2 * offset;
    s.flank[0][offset] = kmer(c->seq, c->length, i + displacement, 2, 0);
    s.flank[1][offset] = kmer(c->seq, c->length, i - displacement - 2, 2, 1);
  }
  memcpy(s.cwt_power, c->cwt_power + (size_t)i * WAVE_COUNT, sizeof(s.cwt_power));
  if (i > 0 && i + 1 < c->length && base_index(c->seq[i - 1]) >= 0 &&
      base >= 0 && base_index(c->seq[i + 1]) >= 0) {
    for (int scale = 0; scale < WAVE_COUNT; scale++) {
      double left = c->cwt_power[(size_t)(i - 1) * WAVE_COUNT + scale];
      double right = c->cwt_power[(size_t)(i + 1) * WAVE_COUNT + scale];
      double slope = (right - left) / (right + left + 1e-12);
      s.derivative[scale] = (int16_t)lrint(slope * INT16_MAX);
    }
  }
  return s;
}

static int boundary_kind(int previous, int current) {
  if ((previous >= 1 && previous <= 3 && current >= 4 && current <= 6) ||
      (previous >= 10 && current >= 7 && current <= 9)) return 0;
  if ((previous >= 4 && previous <= 6 && current >= 1 && current <= 3) ||
      (previous >= 7 && previous <= 9 && current >= 10)) return 1;
  if ((previous == 0 && current == 1) || (previous == 7 && current == 0)) return 2;
  if ((previous == 3 && current == 0) || (previous == 0 && current == 9)) return 3;
  return -1;
}

static int allowed(int a, int b) {
  if (a == 0) return b == 0 || b == 1 || b == 9;
  if (a <= 3) return b == 1 + a % 3 || b == 4 + a % 3 || (a == 3 && b == 0);
  if (a <= 6) return b == a || b == a - 3;
  if (a <= 9) return b == 7 + (a + 1) % 3 || b == 10 + (a + 1) % 3 || (a == 7 && b == 0);
  return b == a || b == a - 3;
}

/* Prefix classes carry start/stop recognition through introns without adding learned features. */
static int grammar_next(int previous, int current, int base) {
  int state = hidden_state[previous], context = hidden_context[previous];
  if (!allowed(state, current)) return -1;
  if (!state) {
    if (!current) return 0;
    if (current == 1) return base == 0 ? hidden_offsets[1] + 2 : -1;
    return base == 3 ? hidden_offsets[9] + 3 : base == 1 ? hidden_offsets[9] + 4 : -1;
  }
  if (!current)
    return (state == 3 || state == 7) && context == 1 ? 0 : -1;
  if (state <= 6) {
    if (current >= 4) {
      if (state == 3 && context == 1) return -1;
      return hidden_offsets[current] + (current == 4 ? 0 : context);
    }
    if (current == 1) {
      if (state == 3 && context == 1) return -1;
      return hidden_offsets[1] + (base == 3);
    }
    if (current == 2) {
      if (context == 2) return base == 3 ? hidden_offsets[2] + 3 : -1;
      int next_context = context == 1 && base == 0 ? 1 : context == 1 && base == 2 ? 2 : 0;
      return hidden_offsets[2] + next_context;
    }
    if (context == 3) return base == 2 ? hidden_offsets[3] : -1;
    int stop = (context == 1 && (base == 0 || base == 2)) || (context == 2 && base == 0);
    return hidden_offsets[3] + stop;
  }
  if (current >= 10) return hidden_offsets[current] + (current == 12 ? 0 : context);
  if (current == 9)
    return hidden_offsets[9] + (base == 3 ? 1 : base == 1 ? 2 : 0);
  if (current == 8) {
    if (context == 3)
      return base == 3 ? hidden_offsets[8] + 5 : base == 1 ? hidden_offsets[8] + 7 : -1;
    if (context == 4) return base == 3 ? hidden_offsets[8] + 6 : -1;
    int next_context = context == 1 && base == 3 ? 1 : context == 2 && base == 3 ? 2 :
                       context == 1 && base == 1 ? 3 : context == 2 && base == 0 ? 4 : 0;
    return hidden_offsets[8] + next_context;
  }
  if (context >= 5) return base == 0 ? hidden_offsets[7] : -1;
  if (context >= 1 && context <= 3 && base == 0) return -1;
  return hidden_offsets[7] + (context == 4 && base == 3);
}

static int initial_allowed(int base, int hidden, int full_start) {
  return full_start ? grammar_next(0, hidden_state[hidden], base) == hidden : initial_states[base][hidden];
}

static int terminal_allowed(int hidden) {
  int state = hidden_state[hidden];
  return !state || ((state == 3 || state == 7) && hidden_context[hidden] == 1);
}

void crf_init(void) {
  ntrans = 0;
  for (int a = 0; a < STATES; a++)
    for (int b = 0; b < STATES; b++)
      if (allowed(a, b)) trans_a[ntrans] = a, trans_b[ntrans++] = b;
  for (int state = 0; state < STATES; state++)
    for (int context = 0; context < hidden_counts[state]; context++) {
      int hidden = hidden_offsets[state] + context;
      hidden_state[hidden] = (uint8_t)state;
      hidden_context[hidden] = (uint8_t)context;
    }
  memset(hidden_edge_count, 0, sizeof(hidden_edge_count));
  memset(initial_states, 0, sizeof(initial_states));
  for (int base = 0; base < KNUC; base++)
    for (int previous = 0; previous < HIDDEN_STATES; previous++)
      for (int edge = 0; edge < ntrans; edge++) {
        if (trans_a[edge] != hidden_state[previous]) continue;
        int current = grammar_next(previous, trans_b[edge], base);
        if (current < 0) continue;
        if (hidden_edge_count[base] >= HIDDEN_STATES * 4) anno_fail("Too many grammar transitions");
        hidden_edges[base][hidden_edge_count[base]++] = (HiddenEdge){(uint8_t)previous, (uint8_t)current, (uint16_t)edge};
        initial_states[base][current] = 1;
      }
}

static void emission_indices(const Site *site, int state, int *indices) {
  indices[0] = OFF_ENUC + state * KNUC + site->nuc;
  indices[1] = OFF_E4 + state * K4 + site->k4[state >= 7];
}

static int transition_indices(const Site *site, int previous, int current, int *indices) {
  int pair = previous * STATES + current, count = 2;
  int reverse = previous >= 7 || current >= 7, kind = boundary_kind(previous, current);
  indices[0] = OFF_TNUC + pair * KNUC + site->nuc;
  indices[1] = OFF_T4 + pair * K4 + site->k4[reverse];
  if (kind >= 0)
    for (int offset = 0; offset < FLANKS; offset++)
      indices[count++] = OFF_SPLICE + (kind * FLANKS + offset) * K2 + site->flank[reverse][offset];
  return count;
}

static double feature_score(const double *weights, const int *indices, int count) {
  double score = 0;
  for (int index = 0; index < count; index++) score += weights[indices[index]];
  return score;
}

static float derivative_value(const Site *site, int scale) {
  return site->derivative[scale] / (float)INT16_MAX;
}

static void add_features(double *gradient, const int *indices, int count, double value) {
  for (int index = 0; index < count; index++) gradient[indices[index]] += value;
}

static double emission(const double *weights, const Site *site, int state) {
  int indices[EMISSION_FEATURES];
  emission_indices(site, state, indices);
  double score = feature_score(weights, indices, EMISSION_FEATURES);
  for (int scale = 0; scale < WAVE_COUNT; scale++) {
    score += weights[OFF_EP + state * WAVE_COUNT + scale] * site->cwt_power[scale];
    score += weights[OFF_ED + state * WAVE_COUNT + scale] * derivative_value(site, scale);
  }
  return score;
}

static double transition(const double *weights, const Site *site, int previous, int current) {
  int indices[TRANSITION_FEATURES];
  int count = transition_indices(site, previous, current, indices);
  int pair = previous * STATES + current;
  double score = feature_score(weights, indices, count);
  for (int scale = 0; scale < WAVE_COUNT; scale++) {
    score += weights[OFF_TP + pair * WAVE_COUNT + scale] * site->cwt_power[scale];
    score += weights[OFF_TD + pair * WAVE_COUNT + scale] * derivative_value(site, scale);
  }
  return score;
}

static void add_emission(double *gradient, const Site *site, int state, double value) {
  int indices[EMISSION_FEATURES];
  emission_indices(site, state, indices);
  add_features(gradient, indices, EMISSION_FEATURES, value);
  for (int scale = 0; scale < WAVE_COUNT; scale++) {
    gradient[OFF_EP + state * WAVE_COUNT + scale] += value * site->cwt_power[scale];
    gradient[OFF_ED + state * WAVE_COUNT + scale] += value * derivative_value(site, scale);
  }
}

static void add_transition(double *gradient, const Site *site, int previous, int current, double value) {
  int indices[TRANSITION_FEATURES];
  int count = transition_indices(site, previous, current, indices);
  add_features(gradient, indices, count, value);
  int pair = previous * STATES + current;
  for (int scale = 0; scale < WAVE_COUNT; scale++) {
    gradient[OFF_TP + pair * WAVE_COUNT + scale] += value * site->cwt_power[scale];
    gradient[OFF_TD + pair * WAVE_COUNT + scale] += value * derivative_value(site, scale);
  }
}

typedef struct { Site *sites; double *emis, *psi, *alpha, *beta, *scale; } Work;

static double log_add(double left, double right) {
  if (left == -INFINITY) return right;
  if (right == -INFINITY) return left;
  double maximum = fmax(left, right);
  return maximum + log1p(exp(-fabs(left - right)));
}

/* Normalized log-space forward-backward on [s, e). */
static double window_posteriors(const double *w, const Contig *c, int s, int e, Work *k) {
  int len = e - s;
  double logz = 0;
  for (int i = 0; i < len; i++) {
    k->sites[i] = site_at(c, s + i);
    for (int st = 0; st < STATES; st++) {
      k->emis[i * STATES + st] = emission(w, k->sites + i, st);
      if (!isfinite(k->emis[i * STATES + st])) anno_fail("Non-finite CRF emission at position %d", s + i + 1);
    }
  }
  for (int i = 0; i < len; i++) {
    const Site *site = k->sites + i;
    const double *em = k->emis + i * STATES;
    double *al = k->alpha + (size_t)i * HIDDEN_STATES, *psi = k->psi + (size_t)i * ntrans;
    if (!i) {
      for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
        al[hidden] = initial_allowed(site->nuc, hidden, s == 0) ? em[hidden_state[hidden]] : -INFINITY;
    } else {
      for (int t = 0; t < ntrans; t++) {
        psi[t] = transition(w, site, trans_a[t], trans_b[t]) + em[trans_b[t]];
        if (!isfinite(psi[t])) anno_fail("Non-finite CRF transition at position %d", s + i + 1);
      }
      for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) al[hidden] = -INFINITY;
      const HiddenEdge *edges = hidden_edges[site->nuc];
      for (int edge = 0; edge < hidden_edge_count[site->nuc]; edge++) {
        const HiddenEdge *step = edges + edge;
        al[step->current] = log_add(al[step->current], al[(int)step->previous - HIDDEN_STATES] + psi[step->transition]);
      }
    }
    if (e == c->length && i == len - 1)
      for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
        if (!terminal_allowed(hidden)) al[hidden] = -INFINITY;
    double total = -INFINITY;
    for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) total = log_add(total, al[hidden]);
    if (!isfinite(total)) anno_fail("Invalid CRF partition at position %d", s + i + 1);
    for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) al[hidden] -= total;
    k->scale[i] = total;
    logz += total;
  }
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
    k->beta[(size_t)(len - 1) * HIDDEN_STATES + hidden] = e < c->length || terminal_allowed(hidden) ? 0 : -INFINITY;
  for (int i = len - 2; i >= 0; i--) {
    double *be = k->beta + (size_t)i * HIDDEN_STATES, *next = be + HIDDEN_STATES, *psi = k->psi + (size_t)(i + 1) * ntrans;
    int base = k->sites[i + 1].nuc;
    const HiddenEdge *edges = hidden_edges[base];
    for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) be[hidden] = -INFINITY;
    for (int edge = 0; edge < hidden_edge_count[base]; edge++) {
      const HiddenEdge *step = edges + edge;
      be[step->previous] = log_add(be[step->previous], psi[step->transition] + next[step->current]);
    }
    for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) be[hidden] -= k->scale[i + 1];
  }
  return logz;
}

static double gold_grammar_mass(const Contig *contig, int start, int end) {
  double previous[HIDDEN_STATES], current[HIDDEN_STATES];
  int base = base_index(contig->seq[start]);
  if (base < 0) base = 4;
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
    previous[hidden] = hidden_state[hidden] == contig->label[start] && initial_allowed(base, hidden, start == 0) ? 0 : -INFINITY;
  for (int position = start + 1; position < end; position++) {
    base = base_index(contig->seq[position]);
    if (base < 0) base = 4;
    for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) current[hidden] = -INFINITY;
    const HiddenEdge *edges = hidden_edges[base];
    for (int edge = 0; edge < hidden_edge_count[base]; edge++) {
      const HiddenEdge *step = edges + edge;
      if (hidden_state[step->current] == contig->label[position])
        current[step->current] = log_add(current[step->current], previous[step->previous]);
    }
    memcpy(previous, current, sizeof(previous));
  }
  double mass = -INFINITY;
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
    if (end < contig->length || terminal_allowed(hidden)) mass = log_add(mass, previous[hidden]);
  return mass;
}

/* Accumulates observed minus expected features for one window.
   Returns the window log-likelihood or NAN when the gold path violates the grammar. */
static double window_gradient(const double *w, const Contig *c, int s, int e, Work *k, double *gradient) {
  int len = e - s;
  double mass = gold_grammar_mass(c, s, e);
  if (mass == -INFINITY) return NAN;
  double gold = mass - window_posteriors(w, c, s, e, k);
  for (int i = 0; i < len; i++) {
    const Site *site = k->sites + i;
    const double *al = k->alpha + (size_t)i * HIDDEN_STATES, *be = k->beta + (size_t)i * HIDDEN_STATES, *psi = k->psi + (size_t)i * ntrans;
    int b = c->label[s + i];
    gold += k->emis[i * STATES + b];
    add_emission(gradient, site, b, 1);
    double state_marginals[STATES] = {0};
    for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
      state_marginals[hidden_state[hidden]] += exp(al[hidden] + be[hidden]);
    for (int st = 0; st < STATES; st++) add_emission(gradient, site, st, -state_marginals[st]);
    if (!i) continue;
    gold += transition(w, site, c->label[s + i - 1], b);
    add_transition(gradient, site, c->label[s + i - 1], b, 1);
    double edge_marginals[STATES * STATES] = {0};
    const HiddenEdge *edges = hidden_edges[site->nuc];
    for (int edge = 0; edge < hidden_edge_count[site->nuc]; edge++) {
      const HiddenEdge *step = edges + edge;
      edge_marginals[step->transition] += exp(al[(int)step->previous - HIDDEN_STATES] + psi[step->transition] + be[step->current] - k->scale[i]);
    }
    for (int t = 0; t < ntrans; t++)
      if (edge_marginals[t] > 0) add_transition(gradient, site, trans_a[t], trans_b[t], -edge_marginals[t]);
  }
  return gold;
}

static int window_count(const Contig *c, int block) {
  return (int)(((long)c->length + block - 1) / block);
}

static Work work_alloc(int block) {
  return (Work){
    .sites = anno_alloc((size_t)block, sizeof(Site)),
    .emis = anno_alloc((size_t)block * STATES, sizeof(double)),
    .psi = anno_alloc((size_t)block * ntrans, sizeof(double)),
    .alpha = anno_alloc((size_t)block * HIDDEN_STATES, sizeof(double)),
    .beta = anno_alloc((size_t)block * HIDDEN_STATES, sizeof(double)),
    .scale = anno_alloc((size_t)block, sizeof(double))
  };
}

static void work_free(Work *k) {
  free(k->sites), free(k->emis), free(k->psi), free(k->alpha), free(k->beta), free(k->scale);
}

/* CRF++-style batch training: gradients of all windows are evaluated in parallel
   at a fixed weight vector, summed, and passed to L-BFGS once per iteration. */
typedef struct {
  Contig *contigs;
  int block, threads, windows;
  long positions;
  int skipped;
} CrfTrain;

typedef struct {
  const CrfTrain *train;
  const double *w;
  int first, last;
  double likelihood;
  long positions;
  int skipped;
  double *gradient;
} CrfWorker;

static void *crf_worker(void *arg) {
  CrfWorker *worker = arg;
  const CrfTrain *train = worker->train;
  Work work = work_alloc(train->block);
  worker->gradient = anno_alloc(WEIGHTS, sizeof(*worker->gradient));
  for (int window = worker->first; window < worker->last; window++) {
    int contig = 0, index = window;
    while (index >= window_count(train->contigs + contig, train->block))
      index -= window_count(train->contigs + contig++, train->block);
    Contig *c = train->contigs + contig;
    int start = index * train->block;
    int end = c->length - start < train->block ? c->length : start + train->block;
    double score = window_gradient(worker->w, c, start, end, &work, worker->gradient);
    if (isnan(score)) { worker->skipped++; continue; }
    if (!isfinite(score)) anno_fail("Non-finite training likelihood");
    worker->likelihood += score;
    worker->positions += end - start;
  }
  work_free(&work);
  return NULL;
}

static double crf_evaluate(const double *w, double *gradient, void *context) {
  CrfTrain *train = context;
  CrfWorker *workers = anno_alloc((size_t)train->threads, sizeof(*workers));
  pthread_t *tid = anno_alloc((size_t)train->threads, sizeof(*tid));
  memset(gradient, 0, WEIGHTS * sizeof(*gradient));
  train->positions = 0;
  train->skipped = 0;
  double likelihood = 0;
  for (int t = 0; t < train->threads; t++) {
    workers[t] = (CrfWorker){.train = train, .w = w,
                             .first = (int)((long)train->windows * t / train->threads),
                             .last = (int)((long)train->windows * (t + 1) / train->threads)};
    if (train->threads > 1 && pthread_create(tid + t, NULL, crf_worker, workers + t))
      anno_fail("Cannot spawn training thread");
  }
  for (int t = 0; t < train->threads; t++) {
    if (train->threads > 1) pthread_join(tid[t], NULL); else crf_worker(workers + t);
    likelihood += workers[t].likelihood;
    train->positions += workers[t].positions;
    train->skipped += workers[t].skipped;
    for (int f = 0; f < WEIGHTS; f++) gradient[f] -= workers[t].gradient[f]; /* minimize -loglik */
    free(workers[t].gradient);
  }
  free(workers);
  free(tid);
  return -likelihood;
}

static void crf_progress(void *context, int iteration, double objective, size_t n, const double *gradient) {
  const CrfTrain *train = context;
  double norm = 0;
  for (size_t i = 0; i < n; i++) norm += gradient[i] * gradient[i];
  fprintf(stderr, "Iteration %d: log-likelihood/bp %.4f, |gradient| %.4g, %d windows skipped\n",
          iteration, train->positions ? -objective / train->positions : 0.0, sqrt(norm), train->skipped);
}

void crf_train(Contig *contigs, int count, double *w, int iterations, int block, int threads) {
  CrfTrain train = {.contigs = contigs, .block = block, .threads = threads < 1 ? 1 : threads};
  for (int c = 0; c < count; c++) train.windows += window_count(contigs + c, block);
  if (!train.windows) anno_fail("No training windows");
  if (train.threads > train.windows) train.threads = train.windows;
  if (lbfgs_minimize(WEIGHTS, w, crf_evaluate, crf_progress, &train, iterations, 10) < 0)
    anno_fail("Out of memory");
}

static void viterbi_step(const double *weights, const Contig *contig, int position,
                         const double *previous, double *current, uint8_t *back) {
  Site site = site_at(contig, position);
  double emissions[STATES], transitions[STATES * STATES];
  for (int state = 0; state < STATES; state++) {
    emissions[state] = emission(weights, &site, state);
    if (!isfinite(emissions[state])) anno_fail("Non-finite prediction emission");
  }
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) {
    current[hidden] = position == 0 && initial_allowed(site.nuc, hidden, 1) ? emissions[hidden_state[hidden]] : -INFINITY;
    if (back) back[hidden] = 0;
  }
  if (position) {
    for (int edge = 0; edge < ntrans; edge++) {
      transitions[edge] = transition(weights, &site, trans_a[edge], trans_b[edge]) + emissions[trans_b[edge]];
      if (!isfinite(transitions[edge])) anno_fail("Non-finite prediction transition");
    }
    const HiddenEdge *edges = hidden_edges[site.nuc];
    for (int edge = 0; edge < hidden_edge_count[site.nuc]; edge++) {
      const HiddenEdge *step = edges + edge;
      double score = previous[step->previous] + transitions[step->transition];
      if (score > current[step->current]) {
        current[step->current] = score;
        if (back) back[step->current] = step->previous;
      }
    }
  }
  double maximum = -INFINITY;
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) maximum = fmax(maximum, current[hidden]);
  if (!isfinite(maximum)) anno_fail("No finite prediction path");
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) current[hidden] -= maximum;
}

/* Exact whole-contig MAP; recomputed checkpoint blocks bound traceback memory. */
void crf_decode(const double *w, const Contig *c, uint8_t *path, int block) {
  int length = c->length, blocks = window_count(c, block);
  double *checkpoints = anno_alloc((size_t)blocks * HIDDEN_STATES, sizeof(*checkpoints));
  int capacity = length < block ? length : block;
  uint8_t *back = anno_alloc((size_t)capacity * HIDDEN_STATES, sizeof(*back));
  double previous[HIDDEN_STATES] = {0}, current[HIDDEN_STATES];
  for (int position = 0; position < length; position++) {
    if (position % block == 0)
      memcpy(checkpoints + (size_t)(position / block) * HIDDEN_STATES, previous, sizeof(previous));
    viterbi_step(w, c, position, previous, current, NULL);
    memcpy(previous, current, sizeof(previous));
  }
  int best = 0;
  for (int hidden = 1; hidden < HIDDEN_STATES; hidden++)
    if (terminal_allowed(hidden) && previous[hidden] > previous[best]) best = hidden;
  if (!isfinite(previous[best])) anno_fail("No complete prediction path");
  for (int chunk = blocks - 1; chunk >= 0; chunk--) {
    int start = chunk * block, end = length - start < block ? length : start + block;
    memcpy(previous, checkpoints + (size_t)chunk * HIDDEN_STATES, sizeof(previous));
    for (int position = start; position < end; position++) {
      viterbi_step(w, c, position, previous, current, back + (size_t)(position - start) * HIDDEN_STATES);
      memcpy(previous, current, sizeof(previous));
    }
    for (int position = end - 1; position >= start; position--) {
      path[position] = hidden_state[best];
      best = back[(size_t)(position - start) * HIDDEN_STATES + best];
    }
  }
  free(checkpoints);
  free(back);
}

#endif
