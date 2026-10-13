#ifndef CRF_H
#define CRF_H

#include "anno.h"
#include "lbfgs.h"
#include <stdint.h>

#define CRF_L2_C2 2.0

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
void crf_train(Contig *contigs, int count, double *weights, int iterations, int threads);
void crf_decode(const double *weights, const Contig *contig, uint8_t *path);

#endif

#ifdef CRF_IMPLEMENTATION

#include <math.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CRF_WINDOW_SIZE 65536

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
  int position = c->offset + i;
  int base = base_index(c->seq[position]);
  Site s = {.k4 = {kmer(c->seq, c->genome_length, position - 2, 4, 0),
                   kmer(c->seq, c->genome_length, position - 2, 4, 1)}, .nuc = base < 0 ? 4 : base};
  for (int offset = 0; offset < FLANKS; offset++) {
    int displacement = -SPLICE_RADIUS + 2 * offset;
    s.flank[0][offset] = kmer(c->seq, c->genome_length, position + displacement, 2, 0);
    s.flank[1][offset] = kmer(c->seq, c->genome_length, position - displacement - 2, 2, 1);
  }
  memcpy(s.cwt_power, c->cwt_power + (size_t)position * WAVE_COUNT, sizeof(s.cwt_power));
  if (position > 0 && position + 1 < c->genome_length &&
      base_index(c->seq[position - 1]) >= 0 && base >= 0 &&
      base_index(c->seq[position + 1]) >= 0) {
    for (int scale = 0; scale < WAVE_COUNT; scale++) {
      double left = c->cwt_power[(size_t)(position - 1) * WAVE_COUNT + scale];
      double right = c->cwt_power[(size_t)(position + 1) * WAVE_COUNT + scale];
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

typedef struct { double *checkpoint, *alpha, *scale; int block_size; } Work;

static double log_add(double left, double right) {
  if (left == -INFINITY) return right;
  if (right == -INFINITY) return left;
  double maximum = fmax(left, right);
  return maximum + log1p(exp(-fabs(left - right)));
}

static Site site_scores(const double *w, const Contig *c, int position,
                        double *emissions, double *psi) {
  Site site = site_at(c, position);
  for (int state = 0; state < STATES; state++) {
    emissions[state] = emission(w, &site, state);
    if (!isfinite(emissions[state]))
      anno_fail("Non-finite CRF emission at position %d", position + 1);
  }
  for (int edge = 0; edge < ntrans; edge++) {
    psi[edge] = transition(w, &site, trans_a[edge], trans_b[edge]) + emissions[trans_b[edge]];
    if (!isfinite(psi[edge]))
      anno_fail("Non-finite CRF transition at position %d", position + 1);
  }
  return site;
}

static double forward_step(const double *w, const Contig *c, int position,
                           const double *previous, double *alpha) {
  Site site;
  double emissions[STATES], psi[STATES * STATES];
  site = site_scores(w, c, position, emissions, psi);
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) alpha[hidden] = -INFINITY;
  if (!position) {
    for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
      if (initial_allowed(site.nuc, hidden, c->offset == 0))
        alpha[hidden] = emissions[hidden_state[hidden]];
  } else {
    const HiddenEdge *edges = hidden_edges[site.nuc];
    for (int edge = 0; edge < hidden_edge_count[site.nuc]; edge++) {
      const HiddenEdge *step = edges + edge;
      alpha[step->current] = log_add(alpha[step->current],
                                     previous[step->previous] + psi[step->transition]);
    }
  }
  if (position == c->length - 1 &&
      c->offset + c->length == c->genome_length)
    for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
      if (!terminal_allowed(hidden)) alpha[hidden] = -INFINITY;

  double total = -INFINITY;
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) total = log_add(total, alpha[hidden]);
  if (!isfinite(total)) anno_fail("Invalid CRF partition at position %d", position + 1);
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) alpha[hidden] -= total;
  return total;
}

static void backward_step(const Site *next_site, const double *next_psi,
                          const double *next_beta, double next_scale, double *beta) {
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) beta[hidden] = -INFINITY;
  const HiddenEdge *edges = hidden_edges[next_site->nuc];
  for (int edge = 0; edge < hidden_edge_count[next_site->nuc]; edge++) {
    const HiddenEdge *step = edges + edge;
    beta[step->previous] = log_add(beta[step->previous],
                                   next_psi[step->transition] + next_beta[step->current]);
  }
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) beta[hidden] -= next_scale;
}

static double gold_grammar_mass(const Contig *contig, int start, int end) {
  double previous[HIDDEN_STATES], current[HIDDEN_STATES];
  int position = contig->offset + start;
  int base = base_index(contig->seq[position]);
  if (base < 0) base = 4;
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
    previous[hidden] = hidden_state[hidden] == contig->label[position] &&
      initial_allowed(base, hidden, contig->offset == 0) ? 0 : -INFINITY;
  for (int position = start + 1; position < end; position++) {
    int genome_position = contig->offset + position;
    base = base_index(contig->seq[genome_position]);
    if (base < 0) base = 4;
    for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) current[hidden] = -INFINITY;
    const HiddenEdge *edges = hidden_edges[base];
    for (int edge = 0; edge < hidden_edge_count[base]; edge++) {
      const HiddenEdge *step = edges + edge;
      if (hidden_state[step->current] == contig->label[genome_position])
        current[step->current] = log_add(current[step->current], previous[step->previous]);
    }
    memcpy(previous, current, sizeof(previous));
  }
  double mass = -INFINITY;
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
    if (contig->offset + end < contig->genome_length || terminal_allowed(hidden))
      mass = log_add(mass, previous[hidden]);
  return mass;
}

/* Accumulates the exact observed-minus-expected gradient for one training window. */
static double contig_gradient(const double *w, const Contig *c, Work *k, double *gradient) {
  double mass = gold_grammar_mass(c, 0, c->length);
  if (mass == -INFINITY) return NAN;

  int block_size = k->block_size;
  double logz = 0, alpha_a[HIDDEN_STATES], alpha_b[HIDDEN_STATES];
  const double *alpha_previous = NULL;
  double *alpha_current = alpha_a;
  for (int position = 0; position < c->length; position++) {
    logz += forward_step(w, c, position, alpha_previous, alpha_current);
    if (position + 1 < c->length && (position + 1) % block_size == 0)
      memcpy(k->checkpoint + (size_t)(position / block_size) * HIDDEN_STATES,
             alpha_current, sizeof(alpha_a));
    alpha_previous = alpha_current;
    alpha_current = alpha_current == alpha_a ? alpha_b : alpha_a;
  }

  double likelihood = mass - logz;
  double beta[HIDDEN_STATES];
  int reaches_contig_end = c->offset + c->length == c->genome_length;
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
    beta[hidden] = !reaches_contig_end || terminal_allowed(hidden) ? 0 : -INFINITY;

  for (int block_start = ((c->length - 1) / block_size) * block_size;
       block_start >= 0;) {
    int block_length = c->length - block_start;
    if (block_length > block_size) block_length = block_size;
    int block_end = block_start + block_length;
    const double *checkpoint = block_start
      ? k->checkpoint + (size_t)(block_start / block_size - 1) * HIDDEN_STATES : NULL;
    for (int position = block_start; position < block_end; position++) {
      int local = position - block_start;
      double *alpha = k->alpha + (size_t)local * HIDDEN_STATES;
      const double *before = position == block_start ? checkpoint : alpha - HIDDEN_STATES;
      k->scale[local] = forward_step(w, c, position, before, alpha);
    }
    for (int position = block_end - 1; position >= block_start; position--) {
      int local = position - block_start;
      Site site;
      double emissions[STATES], psi[STATES * STATES];
      site = site_scores(w, c, position, emissions, psi);
      const double *alpha = k->alpha + (size_t)local * HIDDEN_STATES;
      int state = c->label[c->offset + position];
      likelihood += emissions[state];
      add_emission(gradient, &site, state, 1);

      double state_marginals[STATES] = {0};
      for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
        state_marginals[hidden_state[hidden]] += exp(alpha[hidden] + beta[hidden]);
      for (int st = 0; st < STATES; st++)
        add_emission(gradient, &site, st, -state_marginals[st]);

      if (position) {
        const double *previous_alpha = local
          ? alpha - HIDDEN_STATES : checkpoint;
        int previous_state = c->label[c->offset + position - 1];
        likelihood += transition(w, &site, previous_state, state);
        add_transition(gradient, &site, previous_state, state, 1);

        double edge_marginals[STATES * STATES] = {0};
        const HiddenEdge *edges = hidden_edges[site.nuc];
        for (int edge = 0; edge < hidden_edge_count[site.nuc]; edge++) {
          const HiddenEdge *step = edges + edge;
          edge_marginals[step->transition] += exp(
            previous_alpha[step->previous] + psi[step->transition] +
            beta[step->current] - k->scale[local]);
        }
        for (int edge = 0; edge < ntrans; edge++)
          if (edge_marginals[edge] > 0)
            add_transition(gradient, &site, trans_a[edge], trans_b[edge], -edge_marginals[edge]);

        double beta_previous[HIDDEN_STATES];
        backward_step(&site, psi, beta, k->scale[local], beta_previous);
        memcpy(beta, beta_previous, sizeof(beta));
      }
    }
    if (!block_start) break;
    block_start -= block_size;
  }
  return likelihood;
}

static Work work_alloc(int length) {
  int block_size = (int)ceil(sqrt((double)length));
  int checkpoint_count = (length - 1) / block_size;
  return (Work){
    .checkpoint = anno_alloc((size_t)checkpoint_count, HIDDEN_STATES * sizeof(double)),
    .alpha = anno_alloc((size_t)block_size, HIDDEN_STATES * sizeof(double)),
    .scale = anno_alloc((size_t)block_size, sizeof(double)),
    .block_size = block_size
  };
}

static void work_free(Work *k) {
  free(k->checkpoint), free(k->alpha), free(k->scale);
}

typedef struct {
  Contig *contig;
  int start, length;
} CrfWindow;

static int next_window_end(const Contig *contig, int start) {
  int end = contig->length - start > CRF_WINDOW_SIZE ? start + CRF_WINDOW_SIZE : contig->length;
  if (end == contig->length || !contig->label[end - 1] || !contig->label[end]) return end;
  int gene_start = end - 1, gene_end = end;
  while (gene_start > start && contig->label[gene_start - 1]) gene_start--;
  while (gene_end < contig->length && contig->label[gene_end]) gene_end++;
  return gene_start > start && end - gene_start <= gene_end - end ? gene_start : gene_end;
}

static size_t make_windows(Contig *contigs, int count, CrfWindow **windows) {
  size_t window_count = 0, capacity = 0;
  for (int i = 0; i < count; i++) {
    Contig *contig = contigs + i;
    for (int start = 0; start < contig->length;) {
      int end = next_window_end(contig, start);
      if (end <= start) anno_fail("Cannot split contig into CRF windows: %s", contig->name);
      if (window_count == capacity) {
        if (capacity > SIZE_MAX / 2 / sizeof(**windows))
          anno_fail("Too many CRF windows");
        size_t next_capacity = capacity ? capacity * 2 : 16;
        CrfWindow *grown = realloc(*windows, next_capacity * sizeof(**windows));
        if (!grown) anno_fail("Out of memory");
        *windows = grown;
        capacity = next_capacity;
      }
      (*windows)[window_count++] = (CrfWindow){contig, start, end - start};
      start = end;
    }
  }
  return window_count;
}

typedef struct {
  CrfWindow *windows;
  size_t count;
  int threads;
  atomic_size_t next_instance;
} CrfTrain;

typedef struct {
  CrfTrain *train;
  const double *w;
  double likelihood;
  size_t positions, skipped;
  double *gradient;
} CrfWorker;

static void *crf_worker(void *arg) {
  CrfWorker *worker = arg;
  CrfTrain *train = worker->train;
  worker->likelihood = 0;
  worker->positions = 0;
  worker->skipped = 0;
  memset(worker->gradient, 0, WEIGHTS * sizeof(*worker->gradient));
  size_t slot = atomic_fetch_add_explicit(&train->next_instance, 1, memory_order_relaxed);
  while (slot < train->count) {
    const CrfWindow *window = train->windows + slot;
    Contig instance = *window->contig;
    instance.offset = window->start;
    instance.length = window->length;
    Work work = work_alloc(instance.length);
    double score = contig_gradient(worker->w, &instance, &work, worker->gradient);
    work_free(&work);
    if (isnan(score)) {
      worker->skipped++;
    } else {
      if (!isfinite(score)) anno_fail("Non-finite training likelihood");
      worker->likelihood += score;
      worker->positions += (size_t)instance.length;
    }
    slot = atomic_fetch_add_explicit(&train->next_instance, 1, memory_order_relaxed);
  }
  return NULL;
}

typedef struct {
  CrfTrain train;
  CrfWorker *workers;
  pthread_t *threads;
  double likelihood;
  size_t positions, skipped;
  double *scale, *weights;
} CrfObjective;

static double crf_objective(const double *w, double *gradient, void *context) {
  CrfObjective *objective = context;
  CrfTrain *train = &objective->train;
  atomic_store_explicit(&train->next_instance, 0, memory_order_relaxed);
  memset(gradient, 0, WEIGHTS * sizeof(*gradient));
  for (int t = 0; t < train->threads; t++) {
    objective->workers[t].w = w;
    if (train->threads > 1 &&
        pthread_create(objective->threads + t, NULL, crf_worker, objective->workers + t))
      anno_fail("Cannot spawn training thread");
  }

  double likelihood = 0;
  size_t positions = 0, skipped = 0;
  for (int t = 0; t < train->threads; t++) {
    if (train->threads > 1 && pthread_join(objective->threads[t], NULL))
      anno_fail("Cannot join training thread");
    if (train->threads == 1) crf_worker(objective->workers + t);
    likelihood += objective->workers[t].likelihood;
    positions += objective->workers[t].positions;
    skipped += objective->workers[t].skipped;
    for (int f = 0; f < WEIGHTS; f++) gradient[f] += objective->workers[t].gradient[f];
  }
  if (!positions) anno_fail("No valid CRF training windows (%zu windows skipped)", skipped);
  if (!isfinite(likelihood)) anno_fail("Non-finite training likelihood");

  objective->likelihood = likelihood;
  objective->positions = positions;
  objective->skipped = skipped;
  double penalty = 0;
  for (int f = 0; f < WEIGHTS; f++) {
    penalty += w[f] * w[f];
    gradient[f] = -gradient[f] + 2 * CRF_L2_C2 * w[f];
    if (!isfinite(gradient[f])) anno_fail("Non-finite training gradient");
  }
  double loss = -likelihood + CRF_L2_C2 * penalty;
  if (!isfinite(loss)) anno_fail("Non-finite training objective");
  return loss;
}

static void crf_progress(void *context, int iteration, double objective, double step,
                         int line_search_trials, size_t n, const double *gradient) {
  const CrfObjective *evaluation = context;
  double gradient_norm = 0;
  for (size_t f = 0; f < n; f++) gradient_norm = hypot(gradient_norm, gradient[f]);
  fprintf(stderr,
          "Iteration %d: %zu windows trained, %zu skipped, log-likelihood/bp %.8g, "
          "objective %.8g, gradient-norm %.3g, step %.3g, line-search %d\n",
          iteration, evaluation->train.count - evaluation->skipped, evaluation->skipped,
          evaluation->likelihood / evaluation->positions, objective, gradient_norm,
          step, line_search_trials);
}

/* L-BFGS optimizes v with w = scale * v, so frequent dense and rare sparse features share curvature. */
static double crf_scaled_objective(const double *v, double *gradient, void *context) {
  CrfObjective *objective = context;
  for (int f = 0; f < WEIGHTS; f++) objective->weights[f] = objective->scale[f] * v[f];
  double loss = crf_objective(objective->weights, gradient, context);
  for (int f = 0; f < WEIGHTS; f++) gradient[f] *= objective->scale[f];
  return loss;
}

/* scale_f = 1 / sqrt(max(1, sum of squared gold-path values of feature f)). */
static void feature_scales(const Contig *contigs, int count, double *scale) {
  memset(scale, 0, WEIGHTS * sizeof(*scale));
  for (int c = 0; c < count; c++)
    for (int i = 0; i < contigs[c].length; i++) {
      Site site = site_at(contigs + c, i);
      int state = contigs[c].label[i], indices[TRANSITION_FEATURES];
      emission_indices(&site, state, indices);
      for (int j = 0; j < EMISSION_FEATURES; j++) scale[indices[j]] += 1;
      for (int k = 0; k < WAVE_COUNT; k++) {
        double derivative = derivative_value(&site, k);
        scale[OFF_EP + state * WAVE_COUNT + k] += (double)site.cwt_power[k] * site.cwt_power[k];
        scale[OFF_ED + state * WAVE_COUNT + k] += derivative * derivative;
      }
      if (!i) continue;
      int previous = contigs[c].label[i - 1], pair = previous * STATES + state;
      int n = transition_indices(&site, previous, state, indices);
      for (int j = 0; j < n; j++) scale[indices[j]] += 1;
      for (int k = 0; k < WAVE_COUNT; k++) {
        double derivative = derivative_value(&site, k);
        scale[OFF_TP + pair * WAVE_COUNT + k] += (double)site.cwt_power[k] * site.cwt_power[k];
        scale[OFF_TD + pair * WAVE_COUNT + k] += derivative * derivative;
      }
    }
  for (int f = 0; f < WEIGHTS; f++) scale[f] = 1 / sqrt(fmax(scale[f], 1));
}

void crf_train(Contig *contigs, int count, double *w, int iterations, int threads) {
  if (count < 1 || iterations < 1)
    anno_fail("CRF training requires instances and positive iterations");
  if (threads < 1) anno_fail("CRF training requires at least one worker thread");
  CrfWindow *windows = NULL;
  size_t window_count = make_windows(contigs, count, &windows);
  if (!window_count) anno_fail("CRF training requires instances");
  CrfObjective objective = {
    .train = {.windows = windows, .count = window_count, .threads = threads}
  };
  atomic_init(&objective.train.next_instance, 0);
  fprintf(stderr, "Using %d worker threads for %zu training windows\n",
          objective.train.threads, objective.train.count);
  objective.workers = anno_alloc((size_t)objective.train.threads, sizeof(*objective.workers));
  objective.threads = anno_alloc((size_t)objective.train.threads, sizeof(*objective.threads));
  for (int t = 0; t < objective.train.threads; t++)
    objective.workers[t] = (CrfWorker){
      .train = &objective.train,
      .gradient = anno_alloc(WEIGHTS, sizeof(*objective.workers->gradient))
    };

  int completed = 0;
  objective.scale = anno_alloc(WEIGHTS, sizeof(*objective.scale));
  objective.weights = anno_alloc(WEIGHTS, sizeof(*objective.weights));
  feature_scales(contigs, count, objective.scale);
  for (int f = 0; f < WEIGHTS; f++) w[f] /= objective.scale[f];
  LbfgsStatus status = lbfgs_minimize(WEIGHTS, w, crf_scaled_objective, crf_progress, &objective,
                                      iterations, &completed);
  for (int f = 0; f < WEIGHTS; f++) w[f] *= objective.scale[f];
  free(objective.scale), free(objective.weights);
  for (int t = 0; t < objective.train.threads; t++) free(objective.workers[t].gradient);
  free(objective.workers), free(objective.threads), free(windows);
  if (status == LBFGS_ALLOCATION_FAILURE) anno_fail("Out of memory in L-BFGS optimizer");
  if (status == LBFGS_NONFINITE) anno_fail("Non-finite value in L-BFGS optimizer");
  if (status == LBFGS_LINE_SEARCH_FAILURE)
    anno_fail("L-BFGS line search failed after %d iterations", completed);
  if (status == LBFGS_INVALID_ARGUMENT) anno_fail("Invalid L-BFGS optimizer configuration");
  if (status == LBFGS_CONVERGED)
    fprintf(stderr, "L-BFGS converged after %d iterations\n", completed);
  if (status == LBFGS_STOPPING_CRITERION)
    fprintf(stderr, "L-BFGS stopped by objective history after %d iterations\n", completed);
  if (status == LBFGS_MAX_ITERATIONS)
    fprintf(stderr, "L-BFGS reached the %d-iteration limit\n", completed);
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

/* Exact whole-contig MAP decoding. */
void crf_decode(const double *w, const Contig *c, uint8_t *path) {
  int length = c->length;
  uint8_t *back = anno_alloc((size_t)length, HIDDEN_STATES * sizeof(*back));
  double previous[HIDDEN_STATES] = {0}, current[HIDDEN_STATES];
  for (int position = 0; position < length; position++) {
    viterbi_step(w, c, position, previous, current,
                 back + (size_t)position * HIDDEN_STATES);
    memcpy(previous, current, sizeof(previous));
  }
  int best = 0;
  for (int hidden = 1; hidden < HIDDEN_STATES; hidden++)
    if (terminal_allowed(hidden) && previous[hidden] > previous[best]) best = hidden;
  if (!isfinite(previous[best])) anno_fail("No complete prediction path");
  for (int position = length - 1; position >= 0; position--) {
    path[position] = hidden_state[best];
    best = back[(size_t)position * HIDDEN_STATES + best];
  }
  free(back);
}

#endif
