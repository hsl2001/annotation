#ifndef CRF_H
#define CRF_H

#include "anno.h"
#include "lbfgs.h"
#include <stdint.h>

#ifndef CRF_L2_C2
#define CRF_L2_C2 1.0
#endif
#ifndef CRF_WINDOW_SIZE
#define CRF_WINDOW_SIZE 65536
#endif

/* Linear-chain CRF decoded jointly for both strands.
   0 intergenic; 1-3 CDS(+) codon position; 4-6 intron(+) carrying the codon position of the
   next exon base; 7-9 CDS(-) codon position; 10-12 intron(-). */
enum { STATES = 13, HIDDEN_STATES = 47 };
enum { K2 = 17, K4 = 257, KNUC = 5, RADIUS = 8,
       FLANKS = RADIUS, EMISSION_FEATURES = 2, TRANSITION_FEATURES = 2 + FLANKS };
enum { OFF_ENUC = 0, OFF_E4 = OFF_ENUC + STATES * KNUC,
       OFF_EP = OFF_E4 + STATES * K4, OFF_ED = OFF_EP + STATES * WAVE_COUNT,
       OFF_TNUC = OFF_ED + STATES * WAVE_COUNT,
       OFF_T4 = OFF_TNUC + STATES * STATES * KNUC,
       OFF_TP = OFF_T4 + STATES * STATES * K4, OFF_TD = OFF_TP + STATES * STATES * WAVE_COUNT,
       OFF_SPLICE = OFF_TD + STATES * STATES * WAVE_COUNT,
       OFF_BP = OFF_SPLICE + 4 * FLANKS * K2,
       OFF_BD = OFF_BP + 4 * FLANKS * WAVE_COUNT,
       WEIGHTS = OFF_BD + 4 * FLANKS * WAVE_COUNT };

void crf_init(void);
void crf_train(Contig *contigs, int count, double *weights, int iterations, int threads);
void crf_decode(const double *weights, const Contig *contig, uint8_t *path);

#endif

#ifdef CRF_IMPLEMENTATION

#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  uint16_t k4_emit[2], k4_trans[2];
  float cwt_power[WAVE_COUNT];
  int16_t derivative[WAVE_COUNT];
  uint8_t nuc, flank[2][FLANKS];
} Site;

static const int hidden_offsets[STATES] = {0, 1, 4, 8, 10, 11, 14, 18, 20, 28, 33, 41, 46};
static const int hidden_counts[STATES] = {1, 3, 4, 2, 1, 3, 4, 2, 8, 5, 8, 5, 1};
static int trans_a[STATES * STATES], trans_b[STATES * STATES], ntrans;
static int edge_pair[STATES * STATES], edge_rev[STATES * STATES], edge_kind[STATES * STATES];
static int pair_to_edge[STATES][STATES];
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
  Site s = {
    .k4_emit = {kmer(c->seq, c->genome_length, position - 2, 4, 0),
                kmer(c->seq, c->genome_length, position - 1, 4, 1)},
    .k4_trans = {kmer(c->seq, c->genome_length, position - 2, 4, 0),
                 kmer(c->seq, c->genome_length, position - 2, 4, 1)},
    .nuc = base < 0 ? 4 : base
  };
  for (int offset = 0; offset < FLANKS; offset++) {
    int displacement = -RADIUS + 2 * offset;
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

static void boundary_cwt_at(const Contig *c, int position, int reverse, int offset,
                            float *power, double *derivative) {
  static const int offsets[FLANKS] = {-27, -26, -25, -24, -2, -1, 0, 1};
  int displacement = offsets[offset];
  int g = reverse ? position - displacement - 1 : position + displacement;
  if (g < 0 || g >= c->genome_length || base_index(c->seq[g]) < 0) {
    for (int scale = 0; scale < WAVE_COUNT; scale++) power[scale] = 0, derivative[scale] = 0;
    return;
  }
  memcpy(power, c->cwt_power + (size_t)g * WAVE_COUNT, WAVE_COUNT * sizeof(*power));
  if (g >= 2 && g + 2 < c->genome_length &&
      base_index(c->seq[g - 2]) >= 0 && base_index(c->seq[g + 2]) >= 0) {
    for (int scale = 0; scale < WAVE_COUNT; scale++) {
      double center = power[scale];
      double left = c->cwt_power[(size_t)(g - 2) * WAVE_COUNT + scale];
      double right = c->cwt_power[(size_t)(g + 2) * WAVE_COUNT + scale];
      derivative[scale] = (2.0 * center - left - right) / (2.0 * center + left + right + 1e-12);
    }
  } else {
    for (int scale = 0; scale < WAVE_COUNT; scale++) derivative[scale] = 0;
  }
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
  memset(pair_to_edge, -1, sizeof(pair_to_edge));
  for (int a = 0; a < STATES; a++)
    for (int b = 0; b < STATES; b++)
      if (allowed(a, b)) {
        edge_pair[ntrans] = a * STATES + b;
        edge_rev[ntrans] = a >= 7 || b >= 7;
        edge_kind[ntrans] = boundary_kind(a, b);
        pair_to_edge[a][b] = ntrans;
        trans_a[ntrans] = a, trans_b[ntrans++] = b;
      }
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
  indices[1] = OFF_E4 + state * K4 + site->k4_emit[state >= 7];
}

static int transition_indices(const Site *site, int previous, int current, int *indices) {
  int pair = previous * STATES + current, count = 2;
  int reverse = previous >= 7 || current >= 7, kind = boundary_kind(previous, current);
  indices[0] = OFF_TNUC + pair * KNUC + site->nuc;
  indices[1] = OFF_T4 + pair * K4 + site->k4_trans[reverse];
  if (kind >= 0)
    for (int offset = 0; offset < FLANKS; offset++)
      indices[count++] = OFF_SPLICE + (kind * FLANKS + offset) * K2 + site->flank[reverse][offset];
  return count;
}

static float derivative_value(const Site *site, int scale, int reverse) {
  float value = site->derivative[scale] / (float)INT16_MAX;
  return reverse ? -value : value;
}

static double feature_score(const double *weights, const int *indices, int count) {
  double score = 0;
  for (int index = 0; index < count; index++) score += weights[indices[index]];
  return score;
}

static void add_features(double *gradient, const int *indices, int count, double value) {
  for (int index = 0; index < count; index++) gradient[indices[index]] += value;
}

static double transition(const double *weights, const Contig *c, const Site *site, int position,
                         int previous, int current) {
  int indices[TRANSITION_FEATURES], reverse = previous >= 7 || current >= 7;
  int count = transition_indices(site, previous, current, indices);
  int pair = previous * STATES + current, kind = boundary_kind(previous, current);
  double score = feature_score(weights, indices, count);
  for (int scale = 0; scale < WAVE_COUNT; scale++) {
    score += weights[OFF_TP + pair * WAVE_COUNT + scale] * site->cwt_power[scale];
    score += weights[OFF_TD + pair * WAVE_COUNT + scale] * derivative_value(site, scale, reverse);
  }
  if (kind >= 0) {
    int genome_position = c->offset + position;
    for (int offset = 0; offset < FLANKS; offset++) {
      float power[WAVE_COUNT];
      double derivative[WAVE_COUNT];
      boundary_cwt_at(c, genome_position, reverse, offset, power, derivative);
      int base = (kind * FLANKS + offset) * WAVE_COUNT;
      for (int scale = 0; scale < WAVE_COUNT; scale++) {
        score += weights[OFF_BP + base + scale] * power[scale];
        score += weights[OFF_BD + base + scale] * derivative[scale];
      }
    }
  }
  return score;
}

static void add_transition(double *gradient, const Contig *c, const Site *site, int position,
                           int previous, int current, double value) {
  int indices[TRANSITION_FEATURES], reverse = previous >= 7 || current >= 7;
  int count = transition_indices(site, previous, current, indices);
  add_features(gradient, indices, count, value);
  int pair = previous * STATES + current, kind = boundary_kind(previous, current);
  for (int scale = 0; scale < WAVE_COUNT; scale++) {
    gradient[OFF_TP + pair * WAVE_COUNT + scale] += value * site->cwt_power[scale];
    gradient[OFF_TD + pair * WAVE_COUNT + scale] += value * derivative_value(site, scale, reverse);
  }
  if (kind >= 0) {
    int genome_position = c->offset + position;
    for (int offset = 0; offset < FLANKS; offset++) {
      float power[WAVE_COUNT];
      double derivative[WAVE_COUNT];
      boundary_cwt_at(c, genome_position, reverse, offset, power, derivative);
      int base = (kind * FLANKS + offset) * WAVE_COUNT;
      for (int scale = 0; scale < WAVE_COUNT; scale++) {
        gradient[OFF_BP + base + scale] += value * power[scale];
        gradient[OFF_BD + base + scale] += value * derivative[scale];
      }
    }
  }
}

typedef struct {
  Site *sites;
  double *exp_psi, *alpha, *inv_scale;
  int capacity;
} Work;

static double log_add(double left, double right) {
  if (left == -INFINITY) return right;
  if (right == -INFINITY) return left;
  double maximum = fmax(left, right);
  return maximum + log1p(exp(-fabs(left - right)));
}

static Site site_scores(const double *w, const Contig *c, int position,
                        double *emissions, double *psi) {
  Site site = site_at(c, position);
  float deriv_pos[WAVE_COUNT];
  for (int scale = 0; scale < WAVE_COUNT; scale++)
    deriv_pos[scale] = site.derivative[scale] / (float)INT16_MAX;

  for (int state = 0; state < STATES; state++) {
    int reverse = state >= 7;
    double score = w[OFF_ENUC + state * KNUC + site.nuc] +
                   w[OFF_E4 + state * K4 + site.k4_emit[reverse]];
    const double *wp = w + OFF_EP + state * WAVE_COUNT;
    const double *wd = w + OFF_ED + state * WAVE_COUNT;
    for (int scale = 0; scale < WAVE_COUNT; scale++) {
      score += wp[scale] * site.cwt_power[scale];
      score += wd[scale] * (reverse ? -deriv_pos[scale] : deriv_pos[scale]);
    }
    if (!isfinite(score))
      anno_fail("Non-finite CRF emission at position %d", position + 1);
    emissions[state] = score;
  }

  int genome_position = c->offset + position;
  double b_splice[2][4], b_cwt[2][4];
  for (int rev = 0; rev < 2; rev++) {
    float power[FLANKS][WAVE_COUNT];
    double derivative[FLANKS][WAVE_COUNT];
    for (int offset = 0; offset < FLANKS; offset++)
      boundary_cwt_at(c, genome_position, rev, offset, power[offset], derivative[offset]);
    const uint8_t *fl = site.flank[rev];
    for (int kind = 0; kind < 4; kind++) {
      double sp = 0, cw = 0;
      for (int offset = 0; offset < FLANKS; offset++)
        sp += w[OFF_SPLICE + (kind * FLANKS + offset) * K2 + fl[offset]];
      for (int offset = 0; offset < FLANKS; offset++) {
        int base = (kind * FLANKS + offset) * WAVE_COUNT;
        const double *wbp = w + OFF_BP + base;
        const double *wbd = w + OFF_BD + base;
        for (int scale = 0; scale < WAVE_COUNT; scale++) {
          cw += wbp[scale] * power[offset][scale];
          cw += wbd[scale] * derivative[offset][scale];
        }
      }
      b_splice[rev][kind] = sp;
      b_cwt[rev][kind] = cw;
    }
  }

  for (int edge = 0; edge < ntrans; edge++) {
    int pair = edge_pair[edge], reverse = edge_rev[edge], kind = edge_kind[edge];
    double score = w[OFF_TNUC + pair * KNUC + site.nuc] +
                   w[OFF_T4 + pair * K4 + site.k4_trans[reverse]];
    if (kind >= 0) score += b_splice[reverse][kind];
    const double *wp = w + OFF_TP + pair * WAVE_COUNT;
    const double *wd = w + OFF_TD + pair * WAVE_COUNT;
    for (int scale = 0; scale < WAVE_COUNT; scale++) {
      score += wp[scale] * site.cwt_power[scale];
      score += wd[scale] * (reverse ? -deriv_pos[scale] : deriv_pos[scale]);
    }
    if (kind >= 0) score += b_cwt[reverse][kind];
    double p = score + emissions[trans_b[edge]];
    if (!isfinite(p))
      anno_fail("Non-finite CRF transition at position %d", position + 1);
    psi[edge] = p;
  }
  return site;
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

static void work_ensure(Work *k, int length) {
  if (length <= k->capacity) return;
  free(k->sites), free(k->exp_psi), free(k->alpha), free(k->inv_scale);
  k->sites = anno_alloc((size_t)length, sizeof(Site));
  k->exp_psi = anno_alloc((size_t)length, (size_t)ntrans * sizeof(double));
  k->alpha = anno_alloc((size_t)length, HIDDEN_STATES * sizeof(double));
  k->inv_scale = anno_alloc((size_t)length, sizeof(double));
  k->capacity = length;
}

static void work_free(Work *k) {
  free(k->sites), free(k->exp_psi), free(k->alpha), free(k->inv_scale);
}

/* Accumulates the exact observed-minus-expected gradient for one training window. */
static double contig_gradient(const double *w, const Contig *c, double mass,
                              Work *k, double *gradient) {
  if (mass == -INFINITY) return NAN;

  work_ensure(k, c->length);
  double logz = 0, gold_score = 0;
  for (int position = 0; position < c->length; position++) {
    double emissions[STATES], psi[STATES * STATES];
    Site site = site_scores(w, c, position, emissions, psi);
    k->sites[position] = site;
    int gold_state = c->label[c->offset + position];
    gold_score += emissions[gold_state];

    double *A = k->alpha + (size_t)position * HIDDEN_STATES;
    if (!position) {
      double max_e = -INFINITY;
      for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
        if (initial_allowed(site.nuc, hidden, c->offset == 0)) {
          double e = emissions[hidden_state[hidden]];
          if (e > max_e) max_e = e;
        }
      if (!isfinite(max_e)) anno_fail("Invalid CRF partition at position %d", position + 1);
      double sum = 0;
      for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) {
        double v = initial_allowed(site.nuc, hidden, c->offset == 0)
          ? exp(emissions[hidden_state[hidden]] - max_e) : 0.0;
        if (c->length == 1 && c->offset + c->length == c->genome_length && !terminal_allowed(hidden))
          v = 0.0;
        A[hidden] = v;
        sum += v;
      }
      if (!(sum > 0) || !isfinite(sum)) anno_fail("Invalid CRF partition at position %d", position + 1);
      double inv = 1.0 / sum;
      for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) A[hidden] *= inv;
      k->inv_scale[0] = inv;
      logz += max_e + log(sum);
    } else {
      double max_psi = psi[0];
      for (int edge = 1; edge < ntrans; edge++)
        if (psi[edge] > max_psi) max_psi = psi[edge];

      int prev_gold = c->label[c->offset + position - 1];
      int gold_edge = pair_to_edge[prev_gold][gold_state];
      if (gold_edge >= 0)
        gold_score += psi[gold_edge] - emissions[gold_state];
      else
        gold_score += transition(w, c, &site, position, prev_gold, gold_state);

      double *ep = k->exp_psi + (size_t)position * ntrans;
      for (int edge = 0; edge < ntrans; edge++) ep[edge] = exp(psi[edge] - max_psi);

      const double *prev_A = A - HIDDEN_STATES;
      for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) A[hidden] = 0.0;
      const HiddenEdge *edges = hidden_edges[site.nuc];
      int edge_count = hidden_edge_count[site.nuc];
      for (int edge = 0; edge < edge_count; edge++) {
        const HiddenEdge *step = edges + edge;
        A[step->current] += prev_A[step->previous] * ep[step->transition];
      }
      if (position == c->length - 1 && c->offset + c->length == c->genome_length)
        for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
          if (!terminal_allowed(hidden)) A[hidden] = 0.0;

      double sum = 0;
      for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) sum += A[hidden];
      if (!(sum > 0) || !isfinite(sum))
        anno_fail("Invalid CRF partition at position %d", position + 1);
      double inv = 1.0 / sum;
      for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) A[hidden] *= inv;
      k->inv_scale[position] = inv;
      logz += max_psi + log(sum);
    }
  }

  double likelihood = mass - logz + gold_score;
  double B[HIDDEN_STATES];
  int reaches_contig_end = c->offset + c->length == c->genome_length;
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
    B[hidden] = !reaches_contig_end || terminal_allowed(hidden) ? 1.0 : 0.0;

  for (int position = c->length - 1; position >= 0; position--) {
    const Site *site = k->sites + position;
    const double *A = k->alpha + (size_t)position * HIDDEN_STATES;
    int state = c->label[c->offset + position];

    double state_marginals[STATES] = {0};
    for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
      state_marginals[hidden_state[hidden]] += A[hidden] * B[hidden];
    state_marginals[state] -= 1.0;

    float deriv_pos[WAVE_COUNT];
    for (int scale = 0; scale < WAVE_COUNT; scale++)
      deriv_pos[scale] = site->derivative[scale] / (float)INT16_MAX;

    for (int st = 0; st < STATES; st++) {
      double val = -state_marginals[st];
      int reverse = st >= 7;
      gradient[OFF_ENUC + st * KNUC + site->nuc] += val;
      gradient[OFF_E4 + st * K4 + site->k4_emit[reverse]] += val;
      double *gp = gradient + OFF_EP + st * WAVE_COUNT;
      double *gd = gradient + OFF_ED + st * WAVE_COUNT;
      for (int scale = 0; scale < WAVE_COUNT; scale++) {
        gp[scale] += val * site->cwt_power[scale];
        gd[scale] += val * (reverse ? -deriv_pos[scale] : deriv_pos[scale]);
      }
    }

    if (position) {
      const double *prev_A = A - HIDDEN_STATES;
      const double *ep = k->exp_psi + (size_t)position * ntrans;
      double inv_c = k->inv_scale[position];
      int previous_state = c->label[c->offset + position - 1];
      int gold_edge = pair_to_edge[previous_state][state];

      double edge_marginals[STATES * STATES] = {0};
      double B_prev[HIDDEN_STATES] = {0};
      const HiddenEdge *edges = hidden_edges[site->nuc];
      int edge_count = hidden_edge_count[site->nuc];
      for (int edge = 0; edge < edge_count; edge++) {
        const HiddenEdge *step = edges + edge;
        double w_bq = ep[step->transition] * B[step->current];
        B_prev[step->previous] += w_bq;
        edge_marginals[step->transition] += prev_A[step->previous] * w_bq;
      }
      for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) B[hidden] = B_prev[hidden] * inv_c;

      double b_coeff[4][2] = {{0}};
      int has_boundary = 0;
      for (int edge = 0; edge < ntrans; edge++) edge_marginals[edge] *= inv_c;
      if (gold_edge >= 0) {
        edge_marginals[gold_edge] -= 1.0;
      } else {
        add_transition(gradient, c, site, position, previous_state, state, 1);
      }

      for (int edge = 0; edge < ntrans; edge++) {
        double m = edge_marginals[edge];
        if (m != 0) {
          double val = -m;
          int pair = edge_pair[edge], reverse = edge_rev[edge], kind = edge_kind[edge];
          gradient[OFF_TNUC + pair * KNUC + site->nuc] += val;
          gradient[OFF_T4 + pair * K4 + site->k4_trans[reverse]] += val;
          if (kind >= 0) {
            b_coeff[kind][reverse] += val;
            has_boundary = 1;
          }
          double *gp = gradient + OFF_TP + pair * WAVE_COUNT;
          double *gd = gradient + OFF_TD + pair * WAVE_COUNT;
          for (int scale = 0; scale < WAVE_COUNT; scale++) {
            gp[scale] += val * site->cwt_power[scale];
            gd[scale] += val * (reverse ? -deriv_pos[scale] : deriv_pos[scale]);
          }
        }
      }
      if (has_boundary) {
        int genome_position = c->offset + position;
        for (int rev = 0; rev < 2; rev++) {
          if (b_coeff[0][rev] == 0 && b_coeff[1][rev] == 0 &&
              b_coeff[2][rev] == 0 && b_coeff[3][rev] == 0) continue;
          float power[FLANKS][WAVE_COUNT];
          double derivative[FLANKS][WAVE_COUNT];
          for (int offset = 0; offset < FLANKS; offset++)
            boundary_cwt_at(c, genome_position, rev, offset, power[offset], derivative[offset]);
          for (int kind = 0; kind < 4; kind++) {
            double val = b_coeff[kind][rev];
            if (val == 0) continue;
            const uint8_t *fl = site->flank[rev];
            for (int offset = 0; offset < FLANKS; offset++)
              gradient[OFF_SPLICE + (kind * FLANKS + offset) * K2 + fl[offset]] += val;
            for (int offset = 0; offset < FLANKS; offset++) {
              int base = (kind * FLANKS + offset) * WAVE_COUNT;
              double *gbp = gradient + OFF_BP + base;
              double *gbd = gradient + OFF_BD + base;
              for (int scale = 0; scale < WAVE_COUNT; scale++) {
                gbp[scale] += val * power[offset][scale];
                gbd[scale] += val * derivative[offset][scale];
              }
            }
          }
        }
      }
    }
  }
  return likelihood;
}

typedef struct {
  Contig *contig;
  int start, length;
  double mass;
} CrfWindow;

static int next_window_end(const Contig *contig, int start, int limit) {
  int end = limit - start > CRF_WINDOW_SIZE ? start + CRF_WINDOW_SIZE : limit;
  if (end == limit || !contig->label[end - 1] || !contig->label[end]) return end;
  int gene_start = end - 1, gene_end = end;
  while (gene_start > start && contig->label[gene_start - 1]) gene_start--;
  while (gene_end < limit && contig->label[gene_end]) gene_end++;
  return gene_start > start && end - gene_start <= gene_end - end ? gene_start : gene_end;
}

static size_t make_windows(Contig *contigs, int count, CrfWindow **windows) {
  size_t window_count = 0, capacity = 0;
  for (int i = 0; i < count; i++) {
    Contig *contig = contigs + i;
    int limit = contig->offset + contig->length;
    for (int start = contig->offset; start < limit;) {
      int end = next_window_end(contig, start, limit);
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
      Contig instance = *contig;
      instance.offset = start;
      instance.length = end - start;
      double mass = gold_grammar_mass(&instance, 0, instance.length);
      (*windows)[window_count++] = (CrfWindow){contig, start, end - start, mass};
      start = end;
    }
  }
  return window_count;
}

typedef struct {
  CrfWindow *windows;
  size_t count;
  int threads;
} CrfTrain;

typedef struct {
  CrfTrain *train;
  int thread_id;
  const double *w;
  double likelihood;
  size_t positions, skipped;
  double *gradient;
  Work work;
} CrfWorker;

static void *crf_worker(void *arg) {
  CrfWorker *worker = arg;
  CrfTrain *train = worker->train;
  worker->likelihood = 0;
  worker->positions = 0;
  worker->skipped = 0;
  memset(worker->gradient, 0, WEIGHTS * sizeof(*worker->gradient));
  for (size_t slot = (size_t)worker->thread_id; slot < train->count;
       slot += (size_t)train->threads) {
    const CrfWindow *window = train->windows + slot;
    Contig instance = *window->contig;
    instance.offset = window->start;
    instance.length = window->length;
    double score = contig_gradient(worker->w, &instance, window->mass,
                                   &worker->work, worker->gradient);
    if (isnan(score)) {
      worker->skipped++;
    } else {
      if (!isfinite(score)) anno_fail("Non-finite training likelihood");
      worker->likelihood += score;
      worker->positions += (size_t)instance.length;
    }
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
  double likelihood = 0;
  size_t positions = 0, skipped = 0;

  memset(gradient, 0, WEIGHTS * sizeof(*gradient));
  for (int t = 0; t < train->threads; t++) {
    objective->workers[t].w = w;
    if (train->threads > 1 &&
        pthread_create(objective->threads + t, NULL, crf_worker, objective->workers + t))
      anno_fail("Cannot spawn training thread");
  }
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
      int state = contigs[c].label[contigs[c].offset + i], indices[TRANSITION_FEATURES];
      emission_indices(&site, state, indices);
      for (int j = 0; j < EMISSION_FEATURES; j++) scale[indices[j]] += 1;
      for (int k = 0; k < WAVE_COUNT; k++) {
        double derivative = derivative_value(&site, k, state >= 7);
        scale[OFF_EP + state * WAVE_COUNT + k] += (double)site.cwt_power[k] * site.cwt_power[k];
        scale[OFF_ED + state * WAVE_COUNT + k] += derivative * derivative;
      }
      if (!i) continue;
      int previous = contigs[c].label[contigs[c].offset + i - 1], pair = previous * STATES + state;
      int n = transition_indices(&site, previous, state, indices);
      for (int j = 0; j < n; j++) scale[indices[j]] += 1;
      int reverse = previous >= 7 || state >= 7, kind = boundary_kind(previous, state);
      for (int k = 0; k < WAVE_COUNT; k++) {
        double derivative = derivative_value(&site, k, reverse);
        scale[OFF_TP + pair * WAVE_COUNT + k] += (double)site.cwt_power[k] * site.cwt_power[k];
        scale[OFF_TD + pair * WAVE_COUNT + k] += derivative * derivative;
      }
      if (kind >= 0) {
        int genome_position = contigs[c].offset + i;
        for (int offset = 0; offset < FLANKS; offset++) {
          float power[WAVE_COUNT];
          double derivative[WAVE_COUNT];
          boundary_cwt_at(contigs + c, genome_position, reverse, offset, power, derivative);
          int base = (kind * FLANKS + offset) * WAVE_COUNT;
          for (int k = 0; k < WAVE_COUNT; k++) {
            scale[OFF_BP + base + k] += (double)power[k] * power[k];
            scale[OFF_BD + base + k] += derivative[k] * derivative[k];
          }
        }
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
  fprintf(stderr, "Using %d worker threads for %zu training windows\n",
          objective.train.threads, objective.train.count);
  objective.workers = anno_alloc((size_t)objective.train.threads, sizeof(*objective.workers));
  objective.threads = anno_alloc((size_t)objective.train.threads, sizeof(*objective.threads));
  for (int t = 0; t < objective.train.threads; t++)
    objective.workers[t] = (CrfWorker){
      .train = &objective.train,
      .thread_id = t,
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
  for (int t = 0; t < objective.train.threads; t++) {
    work_free(&objective.workers[t].work);
    free(objective.workers[t].gradient);
  }
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
  double emissions[STATES], transitions[STATES * STATES];
  Site site = site_scores(weights, contig, position, emissions, transitions);
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) {
    current[hidden] = position == 0 && initial_allowed(site.nuc, hidden, 1) ? emissions[hidden_state[hidden]] : -INFINITY;
    if (back) back[hidden] = 0;
  }
  if (position) {
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
