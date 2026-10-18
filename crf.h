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
#ifndef CRF_DERIV_ORDERS
#define CRF_DERIV_ORDERS 1, 2, 3, 4
#endif
enum { DERIV_COUNT = CWT_COUNT(CRF_DERIV_ORDERS), DERIV_RADIUS = 2 };

/* Linear-chain CRF decoded jointly for both strands.
   0 intergenic; 1-3 CDS(+) codon position; 4-6 intron(+) carrying the codon position of the
   next exon base; 7-9 CDS(-) codon position; 10-12 intron(-). */
enum { STATES = 13, HIDDEN_STATES = 47 };
enum { K2 = 17, K3 = 65, KNUC = 5, RADIUS = 8, FLANKS = RADIUS };
enum { OFF_E3 = 0, OFF_SPLICE = OFF_E3 + STATES * K3,
       OFF_BP = OFF_SPLICE + 4 * FLANKS * K2,
       OFF_BD = OFF_BP + 4 * FLANKS * WAVE_COUNT,
       WEIGHTS = OFF_BD + DERIV_COUNT * 4 * FLANKS * WAVE_COUNT };

void crf_init(void);
void crf_train(Contig *contigs, int count, double *weights, int iterations, int threads);
void crf_decode(const double *weights, const Contig *contig, uint8_t *path);
void crf_save_weights_text(const char *path, const double *weights);

#endif

#ifdef CRF_IMPLEMENTATION

#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  uint8_t k3_emit[2];
  uint8_t nuc, flank[2][FLANKS];
} Site;

static const int hidden_offsets[STATES] = {0, 1, 4, 8, 10, 11, 14, 18, 20, 28, 33, 41, 46};
static const int hidden_counts[STATES] = {1, 3, 4, 2, 1, 3, 4, 2, 8, 5, 8, 5, 1};
static int trans_a[STATES * STATES], trans_b[STATES * STATES], ntrans;
static int edge_rev[STATES * STATES], edge_kind[STATES * STATES];
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
    .k3_emit = {kmer(c->seq, c->genome_length, position - 1, 3, 0),
                kmer(c->seq, c->genome_length, position - 1, 3, 1)},
    .nuc = base < 0 ? 4 : base
  };
  for (int offset = 0; offset < FLANKS; offset++) {
    int displacement = -RADIUS + 2 * offset;
    s.flank[0][offset] = kmer(c->seq, c->genome_length, position + displacement, 2, 0);
    s.flank[1][offset] = kmer(c->seq, c->genome_length, position - displacement - 2, 2, 1);
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

static const int boundary_cwt_offsets[FLANKS] = {-27, -26, -25, -24, -2, -1, 0, 1};
static const int derivative_orders[DERIV_COUNT] = {CRF_DERIV_ORDERS};
static double derivative_kernels[DERIV_COUNT][2 * DERIV_RADIUS + 1];

static void derivative_init(void) {
  for (int slot = 0; slot < DERIV_COUNT; slot++) {
    int order = derivative_orders[slot];
    if (order < 1 || order > 4) anno_fail("CWT derivative order must be between 1 and 4");
    for (int previous = 0; previous < slot; previous++)
      if (derivative_orders[previous] == order) anno_fail("Duplicate CWT derivative order: %d", order);
    double *kernel = derivative_kernels[slot];
    memset(kernel, 0, sizeof(derivative_kernels[slot]));
    if (order % 2) kernel[DERIV_RADIUS - 1] = -.5, kernel[DERIV_RADIUS + 1] = .5;
    else kernel[DERIV_RADIUS] = 1;
    for (int degree = 0; degree < order / 2; degree++) {
      double next[2 * DERIV_RADIUS + 1] = {0};
      for (int j = 0; j <= 2 * DERIV_RADIUS; j++) {
        if (j) next[j - 1] += kernel[j];
        next[j] -= 2 * kernel[j];
        if (j < 2 * DERIV_RADIUS) next[j + 1] += kernel[j];
      }
      memcpy(kernel, next, sizeof(next));
    }
  }
}

static int derivative_index(int slot, int kind, int offset, int scale) {
  return OFF_BD + ((slot * 4 + kind) * FLANKS + offset) * WAVE_COUNT + scale;
}

static void boundary_cwt_at(const Contig *c, int position, int reverse, int offset,
                            float *power, double derivative[DERIV_COUNT][WAVE_COUNT]) {
  int displacement = boundary_cwt_offsets[offset];
  int g = reverse ? position - displacement - 1 : position + displacement;
  if (g < 0 || g >= c->genome_length || base_index(c->seq[g]) < 0) {
    memset(power, 0, WAVE_COUNT * sizeof(*power));
    memset(derivative, 0, DERIV_COUNT * WAVE_COUNT * sizeof(double));
    return;
  }
  memcpy(power, c->cwt_power + (size_t)g * WAVE_COUNT, WAVE_COUNT * sizeof(*power));
  for (int slot = 0; slot < DERIV_COUNT; slot++) {
    int order = derivative_orders[slot];
    const double *kernel = derivative_kernels[slot];
    int valid = 1;
    for (int j = -DERIV_RADIUS; j <= DERIV_RADIUS; j++) {
      if (!kernel[j + DERIV_RADIUS]) continue;
      if (j < -g || j >= c->genome_length - g ||
          base_index(c->seq[g + j]) < 0) {
        valid = 0;
        break;
      }
    }
    for (int scale = 0; scale < WAVE_COUNT; scale++) {
      double value = 0;
      if (valid)
        for (int j = -DERIV_RADIUS; j <= DERIV_RADIUS; j++) {
          double coefficient = kernel[j + DERIV_RADIUS];
          if (!coefficient) continue;
          value += coefficient * c->cwt_power[(size_t)(g + j) * WAVE_COUNT + scale];
        }
      derivative[slot][scale] = reverse && order % 2 ? -value : value;
    }
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
  derivative_init();
  ntrans = 0;
  memset(pair_to_edge, -1, sizeof(pair_to_edge));
  for (int a = 0; a < STATES; a++)
    for (int b = 0; b < STATES; b++)
      if (allowed(a, b)) {
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

static int emission_index(const Site *site, int state) {
  return OFF_E3 + state * K3 + site->k3_emit[state >= 7];
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
  for (int state = 0; state < STATES; state++) {
    double score = w[emission_index(&site, state)];
    if (!isfinite(score))
      anno_fail("Non-finite CRF emission at position %d", position + 1);
    emissions[state] = score;
  }

  int genome_position = c->offset + position;
  double b_splice[2][4], b_cwt[2][4];
  for (int rev = 0; rev < 2; rev++) {
    float power[FLANKS][WAVE_COUNT];
    double derivative[FLANKS][DERIV_COUNT][WAVE_COUNT];
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
        for (int scale = 0; scale < WAVE_COUNT; scale++) {
          cw += wbp[scale] * power[offset][scale];
          for (int slot = 0; slot < DERIV_COUNT; slot++)
            cw += w[derivative_index(slot, kind, offset, scale)] * derivative[offset][slot][scale];
        }
      }
      b_splice[rev][kind] = sp;
      b_cwt[rev][kind] = cw;
    }
  }

  for (int edge = 0; edge < ntrans; edge++) {
    int reverse = edge_rev[edge], kind = edge_kind[edge];
    double score = 0;
    if (kind >= 0) {
      score += b_splice[reverse][kind];
      score += b_cwt[reverse][kind];
    }
    double p = score + emissions[trans_b[edge]];
    if (!isfinite(p))
      anno_fail("Non-finite CRF transition at position %d", position + 1);
    psi[edge] = p;
  }
  return site;
}

static double gold_grammar_mass(const Contig *contig) {
  double previous[HIDDEN_STATES], current[HIDDEN_STATES];
  int position = contig->offset;
  int base = base_index(contig->seq[position]);
  if (base < 0) base = 4;
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
    previous[hidden] = hidden_state[hidden] == contig->label[position] &&
      initial_allowed(base, hidden, contig->offset == 0) ? 0 : -INFINITY;
  for (int position = 1; position < contig->length; position++) {
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
    if (contig->offset + contig->length < contig->genome_length || terminal_allowed(hidden))
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
        if (initial_allowed(site.nuc, hidden, c->offset == 0) &&
            (c->length > 1 || c->offset + c->length < c->genome_length || terminal_allowed(hidden))) {
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
      int prev_gold = c->label[c->offset + position - 1];
      int gold_edge = pair_to_edge[prev_gold][gold_state];
      if (gold_edge < 0) anno_fail("Invalid transition in validated CRF training path");
      gold_score += psi[gold_edge] - emissions[gold_state];

      const double *prev_A = A - HIDDEN_STATES;
      const HiddenEdge *edges = hidden_edges[site.nuc];
      int edge_count = hidden_edge_count[site.nuc];
      double max_psi = -INFINITY;
      for (int edge = 0; edge < edge_count; edge++) {
        const HiddenEdge *step = edges + edge;
        if (prev_A[step->previous] > 0.0 &&
            (position < c->length - 1 || c->offset + c->length < c->genome_length ||
             terminal_allowed(step->current))) {
          double p = psi[step->transition];
          if (p > max_psi) max_psi = p;
        }
      }
      if (!isfinite(max_psi))
        anno_fail("Invalid CRF partition at position %d", position + 1);

      double *ep = k->exp_psi + (size_t)position * ntrans;
      for (int edge = 0; edge < ntrans; edge++) ep[edge] = exp(psi[edge] - max_psi);

      for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) A[hidden] = 0.0;
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
  const double *last_A = k->alpha + (size_t)(c->length - 1) * HIDDEN_STATES;
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
    B[hidden] = last_A[hidden];

  for (int position = c->length - 1; position >= 0; position--) {
    const Site *site = k->sites + position;
    const double *A = k->alpha + (size_t)position * HIDDEN_STATES;
    int state = c->label[c->offset + position];

    double state_marginals[STATES] = {0};
    for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
      state_marginals[hidden_state[hidden]] += B[hidden];
    state_marginals[state] -= 1.0;

    for (int st = 0; st < STATES; st++) {
      double val = -state_marginals[st];
      gradient[emission_index(site, st)] += val;
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
        double a_cur = A[step->current], b_cur = B[step->current];
        if (a_cur > 0.0 && b_cur > 0.0) {
          double edge_prob = ((prev_A[step->previous] * ep[step->transition] * inv_c) / a_cur) * b_cur;
          B_prev[step->previous] += edge_prob;
          edge_marginals[step->transition] += edge_prob;
        }
      }
      double sum_b = 0;
      for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) sum_b += B_prev[hidden];
      double inv_b = sum_b > 0.0 ? 1.0 / sum_b : 0.0;
      for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) B[hidden] = B_prev[hidden] * inv_b;

      double b_coeff[4][2] = {{0}};
      int has_boundary = 0;
      edge_marginals[gold_edge] -= 1.0;

      for (int edge = 0; edge < ntrans; edge++) {
        double m = edge_marginals[edge];
        if (m != 0) {
          double val = -m;
          int reverse = edge_rev[edge], kind = edge_kind[edge];
          if (kind >= 0) {
            b_coeff[kind][reverse] += val;
            has_boundary = 1;
          }
        }
      }
      if (has_boundary) {
        int genome_position = c->offset + position;
        for (int rev = 0; rev < 2; rev++) {
          if (b_coeff[0][rev] == 0 && b_coeff[1][rev] == 0 &&
              b_coeff[2][rev] == 0 && b_coeff[3][rev] == 0) continue;
          float power[FLANKS][WAVE_COUNT];
          double derivative[FLANKS][DERIV_COUNT][WAVE_COUNT];
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
              for (int scale = 0; scale < WAVE_COUNT; scale++) {
                gbp[scale] += val * power[offset][scale];
                for (int slot = 0; slot < DERIV_COUNT; slot++)
                  gradient[derivative_index(slot, kind, offset, scale)] +=
                    val * derivative[offset][slot][scale];
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
      double mass = gold_grammar_mass(&instance);
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
      int state = contigs[c].label[contigs[c].offset + i];
      scale[emission_index(&site, state)] += 1;
      if (!i) continue;
      int previous = contigs[c].label[contigs[c].offset + i - 1];
      int reverse = previous >= 7 || state >= 7, kind = boundary_kind(previous, state);
      if (kind >= 0) {
        for (int offset = 0; offset < FLANKS; offset++)
          scale[OFF_SPLICE + (kind * FLANKS + offset) * K2 + site.flank[reverse][offset]] += 1;
        int genome_position = contigs[c].offset + i;
        for (int offset = 0; offset < FLANKS; offset++) {
          float power[WAVE_COUNT];
          double derivative[DERIV_COUNT][WAVE_COUNT];
          boundary_cwt_at(contigs + c, genome_position, reverse, offset, power, derivative);
          int base = (kind * FLANKS + offset) * WAVE_COUNT;
          for (int k = 0; k < WAVE_COUNT; k++) {
            scale[OFF_BP + base + k] += (double)power[k] * power[k];
            for (int slot = 0; slot < DERIV_COUNT; slot++)
              scale[derivative_index(slot, kind, offset, k)] +=
                derivative[slot][k] * derivative[slot][k];
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
    back[hidden] = 0;
  }
  if (position) {
    const HiddenEdge *edges = hidden_edges[site.nuc];
    for (int edge = 0; edge < hidden_edge_count[site.nuc]; edge++) {
      const HiddenEdge *step = edges + edge;
      double score = previous[step->previous] + transitions[step->transition];
      if (score > current[step->current]) {
        current[step->current] = score;
        back[step->current] = step->previous;
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

static void format_kmer(int code, int size, char *out) {
  static const char bases[4] = {'A', 'C', 'G', 'T'};
  if (code >= (1 << (2 * size))) {
    out[0] = 'N';
    out[1] = '\0';
    return;
  }
  for (int i = size - 1; i >= 0; i--) {
    out[i] = bases[code & 3];
    code >>= 2;
  }
  out[size] = '\0';
}

void crf_save_weights_text(const char *path, const double *w) {
  static const char *const state_names[STATES] = {
    "intergenic",
    "CDS(+):0", "CDS(+):1", "CDS(+):2",
    "intron(+):0", "intron(+):1", "intron(+):2",
    "CDS(-):0", "CDS(-):1", "CDS(-):2",
    "intron(-):0", "intron(-):1", "intron(-):2"
  };
  static const char *const kind_names[4] = {
    "5ss", "3ss", "start_codon", "stop_codon"
  };

  FILE *file = fopen(path, "w");
  if (!file) anno_fail("Cannot write weights text: %s", path);
  fprintf(file, "# index\tgroup\ttarget\tfeature\tweight\n");

  char kbuf[8];
  for (int s = 0; s < STATES; s++) {
    for (int k = 0; k < K3; k++) {
      int idx = OFF_E3 + s * K3 + k;
      format_kmer(k, 3, kbuf);
      fprintf(file, "%d\temit_k3\t%s\tk3=%s\t%.9g\n",
              idx, state_names[s], kbuf, (float)w[idx]);
    }
  }

  for (int kind = 0; kind < 4; kind++) {
    for (int offset = 0; offset < FLANKS; offset++) {
      int disp = -RADIUS + 2 * offset;
      for (int k = 0; k < K2; k++) {
        int idx = OFF_SPLICE + (kind * FLANKS + offset) * K2 + k;
        format_kmer(k, 2, kbuf);
        fprintf(file, "%d\tboundary_k2\t%s\toffset=%+d,k2=%s\t%.9g\n",
                idx, kind_names[kind], disp, kbuf, (float)w[idx]);
      }
    }
  }
  for (int kind = 0; kind < 4; kind++) {
    for (int offset = 0; offset < FLANKS; offset++) {
      int disp = boundary_cwt_offsets[offset];
      int base = (kind * FLANKS + offset) * WAVE_COUNT;
      for (int sc = 0; sc < WAVE_COUNT; sc++) {
        int idx = OFF_BP + base + sc;
        fprintf(file, "%d\tboundary_cwt_power\t%s\toffset=%+d,scale=%d\t%.9g\n",
                idx, kind_names[kind], disp, wave_sizes[sc], (float)w[idx]);
      }
    }
  }
  for (int slot = 0; slot < DERIV_COUNT; slot++) {
    for (int kind = 0; kind < 4; kind++) {
      for (int offset = 0; offset < FLANKS; offset++) {
        int disp = boundary_cwt_offsets[offset];
        for (int sc = 0; sc < WAVE_COUNT; sc++) {
          int idx = derivative_index(slot, kind, offset, sc);
          fprintf(file, "%d\tboundary_cwt_deriv%d\t%s\toffset=%+d,scale=%d\t%.9g\n",
                  idx, derivative_orders[slot], kind_names[kind], disp, wave_sizes[sc], (float)w[idx]);
        }
      }
    }
  }

  if (fclose(file)) anno_fail("Cannot finish writing weights text: %s", path);
}

#endif
