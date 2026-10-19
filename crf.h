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
#define CRF_DERIV_ORDERS 1, 2, 3
#endif
#ifndef CRF_EXON_LEN
#define CRF_EXON_LEN 512
#endif
#ifndef CRF_INTRON_LEN
#define CRF_INTRON_LEN 512
#endif
enum { DERIV_COUNT = CWT_COUNT(CRF_DERIV_ORDERS), DERIV_RADIUS = 2,
       INTRON_CHANNELS = 18, SM_BLK_SIZE = 256 };

/* Semi-Markov CRF decoded jointly for both strands.
   0 intergenic (1st-order Markov); 1-3 CDS(+) codon position (exact stop-bounded smCRF);
   4-6 intron(+) carrying the codon position of the next exon base (Method B smCRF);
   7-9 CDS(-) codon position; 10-12 intron(-). */
enum { STATES = 13, HIDDEN_STATES = 47 };
static const int boundary_kmer_offsets[] = {-8, -6, -4, -2, 0, 2, 4, 6};
static const int boundary_cwt_offsets[] = {-27, -26, -25, -24, -2, -1, 0, 1, 2};
enum {
  K2 = 17, K3 = 65, KNUC = 5,
  BOUNDARY_KMER_COUNT = sizeof(boundary_kmer_offsets) / sizeof(boundary_kmer_offsets[0]),
  BOUNDARY_CWT_COUNT = sizeof(boundary_cwt_offsets) / sizeof(boundary_cwt_offsets[0])
};
enum { OFF_E3 = 0, OFF_SPLICE_START_STOP = OFF_E3 + STATES * K3,
       OFF_BP = OFF_SPLICE_START_STOP + 4 * BOUNDARY_KMER_COUNT * K2,
       OFF_BD = OFF_BP + BOUNDARY_CWT_COUNT * WAVE_COUNT,
       OFF_EXON_LEN = OFF_BD + DERIV_COUNT * BOUNDARY_CWT_COUNT * WAVE_COUNT,
       OFF_EXON_TAIL = OFF_EXON_LEN + CRF_EXON_LEN,
       OFF_INTRON_LEN = OFF_EXON_TAIL + 1,
       OFF_INTRON_TAIL = OFF_INTRON_LEN + CRF_INTRON_LEN,
       WEIGHTS = OFF_INTRON_TAIL + 1 };

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
  uint8_t nuc, boundary_k2[2][BOUNDARY_KMER_COUNT];
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
  for (int offset = 0; offset < BOUNDARY_KMER_COUNT; offset++) {
    int displacement = boundary_kmer_offsets[offset];
    s.boundary_k2[0][offset] = kmer(c->seq, c->genome_length, position + displacement, 2, 0);
    s.boundary_k2[1][offset] = kmer(c->seq, c->genome_length, position - displacement - 2, 2, 1);
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

static const int derivative_orders[DERIV_COUNT] = {CRF_DERIV_ORDERS};
static double derivative_kernels[DERIV_COUNT][2 * DERIV_RADIUS + 1];

static void derivative_init(void) {
  for (int slot = 0; slot < DERIV_COUNT; slot++) {
    int order = derivative_orders[slot];
    if (order < 1 || order > 3) anno_fail("CWT derivative order must be between 1 and 3");
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

static int boundary_kmer_index(const Site *site, int reverse, int kind, int offset) {
  return OFF_SPLICE_START_STOP + (kind * BOUNDARY_KMER_COUNT + offset) * K2 + site->boundary_k2[reverse][offset];
}

static int derivative_index(int slot, int offset, int scale) {
  return OFF_BD + (slot * BOUNDARY_CWT_COUNT + offset) * WAVE_COUNT + scale;
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

static inline int is_cds_state(int s) { return (s >= 1 && s <= 3) || (s >= 7 && s <= 9); }
static inline int next_cds_state(int s) { return s <= 3 ? (1 + s % 3) : (7 + (s + 1) % 3); }
static inline int cds_state_at(int pos, int ch) {
  if (ch < 3) return 1 + (pos - ch + 3) % 3;
  return 7 + ((ch - 3) - pos % 3 + 3) % 3;
}
static inline int exon_start_channel(int t, int s0) {
  if (s0 >= 1 && s0 <= 3) return (t - (s0 - 1) + 3) % 3;
  if (s0 >= 7 && s0 <= 9) return 3 + (t + (s0 - 7)) % 3;
  return -1;
}
static inline int intron_h_to_ch(int h) {
  if (h >= 10 && h <= 17) return h - 10;
  if (h == 33) return 8;
  if (h >= 34 && h <= 36) return 9;
  if (h == 37) return 10;
  if (h >= 38 && h <= 40) return 11;
  if (h >= 41 && h <= 46) return 12 + (h - 41);
  return -1;
}
static const int intron_ch_rep[INTRON_CHANNELS] = {
  10, 11, 12, 13, 14, 15, 16, 17,
  33, 34, 37, 38, 41, 42, 43, 44, 45, 46
};
static const int intron_ch_to_iq[INTRON_CHANNELS] = {
  0, 1, 1, 1, 2, 2, 2, 2,
  3, 3, 3, 3, 4, 4, 4, 4, 4, 5
};
static uint8_t intron_can_exit[KNUC][INTRON_CHANNELS];

static inline double phi_exon(const double *w, int d) {
  if (d <= CRF_EXON_LEN) return w[OFF_EXON_LEN + d - 1];
  return w[OFF_EXON_LEN + CRF_EXON_LEN - 1] + w[OFF_EXON_TAIL] * (d - CRF_EXON_LEN);
}

static inline double phi_intron(const double *w, int d) {
  if (d <= CRF_INTRON_LEN) return w[OFF_INTRON_LEN + d - 1];
  return w[OFF_INTRON_LEN + CRF_INTRON_LEN - 1] + w[OFF_INTRON_TAIL] * (d - CRF_INTRON_LEN);
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
  memset(intron_can_exit, 0, sizeof(intron_can_exit));
  for (int base = 0; base < KNUC; base++)
    for (int previous = 0; previous < HIDDEN_STATES; previous++)
      for (int edge = 0; edge < ntrans; edge++) {
        if (trans_a[edge] != hidden_state[previous]) continue;
        int current = grammar_next(previous, trans_b[edge], base);
        if (current < 0) continue;
        if (hidden_edge_count[base] >= HIDDEN_STATES * 4) anno_fail("Too many grammar transitions");
        hidden_edges[base][hidden_edge_count[base]++] = (HiddenEdge){(uint8_t)previous, (uint8_t)current, (uint16_t)edge};
        initial_states[base][current] = 1;
        if (edge == 0 || edge_kind[edge] >= 0) {
          int ic = intron_h_to_ch(previous);
          if (ic >= 0) intron_can_exit[base][ic] = 1;
        }
      }
}

static int emission_index(const Site *site, int state) {
  return OFF_E3 + state * K3 + site->k3_emit[state >= 7];
}

typedef struct {
  Site *sites;
  double *psi;
  double *e_exon, *e_intron, *S_exon, *S_intron;
  int8_t *canon_h;
  int32_t *min_t;
  double *log_a_start, *log_a_end, *log_exon_merged;
  double *X_exon, *X_intron, *exp_X_exon, *exp_X_intron;
  double *Z_exon, *Z_intron;
  double *log_tail, *t_enter;
  double *P_end, *P_start, *P_tail;
  double *K_exon_rev, *G_exon_rev;
  int capacity, exon_k_cap;
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
    float power[BOUNDARY_CWT_COUNT][WAVE_COUNT];
    double derivative[BOUNDARY_CWT_COUNT][DERIV_COUNT][WAVE_COUNT];
    for (int offset = 0; offset < BOUNDARY_CWT_COUNT; offset++)
      boundary_cwt_at(c, genome_position, rev, offset, power[offset], derivative[offset]);
    for (int kind = 0; kind < 4; kind++) {
      double sp = 0, cw = 0;
      for (int offset = 0; offset < BOUNDARY_KMER_COUNT; offset++)
        sp += w[boundary_kmer_index(&site, rev, kind, offset)];
      for (int offset = 0; kind == 1 && offset < BOUNDARY_CWT_COUNT; offset++) {
        int base = offset * WAVE_COUNT;
        const double *wbp = w + OFF_BP + base;
        for (int scale = 0; scale < WAVE_COUNT; scale++) {
          cw += wbp[scale] * power[offset][scale];
          for (int slot = 0; slot < DERIV_COUNT; slot++)
            cw += w[derivative_index(slot, offset, scale)] * derivative[offset][slot][scale];
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
  size_t L = (size_t)length;
  size_t e_stride = L + 8, i_stride = L + 8;
  free(k->sites), free(k->psi);
  free(k->e_exon), free(k->e_intron), free(k->S_exon), free(k->S_intron);
  free(k->canon_h), free(k->min_t);
  free(k->log_a_start), free(k->log_a_end), free(k->log_exon_merged);
  free(k->X_exon), free(k->X_intron), free(k->exp_X_exon), free(k->exp_X_intron);
  free(k->Z_exon), free(k->Z_intron);
  free(k->log_tail), free(k->t_enter);
  free(k->P_end), free(k->P_start), free(k->P_tail);

  k->sites = anno_alloc(L, sizeof(Site));
  k->psi = anno_alloc(L, (size_t)ntrans * sizeof(double));
  k->e_exon = anno_alloc(L, 6 * sizeof(double));
  k->e_intron = anno_alloc(L, 6 * sizeof(double));
  k->S_exon = anno_alloc(L, 6 * sizeof(double));
  k->S_intron = anno_alloc(L, 6 * sizeof(double));
  k->canon_h = anno_alloc(L, 6 * sizeof(int8_t));
  k->min_t = anno_alloc(L, 6 * sizeof(int32_t));
  k->log_a_start = anno_alloc(L, HIDDEN_STATES * sizeof(double));
  k->log_a_end = anno_alloc(L, HIDDEN_STATES * sizeof(double));
  k->log_exon_merged = anno_alloc(L, 6 * sizeof(double));
  k->X_exon = anno_alloc(L, 6 * sizeof(double));
  k->X_intron = anno_alloc(L, INTRON_CHANNELS * sizeof(double));
  k->exp_X_exon = anno_alloc(6, e_stride * sizeof(double));
  k->exp_X_intron = anno_alloc(INTRON_CHANNELS, i_stride * sizeof(double));
  k->Z_exon = anno_alloc(6, e_stride * sizeof(double));
  k->Z_intron = anno_alloc(INTRON_CHANNELS, i_stride * sizeof(double));
  k->log_tail = anno_alloc(L, INTRON_CHANNELS * sizeof(double));
  k->t_enter = anno_alloc(L, INTRON_CHANNELS * sizeof(double));
  k->P_end = anno_alloc(L, HIDDEN_STATES * sizeof(double));
  k->P_start = anno_alloc(L, HIDDEN_STATES * sizeof(double));
  k->P_tail = anno_alloc(L, INTRON_CHANNELS * sizeof(double));
  k->capacity = length;
}

static void work_ensure_exon_k(Work *k, int max_d) {
  if (max_d <= k->exon_k_cap) return;
  free(k->K_exon_rev), free(k->G_exon_rev);
  k->K_exon_rev = anno_alloc((size_t)max_d, sizeof(double));
  k->G_exon_rev = anno_alloc((size_t)max_d, sizeof(double));
  k->exon_k_cap = max_d;
}

static void work_free(Work *k) {
  free(k->sites), free(k->psi);
  free(k->e_exon), free(k->e_intron), free(k->S_exon), free(k->S_intron);
  free(k->canon_h), free(k->min_t);
  free(k->log_a_start), free(k->log_a_end), free(k->log_exon_merged);
  free(k->X_exon), free(k->X_intron), free(k->exp_X_exon), free(k->exp_X_intron);
  free(k->Z_exon), free(k->Z_intron);
  free(k->log_tail), free(k->t_enter);
  free(k->P_end), free(k->P_start), free(k->P_tail);
  free(k->K_exon_rev), free(k->G_exon_rev);
}

/* Accumulates the exact observed-minus-expected smCRF gradient for one training window. */
static double contig_gradient(const double *w, const Contig *c, double mass,
                              Work *k, double *gradient) {
  if (mass == -INFINITY) return NAN;

  int length = c->length;
  int full_start = (c->offset == 0);
  int full_end = (c->offset + c->length == c->genome_length);
  work_ensure(k, length);
  size_t e_stride = (size_t)length + 8, i_stride = (size_t)length + 8;

  Site *sites = k->sites;
  double *psi = k->psi;
  double *e_exon = k->e_exon, *e_intron = k->e_intron;
  double *S_exon = k->S_exon, *S_intron = k->S_intron;
  int8_t *canon_h = k->canon_h;
  int32_t *min_t = k->min_t;
  double *log_a_start = k->log_a_start, *log_a_end = k->log_a_end;
  double *log_exon_merged = k->log_exon_merged;
  double *X_exon = k->X_exon, *X_intron = k->X_intron;
  double *exp_X_exon = k->exp_X_exon, *exp_X_intron = k->exp_X_intron;
  double *Z_exon = k->Z_exon, *Z_intron = k->Z_intron;
  double *log_tail = k->log_tail, *t_enter = k->t_enter;
  double *P_end = k->P_end, *P_start = k->P_start, *P_tail = k->P_tail;

  double gold_score = 0;
  for (int i = 0; i < length; i++) {
    double emissions[STATES], tr[STATES * STATES];
    sites[i] = site_scores(w, c, i, emissions, tr);
    memcpy(psi + (size_t)i * ntrans, tr, (size_t)ntrans * sizeof(double));
    for (int ch = 0; ch < 6; ch++) {
      int s = cds_state_at(i, ch);
      e_exon[i * 6 + ch] = emissions[s];
      S_exon[i * 6 + ch] = (i > 0 ? S_exon[(i - 1) * 6 + ch] : 0.0) + emissions[s];
    }
    for (int iq = 0; iq < 6; iq++) {
      int q = iq < 3 ? 4 + iq : 10 + (iq - 3);
      e_intron[i * 6 + iq] = emissions[q];
      S_intron[i * 6 + iq] = (i > 0 ? S_intron[(i - 1) * 6 + iq] : 0.0) + emissions[q];
    }
    int gold_state = c->label[c->offset + i];
    gold_score += emissions[gold_state];
    gradient[emission_index(&sites[i], gold_state)] += 1.0;
    if (i > 0) {
      int prev_gold = c->label[c->offset + i - 1];
      int gold_edge = pair_to_edge[prev_gold][gold_state];
      if (gold_edge < 0) anno_fail("Invalid transition in validated CRF training path");
      gold_score += tr[gold_edge] - emissions[gold_state];
      int kind = edge_kind[gold_edge], rev = edge_rev[gold_edge];
      if (kind >= 0) {
        for (int off = 0; off < BOUNDARY_KMER_COUNT; off++)
          gradient[boundary_kmer_index(&sites[i], rev, kind, off)] += 1.0;
        if (kind == 1) {
          int gpos = c->offset + i;
          float power[WAVE_COUNT];
          double deriv[DERIV_COUNT][WAVE_COUNT];
          for (int off = 0; off < BOUNDARY_CWT_COUNT; off++) {
            boundary_cwt_at(c, gpos, rev, off, power, deriv);
            int base = off * WAVE_COUNT;
            for (int sc = 0; sc < WAVE_COUNT; sc++) {
              gradient[OFF_BP + base + sc] += power[sc];
              for (int slot = 0; slot < DERIV_COUNT; slot++)
                gradient[derivative_index(slot, off, sc)] += deriv[slot][sc];
            }
          }
        }
      }
    }
  }
  for (int t = 0; t < length;) {
    int gs = c->label[c->offset + t];
    if (gs == 0) { t++; continue; }
    int u = t;
    if (is_cds_state(gs)) {
      int minus = gs >= 7;
      while (u + 1 < length && is_cds_state(c->label[c->offset + u + 1]) &&
             (c->label[c->offset + u + 1] >= 7) == minus) u++;
      int d = u - t + 1;
      gold_score += phi_exon(w, d);
      if (d <= CRF_EXON_LEN) gradient[OFF_EXON_LEN + d - 1] += 1.0;
      else {
        gradient[OFF_EXON_LEN + CRF_EXON_LEN - 1] += 1.0;
        gradient[OFF_EXON_TAIL] += (double)(d - CRF_EXON_LEN);
      }
    } else {
      while (u + 1 < length && c->label[c->offset + u + 1] == gs) u++;
      int d = u - t + 1;
      gold_score += phi_intron(w, d);
      if (d <= CRF_INTRON_LEN) gradient[OFF_INTRON_LEN + d - 1] += 1.0;
      else {
        gradient[OFF_INTRON_LEN + CRF_INTRON_LEN - 1] += 1.0;
        gradient[OFF_INTRON_TAIL] += (double)(d - CRF_INTRON_LEN);
      }
    }
    t = u + 1;
  }

  int max_exon_d = 3;
  for (int ch = 0; ch < 6; ch++) {
    int cur_min_t = 0;
    for (int u = 0; u < length; u++) {
      if (u < 2) {
        canon_h[u * 6 + ch] = -1;
        min_t[u * 6 + ch] = 0;
        continue;
      }
      int s0 = cds_state_at(u - 2, ch);
      int s1 = next_cds_state(s0);
      int s2 = next_cds_state(s1);
      int b2 = sites[u - 2].nuc, b1 = sites[u - 1].nuc, b0 = sites[u].nuc;
      int h2_found = -1;
      for (int idx0 = 0; idx0 < hidden_counts[s0]; idx0++) {
        int h0 = hidden_offsets[s0] + idx0;
        if (!initial_states[b2][h0]) continue;
        int h1 = grammar_next(h0, s1, b1);
        if (h1 < 0) continue;
        int h2 = grammar_next(h1, s2, b0);
        if (h2 >= 0) { h2_found = h2; break; }
      }
      canon_h[u * 6 + ch] = (int8_t)h2_found;
      if (u > 2) {
        int h_prev = canon_h[(u - 1) * 6 + ch];
        if (h_prev < 0 || grammar_next(h_prev, s2, b0) < 0)
          cur_min_t = u - 2;
      }
      if (h2_found < 0) cur_min_t = u - 1;
      min_t[u * 6 + ch] = cur_min_t;
      int d_span = u - cur_min_t + 1;
      if (d_span > max_exon_d) max_exon_d = d_span;
    }
  }
  work_ensure_exon_k(k, max_exon_d);
  double *K_exon_rev = k->K_exon_rev, *G_exon_rev = k->G_exon_rev;

  double max_phi_ex = -INFINITY;
  for (int d = 1; d <= max_exon_d; d++) {
    double v = phi_exon(w, d);
    if (!isfinite(v)) anno_fail("Non-finite exon duration weight");
    if (v > max_phi_ex) max_phi_ex = v;
  }
  for (int d = 1; d <= max_exon_d; d++)
    K_exon_rev[max_exon_d - d] = exp(phi_exon(w, d) - max_phi_ex);

  double max_phi_int = -INFINITY;
  for (int d = 1; d <= CRF_INTRON_LEN; d++) {
    double v = phi_intron(w, d);
    if (!isfinite(v)) anno_fail("Non-finite intron duration weight");
    if (v > max_phi_int) max_phi_int = v;
  }
  double K_intron_rev[CRF_INTRON_LEN];
  for (int d = 1; d <= CRF_INTRON_LEN; d++)
    K_intron_rev[CRF_INTRON_LEN - d] = exp(phi_intron(w, d) - max_phi_int);

  double phi_ex1 = phi_exon(w, 1), phi_ex2 = phi_exon(w, 2);
  double w_intron_max = w[OFF_INTRON_LEN + CRF_INTRON_LEN - 1];
  double w_intron_tail = w[OFF_INTRON_TAIL];
  if (!isfinite(w_intron_tail)) anno_fail("Non-finite intron tail weight");

  for (size_t i = 0; i < (size_t)length * HIDDEN_STATES; i++) {
    log_a_start[i] = -INFINITY;
    log_a_end[i] = -INFINITY;
  }
  for (size_t i = 0; i < (size_t)length * 6; i++) {
    log_exon_merged[i] = -INFINITY;
    X_exon[i] = -INFINITY;
  }
  for (size_t i = 0; i < (size_t)length * INTRON_CHANNELS; i++) {
    X_intron[i] = -INFINITY;
    log_tail[i] = -INFINITY;
    t_enter[i] = -INFINITY;
  }

  for (int b0 = 0; b0 < length; b0 += SM_BLK_SIZE) {
    int b_end = b0 + SM_BLK_SIZE < length ? b0 + SM_BLK_SIZE : length;
    double anchor_ex[6], anchor_int[6];
    for (int ch = 0; ch < 6; ch++) anchor_ex[ch] = S_exon[b0 * 6 + ch];
    for (int iq = 0; iq < 6; iq++) anchor_int[iq] = S_intron[b0 * 6 + iq];

    double m_ex[6], m_int[6];
    for (int ch = 0; ch < 6; ch++) {
      int t_earliest = b0;
      for (int u = b0; u < b_end; u++)
        if (u >= 2 && canon_h[u * 6 + ch] >= 0 && min_t[u * 6 + ch] < t_earliest)
          t_earliest = min_t[u * 6 + ch];
      double m = -INFINITY;
      if (b0 > 0) {
        double *le_prev = log_a_end + (size_t)(b0 - 1) * HIDDEN_STATES;
        for (int h = 0; h < HIDDEN_STATES; h++) if (le_prev[h] > m) m = le_prev[h];
      } else {
        m = 0.0;
      }
      for (int t = t_earliest; t < b0; t++) {
        double v = X_exon[t * 6 + ch] + anchor_ex[ch];
        if (v > m) m = v;
      }
      if (!isfinite(m)) m = 0.0;
      m_ex[ch] = m - anchor_ex[ch];
      for (int t = t_earliest; t < b0; t++) {
        double v = X_exon[t * 6 + ch] - m_ex[ch];
        exp_X_exon[(size_t)ch * e_stride + (size_t)t] = v > -745.0 ? exp(v) : 0.0;
      }
    }

    int t_int_earliest = b0 - CRF_INTRON_LEN + 1 > 0 ? b0 - CRF_INTRON_LEN + 1 : 0;
    for (int iq = 0; iq < 6; iq++) {
      double m = -INFINITY;
      if (b0 > 0) {
        double *le_prev = log_a_end + (size_t)(b0 - 1) * HIDDEN_STATES;
        for (int h = 0; h < HIDDEN_STATES; h++) if (le_prev[h] > m) m = le_prev[h];
      } else {
        m = 0.0;
      }
      for (int t = t_int_earliest; t < b0; t++) {
        for (int ic = 0; ic < INTRON_CHANNELS; ic++) {
          if (intron_ch_to_iq[ic] != iq) continue;
          double v = X_intron[t * INTRON_CHANNELS + ic] + anchor_int[iq];
          if (v > m) m = v;
        }
      }
      if (!isfinite(m)) m = 0.0;
      m_int[iq] = m - anchor_int[iq];
    }
    for (int ic = 0; ic < INTRON_CHANNELS; ic++) {
      int iq = intron_ch_to_iq[ic];
      double mi = m_int[iq];
      double * restrict ex_ptr = exp_X_intron + (size_t)ic * i_stride;
      for (int t = t_int_earliest; t < b0; t++) {
        double v = X_intron[t * INTRON_CHANNELS + ic] - mi;
        ex_ptr[t] = v > -745.0 ? exp(v) : 0.0;
      }
    }

    for (int u = b0; u < b_end; u++) {
      double *ls_u = log_a_start + (size_t)u * HIDDEN_STATES;
      double *le_u = log_a_end + (size_t)u * HIDDEN_STATES;

      if (u == 0) {
        int b = sites[0].nuc;
        for (int h = 0; h < HIDDEN_STATES; h++) {
          if (initial_allowed(b, h, full_start)) {
            int s = hidden_state[h];
            ls_u[h] = s == 0 ? w[emission_index(&sites[0], 0)]
                     : is_cds_state(s) ? e_exon[exon_start_channel(0, s)]
                     : e_intron[s <= 6 ? s - 4 : s - 10 + 3];
          }
        }
      } else {
        double *le_prev = log_a_end + (size_t)(u - 1) * HIDDEN_STATES;
        const double *psi_u = psi + (size_t)u * ntrans;
        int b = sites[u].nuc;
        const HiddenEdge *edges = hidden_edges[b];
        int ec = hidden_edge_count[b];
        for (int e = 0; e < ec; e++) {
          const HiddenEdge *step = edges + e;
          if (step->transition != 0 && edge_kind[step->transition] < 0) continue;
          double prev_val = le_prev[step->previous];
          if (prev_val == -INFINITY) continue;
          ls_u[step->current] = log_add(ls_u[step->current], prev_val + psi_u[step->transition]);
        }
      }

      if (u + 2 < length) {
        int b1 = sites[u + 1].nuc, b2 = sites[u + 2].nuc;
        for (int ch = 0; ch < 6; ch++) {
          int s0 = cds_state_at(u, ch);
          int s1 = next_cds_state(s0);
          int s2 = next_cds_state(s1);
          double m = -INFINITY;
          for (int idx0 = 0; idx0 < hidden_counts[s0]; idx0++) {
            int h0 = hidden_offsets[s0] + idx0;
            double v0 = ls_u[h0];
            if (v0 == -INFINITY) continue;
            int h1 = grammar_next(h0, s1, b1);
            if (h1 < 0) continue;
            int h2 = grammar_next(h1, s2, b2);
            if (h2 < 0) continue;
            m = log_add(m, v0);
          }
          log_exon_merged[u * 6 + ch] = m;
          double x = m == -INFINITY ? -INFINITY : m - S_exon[u * 6 + ch];
          X_exon[u * 6 + ch] = x;
          double v = x - m_ex[ch];
          exp_X_exon[(size_t)ch * e_stride + (size_t)u] = v > -745.0 ? exp(v) : 0.0;
        }
      } else {
        for (int ch = 0; ch < 6; ch++) exp_X_exon[(size_t)ch * e_stride + (size_t)u] = 0.0;
      }

      for (int h = 10; h < HIDDEN_STATES; h++) {
        int ic = intron_h_to_ch(h);
        if (ic < 0) continue;
        double v0 = ls_u[h];
        if (v0 == -INFINITY) continue;
        int iq = intron_ch_to_iq[ic];
        X_intron[u * INTRON_CHANNELS + ic] =
          log_add(X_intron[u * INTRON_CHANNELS + ic], v0 - S_intron[u * 6 + iq]);
      }
      for (int ic = 0; ic < INTRON_CHANNELS; ic++) {
        int iq = intron_ch_to_iq[ic];
        double x = X_intron[u * INTRON_CHANNELS + ic];
        double v = x - m_int[iq];
        exp_X_intron[(size_t)ic * i_stride + (size_t)u] = v > -745.0 ? exp(v) : 0.0;
      }

      le_u[0] = ls_u[0];

      for (int h0 = 1; h0 < HIDDEN_STATES; h0++) {
        int s0 = hidden_state[h0];
        if (!is_cds_state(s0)) continue;
        double v0 = ls_u[h0];
        if (v0 != -INFINITY) le_u[h0] = log_add(le_u[h0], v0 + phi_ex1);
      }
      if (u >= 1) {
        double *ls_prev = log_a_start + (size_t)(u - 1) * HIDDEN_STATES;
        int b_u = sites[u].nuc;
        for (int h0 = 1; h0 < HIDDEN_STATES; h0++) {
          int s0 = hidden_state[h0];
          if (!is_cds_state(s0)) continue;
          double v0 = ls_prev[h0];
          if (v0 == -INFINITY) continue;
          int s1 = next_cds_state(s0);
          int h1 = grammar_next(h0, s1, b_u);
          if (h1 < 0) continue;
          int ch = exon_start_channel(u - 1, s0);
          le_u[h1] = log_add(le_u[h1], v0 + e_exon[u * 6 + ch] + phi_ex2);
        }
      }
      if (u >= 2) {
        for (int ch = 0; ch < 6; ch++) {
          int h_end = canon_h[u * 6 + ch];
          if (h_end < 0) continue;
          int t_min = min_t[u * 6 + ch];
          int t_max = u - 2;
          if (t_min > t_max) continue;
          int n_len = t_max - t_min + 1;
          int d_at_tmin = u - t_min + 1;
          const double * restrict xs = exp_X_exon + ((size_t)ch * e_stride + (size_t)t_min);
          const double * restrict ks = K_exon_rev + (max_exon_d - d_at_tmin);
          double Z = 0.0;
          for (int i = 0; i < n_len; i++) Z += xs[i] * ks[i];
          if (Z > 0.0) {
            double log_v = m_ex[ch] + max_phi_ex + S_exon[u * 6 + ch] + log(Z);
            le_u[h_end] = log_add(le_u[h_end], log_v);
          }
        }
      }

      int max_d = u + 1 < CRF_INTRON_LEN ? u + 1 : CRF_INTRON_LEN;
      int t_min_int = u - max_d + 1;
      int n_len_int = max_d;
      int b_next = u + 1 < length ? sites[u + 1].nuc : -1;
      const double * restrict ks_int = K_intron_rev + (CRF_INTRON_LEN - max_d);
      for (int ic = 0; ic < INTRON_CHANNELS; ic++) {
        int iq = intron_ch_to_iq[ic];
        if (u >= CRF_INTRON_LEN) {
          int t_ent = u - CRF_INTRON_LEN;
          double x_ent = X_intron[t_ent * INTRON_CHANNELS + ic];
          double te = x_ent == -INFINITY ? -INFINITY
                      : x_ent + S_intron[(u - 1) * 6 + iq] + w_intron_max;
          t_enter[(u - 1) * INTRON_CHANNELS + ic] = te;
          double lt_prev = log_add(te, log_tail[(u - 1) * INTRON_CHANNELS + ic]);
          log_tail[u * INTRON_CHANNELS + ic] =
            lt_prev == -INFINITY ? -INFINITY : lt_prev + e_intron[u * 6 + iq] + w_intron_tail;
        }
        if (u < length - 1) {
          if (!intron_can_exit[b_next][ic]) continue;
        } else if (full_end) {
          continue;
        }
        const double * restrict xs = exp_X_intron + ((size_t)ic * i_stride + (size_t)t_min_int);
        double Z = 0.0;
        for (int i = 0; i < n_len_int; i++) Z += xs[i] * ks_int[i];
        double log_exp = Z > 0.0 ? m_int[iq] + max_phi_int + S_intron[u * 6 + iq] + log(Z) : -INFINITY;
        int h = intron_ch_rep[ic];
        le_u[h] = log_add(log_exp, log_tail[u * INTRON_CHANNELS + ic]);
      }

      if (u == length - 1 && full_end) {
        for (int h = 0; h < HIDDEN_STATES; h++)
          if (!terminal_allowed(h)) le_u[h] = -INFINITY;
      }
    }
  }

  double logz = -INFINITY;
  double *le_last = log_a_end + (size_t)(length - 1) * HIDDEN_STATES;
  for (int h = 0; h < HIDDEN_STATES; h++) logz = log_add(logz, le_last[h]);
  if (!isfinite(logz)) anno_fail("Invalid CRF partition at position %d", length);
  double likelihood = mass - logz + gold_score;

  memset(P_end, 0, (size_t)length * HIDDEN_STATES * sizeof(double));
  memset(P_start, 0, (size_t)length * HIDDEN_STATES * sizeof(double));
  memset(P_tail, 0, (size_t)length * INTRON_CHANNELS * sizeof(double));
  memset(Z_exon, 0, 6 * e_stride * sizeof(double));
  memset(Z_intron, 0, INTRON_CHANNELS * i_stride * sizeof(double));

  for (int h = 0; h < HIDDEN_STATES; h++)
    if (le_last[h] > -INFINITY)
      P_end[(size_t)(length - 1) * HIDDEN_STATES + h] = exp(le_last[h] - logz);

  int n_blocks = (length + SM_BLK_SIZE - 1) / SM_BLK_SIZE;
  for (int blk = n_blocks - 1; blk >= 0; blk--) {
    int b0 = blk * SM_BLK_SIZE;
    int b_end = b0 + SM_BLK_SIZE < length ? b0 + SM_BLK_SIZE : length;

    int t_ex_earliest[6];
    double m_ex[6], m_int[6];
    for (int ch = 0; ch < 6; ch++) {
      int te = b0;
      for (int u = b0; u < b_end; u++)
        if (u >= 2 && canon_h[u * 6 + ch] >= 0 && min_t[u * 6 + ch] < te)
          te = min_t[u * 6 + ch];
      t_ex_earliest[ch] = te;
      double m = -INFINITY;
      for (int t = te; t < b_end; t++)
        if (X_exon[t * 6 + ch] > m) m = X_exon[t * 6 + ch];
      if (!isfinite(m)) m = 0.0;
      m_ex[ch] = m;
      for (int t = te; t < b_end; t++) {
        double v = X_exon[t * 6 + ch] - m_ex[ch];
        exp_X_exon[(size_t)ch * e_stride + (size_t)t] = v > -745.0 ? exp(v) : 0.0;
        Z_exon[(size_t)ch * e_stride + (size_t)t] = 0.0;
      }
    }

    int t_int_earliest = b0 - CRF_INTRON_LEN + 1 > 0 ? b0 - CRF_INTRON_LEN + 1 : 0;
    for (int iq = 0; iq < 6; iq++) {
      double m = -INFINITY;
      for (int ic = 0; ic < INTRON_CHANNELS; ic++) {
        if (intron_ch_to_iq[ic] != iq) continue;
        for (int t = t_int_earliest; t < b_end; t++)
          if (X_intron[t * INTRON_CHANNELS + ic] > m) m = X_intron[t * INTRON_CHANNELS + ic];
      }
      if (!isfinite(m)) m = 0.0;
      m_int[iq] = m;
    }
    for (int ic = 0; ic < INTRON_CHANNELS; ic++) {
      int iq = intron_ch_to_iq[ic];
      double mi = m_int[iq];
      double * restrict ex_ptr = exp_X_intron + (size_t)ic * i_stride;
      double * restrict z_ptr = Z_intron + (size_t)ic * i_stride;
      for (int t = t_int_earliest; t < b_end; t++) {
        double v = X_intron[t * INTRON_CHANNELS + ic] - mi;
        ex_ptr[t] = v > -745.0 ? exp(v) : 0.0;
        z_ptr[t] = 0.0;
      }
    }

    double G_intron_rev[CRF_INTRON_LEN] = {0};
    memset(G_exon_rev, 0, (size_t)max_exon_d * sizeof(double));

    for (int u = b_end - 1; u >= b0; u--) {
      double *pe_u = P_end + (size_t)u * HIDDEN_STATES;
      double *ps_u = P_start + (size_t)u * HIDDEN_STATES;
      double *le_u = log_a_end + (size_t)u * HIDDEN_STATES;
      double *ls_u = log_a_start + (size_t)u * HIDDEN_STATES;

      if (u + 1 < length) {
        double *ps_next = P_start + (size_t)(u + 1) * HIDDEN_STATES;
        double *ls_next = log_a_start + (size_t)(u + 1) * HIDDEN_STATES;
        const double *psi_next = psi + (size_t)(u + 1) * ntrans;
        int b_next = sites[u + 1].nuc;
        const HiddenEdge *edges = hidden_edges[b_next];
        int ec = hidden_edge_count[b_next];
        double edge_marg[STATES * STATES] = {0};
        for (int e = 0; e < ec; e++) {
          const HiddenEdge *step = edges + e;
          if (step->transition != 0 && edge_kind[step->transition] < 0) continue;
          double ps_n = ps_next[step->current];
          if (ps_n <= 0.0) continue;
          double le_p = le_u[step->previous];
          if (le_p == -INFINITY) continue;
          double w_edge = ps_n * exp(le_p + psi_next[step->transition] - ls_next[step->current]);
          pe_u[step->previous] += w_edge;
          edge_marg[step->transition] += w_edge;
        }
        double b_coeff[4][2] = {{0}};
        int has_boundary = 0;
        for (int edge = 0; edge < ntrans; edge++) {
          double m = edge_marg[edge];
          if (m != 0.0) {
            int kind = edge_kind[edge], rev = edge_rev[edge];
            if (kind >= 0) { b_coeff[kind][rev] -= m; has_boundary = 1; }
          }
        }
        if (has_boundary) {
          int genome_position = c->offset + u + 1;
          for (int rev = 0; rev < 2; rev++) {
            if (b_coeff[0][rev] == 0 && b_coeff[1][rev] == 0 &&
                b_coeff[2][rev] == 0 && b_coeff[3][rev] == 0) continue;
            float power[BOUNDARY_CWT_COUNT][WAVE_COUNT];
            double derivative[BOUNDARY_CWT_COUNT][DERIV_COUNT][WAVE_COUNT];
            for (int offset = 0; b_coeff[1][rev] != 0 && offset < BOUNDARY_CWT_COUNT; offset++)
              boundary_cwt_at(c, genome_position, rev, offset, power[offset], derivative[offset]);
            for (int kind = 0; kind < 4; kind++) {
              double val = b_coeff[kind][rev];
              if (val == 0.0) continue;
              for (int offset = 0; offset < BOUNDARY_KMER_COUNT; offset++)
                gradient[boundary_kmer_index(&sites[u + 1], rev, kind, offset)] += val;
              for (int offset = 0; kind == 1 && offset < BOUNDARY_CWT_COUNT; offset++) {
                int base = offset * WAVE_COUNT;
                double *gbp = gradient + OFF_BP + base;
                for (int scale = 0; scale < WAVE_COUNT; scale++) {
                  gbp[scale] += val * power[offset][scale];
                  for (int slot = 0; slot < DERIV_COUNT; slot++)
                    gradient[derivative_index(slot, offset, scale)] +=
                      val * derivative[offset][slot][scale];
                }
              }
            }
          }
        }
      }

      ps_u[0] += pe_u[0];

      for (int h0 = 1; h0 < HIDDEN_STATES; h0++) {
        int s0 = hidden_state[h0];
        if (!is_cds_state(s0)) continue;
        double pe = pe_u[h0], v0 = ls_u[h0];
        if (pe <= 0.0 || v0 == -INFINITY) continue;
        double post = pe * exp(v0 + phi_ex1 - le_u[h0]);
        ps_u[h0] += post;
        gradient[OFF_EXON_LEN + 0] -= post;
      }

      if (u >= 1) {
        double *ls_prev = log_a_start + (size_t)(u - 1) * HIDDEN_STATES;
        double *ps_prev = P_start + (size_t)(u - 1) * HIDDEN_STATES;
        int b_u = sites[u].nuc;
        for (int h0 = 1; h0 < HIDDEN_STATES; h0++) {
          int s0 = hidden_state[h0];
          if (!is_cds_state(s0)) continue;
          double v0 = ls_prev[h0];
          if (v0 == -INFINITY) continue;
          int s1 = next_cds_state(s0);
          int h1 = grammar_next(h0, s1, b_u);
          if (h1 < 0) continue;
          double pe = pe_u[h1];
          if (pe <= 0.0) continue;
          int ch = exon_start_channel(u - 1, s0);
          double post = pe * exp(v0 + e_exon[u * 6 + ch] + phi_ex2 - le_u[h1]);
          ps_prev[h0] += post;
          gradient[OFF_EXON_LEN + 1] -= post;
        }
      }

      if (u >= 2) {
        for (int ch = 0; ch < 6; ch++) {
          int h_end = canon_h[u * 6 + ch];
          if (h_end < 0) continue;
          double pe = pe_u[h_end];
          if (pe <= 0.0) continue;
          int t_min = min_t[u * 6 + ch];
          int t_max = u - 2;
          if (t_min > t_max) continue;
          int n_len = t_max - t_min + 1;
          int d_at_tmin = u - t_min + 1;
          double coeff = pe * exp(m_ex[ch] + max_phi_ex + S_exon[u * 6 + ch] - le_u[h_end]);
          const double * restrict xs = exp_X_exon + ((size_t)ch * e_stride + (size_t)t_min);
          double * restrict zs = Z_exon + ((size_t)ch * e_stride + (size_t)t_min);
          const double * restrict ks = K_exon_rev + (max_exon_d - d_at_tmin);
          double * restrict gs = G_exon_rev + (max_exon_d - d_at_tmin);
          for (int i = 0; i < n_len; i++) {
            zs[i] += coeff * ks[i];
            gs[i] += coeff * xs[i];
          }
        }
      }

      int max_d = u + 1 < CRF_INTRON_LEN ? u + 1 : CRF_INTRON_LEN;
      int t_min_int = u - max_d + 1;
      int n_len_int = max_d;
      const double * restrict ks_int = K_intron_rev + (CRF_INTRON_LEN - max_d);
      double * restrict gs_int = G_intron_rev + (CRF_INTRON_LEN - max_d);

      for (int ic = 0; ic < INTRON_CHANNELS; ic++) {
        int h = intron_ch_rep[ic];
        int iq = intron_ch_to_iq[ic];
        double p_end = pe_u[h];
        if (p_end > 0.0 && le_u[h] > -INFINITY) {
          double lt = log_tail[u * INTRON_CHANNELS + ic];
          if (lt > -INFINITY)
            P_tail[u * INTRON_CHANNELS + ic] += p_end * exp(lt - le_u[h]);
          double coeff = p_end * exp(m_int[iq] + max_phi_int + S_intron[u * 6 + iq] - le_u[h]);
          if (coeff > 0.0) {
            const double * restrict xs = exp_X_intron + ((size_t)ic * i_stride + (size_t)t_min_int);
            double * restrict zs = Z_intron + ((size_t)ic * i_stride + (size_t)t_min_int);
            for (int i = 0; i < n_len_int; i++) {
              zs[i] += coeff * ks_int[i];
              gs_int[i] += coeff * xs[i];
            }
          }
        }

        double pt = P_tail[u * INTRON_CHANNELS + ic];
        if (pt > 0.0 && u >= CRF_INTRON_LEN) {
          gradient[OFF_INTRON_TAIL] -= pt;
          double te = t_enter[(u - 1) * INTRON_CHANNELS + ic];
          double lt_prev = log_tail[(u - 1) * INTRON_CHANNELS + ic];
          double lt_tot = log_add(te, lt_prev);
          if (lt_tot > -INFINITY) {
            double p_enter = te > -INFINITY ? pt * exp(te - lt_tot) : 0.0;
            double p_cont  = lt_prev > -INFINITY ? pt * exp(lt_prev - lt_tot) : 0.0;
            P_tail[(u - 1) * INTRON_CHANNELS + ic] += p_cont;
            if (p_enter > 0.0) {
              gradient[OFF_INTRON_LEN + CRF_INTRON_LEN - 1] -= p_enter;
              int t_ent = u - CRF_INTRON_LEN;
              double x_ent = X_intron[t_ent * INTRON_CHANNELS + ic];
              for (int h0 = 10; h0 < HIDDEN_STATES; h0++) {
                if (intron_h_to_ch(h0) != ic) continue;
                double v0 = log_a_start[(size_t)t_ent * HIDDEN_STATES + h0];
                if (v0 == -INFINITY) continue;
                P_start[(size_t)t_ent * HIDDEN_STATES + h0] +=
                  p_enter * exp(v0 - S_intron[t_ent * 6 + iq] - x_ent);
              }
            }
          }
        }
      }

      if (u + 2 < length) {
        int b1 = sites[u + 1].nuc, b2 = sites[u + 2].nuc;
        for (int ch = 0; ch < 6; ch++) {
          double z_val = Z_exon[(size_t)ch * e_stride + (size_t)u];
          double p_m = exp_X_exon[(size_t)ch * e_stride + (size_t)u] * z_val;
          if (p_m <= 0.0) continue;
          double lm = log_exon_merged[u * 6 + ch];
          if (lm == -INFINITY) continue;
          int s0 = cds_state_at(u, ch);
          int s1 = next_cds_state(s0);
          int s2 = next_cds_state(s1);
          for (int idx0 = 0; idx0 < hidden_counts[s0]; idx0++) {
            int h0 = hidden_offsets[s0] + idx0;
            double v0 = ls_u[h0];
            if (v0 == -INFINITY) continue;
            int h1 = grammar_next(h0, s1, b1);
            if (h1 < 0) continue;
            int h2 = grammar_next(h1, s2, b2);
            if (h2 < 0) continue;
            ps_u[h0] += p_m * exp(v0 - lm);
          }
        }
      }
      for (int ic = 0; ic < INTRON_CHANNELS; ic++) {
        double z_val = Z_intron[(size_t)ic * i_stride + (size_t)u];
        double p_tot = exp_X_intron[(size_t)ic * i_stride + (size_t)u] * z_val;
        if (p_tot <= 0.0) continue;
        int iq = intron_ch_to_iq[ic];
        double x_tot = X_intron[u * INTRON_CHANNELS + ic];
        for (int h0 = 10; h0 < HIDDEN_STATES; h0++) {
          if (intron_h_to_ch(h0) != ic) continue;
          double v0 = ls_u[h0];
          if (v0 == -INFINITY) continue;
          ps_u[h0] += p_tot * exp(v0 - S_intron[u * 6 + iq] - x_tot);
        }
      }
    }

    for (int d = 1; d <= CRF_INTRON_LEN; d++) {
      int j = CRF_INTRON_LEN - d;
      gradient[OFF_INTRON_LEN + d - 1] -= G_intron_rev[j] * K_intron_rev[j];
    }
    for (int d = 3; d <= max_exon_d; d++) {
      int j = max_exon_d - d;
      double post = G_exon_rev[j] * K_exon_rev[j];
      if (d <= CRF_EXON_LEN) {
        gradient[OFF_EXON_LEN + d - 1] -= post;
      } else {
        gradient[OFF_EXON_LEN + CRF_EXON_LEN - 1] -= post;
        gradient[OFF_EXON_TAIL] -= post * (double)(d - CRF_EXON_LEN);
      }
    }

    for (int ch = 0; ch < 6; ch++) {
      for (int t = t_ex_earliest[ch]; t < b0; t++) {
        double z_val = Z_exon[(size_t)ch * e_stride + (size_t)t];
        double p_m = exp_X_exon[(size_t)ch * e_stride + (size_t)t] * z_val;
        if (p_m <= 0.0) continue;
        double lm = log_exon_merged[t * 6 + ch];
        if (lm == -INFINITY) continue;
        int b1 = sites[t + 1].nuc, b2 = sites[t + 2].nuc;
        int s0 = cds_state_at(t, ch);
        int s1 = next_cds_state(s0);
        int s2 = next_cds_state(s1);
        for (int idx0 = 0; idx0 < hidden_counts[s0]; idx0++) {
          int h0 = hidden_offsets[s0] + idx0;
          double v0 = log_a_start[(size_t)t * HIDDEN_STATES + h0];
          if (v0 == -INFINITY) continue;
          int h1 = grammar_next(h0, s1, b1);
          if (h1 < 0) continue;
          int h2 = grammar_next(h1, s2, b2);
          if (h2 < 0) continue;
          P_start[(size_t)t * HIDDEN_STATES + h0] += p_m * exp(v0 - lm);
        }
      }
    }
    for (int ic = 0; ic < INTRON_CHANNELS; ic++) {
      int iq = intron_ch_to_iq[ic];
      for (int t = t_int_earliest; t < b0; t++) {
        double z_val = Z_intron[(size_t)ic * i_stride + (size_t)t];
        double p_tot = exp_X_intron[(size_t)ic * i_stride + (size_t)t] * z_val;
        if (p_tot <= 0.0) continue;
        double x_tot = X_intron[t * INTRON_CHANNELS + ic];
        for (int h0 = 10; h0 < HIDDEN_STATES; h0++) {
          if (intron_h_to_ch(h0) != ic) continue;
          double v0 = log_a_start[(size_t)t * HIDDEN_STATES + h0];
          if (v0 == -INFINITY) continue;
          P_start[(size_t)t * HIDDEN_STATES + h0] += p_tot * exp(v0 - S_intron[t * 6 + iq] - x_tot);
        }
      }
    }
  }

  double gamma_ex[6] = {0}, gamma_int[6] = {0};
  for (int i = length - 1; i >= 0; i--) {
    const double *pe_i = P_end + (size_t)i * HIDDEN_STATES;
    for (int h = 1; h < HIDDEN_STATES; h++) {
      double pe = pe_i[h];
      if (pe == 0.0) continue;
      int s = hidden_state[h];
      if (is_cds_state(s)) {
        int ch = s <= 3 ? (i - (s - 1) + 3) % 3 : 3 + (i + (s - 7)) % 3;
        gamma_ex[ch] += pe;
      } else {
        int iq = s <= 6 ? s - 4 : s - 10 + 3;
        gamma_int[iq] += pe;
      }
    }
    if (i + 1 < length) {
      const double *ps_next = P_start + (size_t)(i + 1) * HIDDEN_STATES;
      for (int h = 1; h < HIDDEN_STATES; h++) {
        double ps = ps_next[h];
        if (ps == 0.0) continue;
        int s = hidden_state[h];
        if (is_cds_state(s)) {
          int ch = exon_start_channel(i + 1, s);
          gamma_ex[ch] -= ps;
        } else {
          int iq = s <= 6 ? s - 4 : s - 10 + 3;
          gamma_int[iq] -= ps;
        }
      }
    }
    gradient[emission_index(&sites[i], 0)] -= pe_i[0];
    for (int ch = 0; ch < 6; ch++) {
      int s = cds_state_at(i, ch);
      gradient[emission_index(&sites[i], s)] -= gamma_ex[ch];
    }
    for (int iq = 0; iq < 6; iq++) {
      int q = iq < 3 ? 4 + iq : 10 + (iq - 3);
      gradient[emission_index(&sites[i], q)] -= gamma_int[iq];
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
  for (int c = 0; c < count; c++) {
    for (int i = 0; i < contigs[c].length; i++) {
      Site site = site_at(contigs + c, i);
      int state = contigs[c].label[contigs[c].offset + i];
      scale[emission_index(&site, state)] += 1;
      if (!i) continue;
      int previous = contigs[c].label[contigs[c].offset + i - 1];
      int reverse = previous >= 7 || state >= 7, kind = boundary_kind(previous, state);
      if (kind >= 0) {
        for (int offset = 0; offset < BOUNDARY_KMER_COUNT; offset++)
          scale[boundary_kmer_index(&site, reverse, kind, offset)] += 1;
        int genome_position = contigs[c].offset + i;
        for (int offset = 0; kind == 1 && offset < BOUNDARY_CWT_COUNT; offset++) {
          float power[WAVE_COUNT];
          double derivative[DERIV_COUNT][WAVE_COUNT];
          boundary_cwt_at(contigs + c, genome_position, reverse, offset, power, derivative);
          int base = offset * WAVE_COUNT;
          for (int k = 0; k < WAVE_COUNT; k++) {
            scale[OFF_BP + base + k] += (double)power[k] * power[k];
            for (int slot = 0; slot < DERIV_COUNT; slot++)
              scale[derivative_index(slot, offset, k)] +=
                derivative[slot][k] * derivative[slot][k];
          }
        }
      }
    }
    const uint8_t *labels = contigs[c].label + contigs[c].offset;
    for (int t = 0; t < contigs[c].length;) {
      int s0 = labels[t];
      if (!s0) { t++; continue; }
      int u = t;
      if (is_cds_state(s0)) {
        while (u + 1 < contigs[c].length && is_cds_state(labels[u + 1])) u++;
        int d = u - t + 1;
        if (d <= CRF_EXON_LEN) {
          scale[OFF_EXON_LEN + d - 1] += 1.0;
        } else {
          scale[OFF_EXON_LEN + CRF_EXON_LEN - 1] += 1.0;
          double tail = (double)(d - CRF_EXON_LEN);
          scale[OFF_EXON_TAIL] += tail * tail;
        }
      } else {
        while (u + 1 < contigs[c].length && labels[u + 1] == s0) u++;
        int d = u - t + 1;
        if (d <= CRF_INTRON_LEN) {
          scale[OFF_INTRON_LEN + d - 1] += 1.0;
        } else {
          scale[OFF_INTRON_LEN + CRF_INTRON_LEN - 1] += 1.0;
          double tail = (double)(d - CRF_INTRON_LEN);
          scale[OFF_INTRON_TAIL] += tail * tail;
        }
      }
      t = u + 1;
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

/* Exact whole-contig semi-Markov MAP decoding. */
void crf_decode(const double *w, const Contig *c, uint8_t *path) {
  int length = c->length;
  int intron_ring = CRF_INTRON_LEN + 2;

  double *X_exon = anno_alloc((size_t)length * 6, sizeof(double));
  uint8_t *bp_exon_prev_h = anno_alloc((size_t)length * 6, sizeof(uint8_t));
  double *X_intron_ring = anno_alloc((size_t)intron_ring * INTRON_CHANNELS, sizeof(double));
  uint8_t *bp_intron_prev_h = anno_alloc((size_t)intron_ring * INTRON_CHANNELS, sizeof(uint8_t));
  int32_t *bp_t = anno_alloc((size_t)length * HIDDEN_STATES, sizeof(int32_t));
  uint8_t *bp_prev_h = anno_alloc((size_t)length * HIDDEN_STATES, sizeof(uint8_t));

  int nuc_ring[4] = {0};
  double emit_ring[4][STATES];
  double trans_ring[4][STATES * STATES];
  double v_start_ring[4][HIDDEN_STATES];
  uint8_t vsp_ring[4][HIDDEN_STATES];
  double v_end[HIDDEN_STATES];
  for (int h = 0; h < HIDDEN_STATES; h++) v_end[h] = -INFINITY;

  double S_exon[6] = {0}, S_intron[6] = {0};
  int min_t[6] = {0};
  double val_tail[INTRON_CHANNELS];
  int32_t tail_t[INTRON_CHANNELS];
  uint8_t tail_prev_h[INTRON_CHANNELS];
  for (int ic = 0; ic < INTRON_CHANNELS; ic++) {
    val_tail[ic] = -INFINITY;
    tail_t[ic] = 0;
    tail_prev_h[ic] = 0;
  }

  double phi_ex1 = phi_exon(w, 1), phi_ex2 = phi_exon(w, 2);

  for (int u = 0; u < length; u++) {
    int r0 = u & 3;
    Site site = site_scores(w, c, u, emit_ring[r0], trans_ring[r0]);
    int b0 = site.nuc;
    nuc_ring[r0] = b0;
    const double *emit_u = emit_ring[r0];
    const double *trans_u = trans_ring[r0];
    double *vs_u = v_start_ring[r0];
    uint8_t *vs_prev_h = vsp_ring[r0];
    uint8_t *bp_ph_u = bp_prev_h + (size_t)u * HIDDEN_STATES;
    int32_t *bt_u = bp_t + (size_t)u * HIDDEN_STATES;

    for (int h = 0; h < HIDDEN_STATES; h++) {
      vs_u[h] = -INFINITY;
      vs_prev_h[h] = 0;
      bt_u[h] = u;
    }

    if (u == 0) {
      for (int h = 0; h < HIDDEN_STATES; h++)
        if (initial_allowed(b0, h, 1))
          vs_u[h] = emit_u[hidden_state[h]];
    } else {
      const HiddenEdge *edges = hidden_edges[b0];
      int n_edges = hidden_edge_count[b0];
      for (int e = 0; e < n_edges; e++) {
        int prev_h = edges[e].previous;
        int cur_h = edges[e].current;
        int edge_idx = edges[e].transition;
        if (edge_idx != 0 && edge_kind[edge_idx] < 0) continue;
        double pv = v_end[prev_h];
        if (pv == -INFINITY) continue;
        double cand = pv + trans_u[edge_idx] + emit_u[hidden_state[cur_h]];
        if (cand > vs_u[cur_h]) {
          vs_u[cur_h] = cand;
          vs_prev_h[cur_h] = (uint8_t)prev_h;
        }
      }
    }

    for (int ch = 0; ch < 6; ch++) {
      int s = cds_state_at(u, ch);
      if (u == 0) {
        S_exon[ch] = emit_u[s];
      } else {
        int s_prev = cds_state_at(u - 1, ch);
        S_exon[ch] += emit_u[s] + trans_u[pair_to_edge[s_prev][s]];
      }
    }

    double step_int[6];
    for (int iq = 0; iq < 6; iq++) {
      int q = iq < 3 ? 4 + iq : 10 + (iq - 3);
      step_int[iq] = u == 0 ? emit_u[q] : emit_u[q] + trans_u[pair_to_edge[q][q]];
      S_intron[iq] += step_int[iq];
    }

    if (u >= 2) {
      int t = u - 2;
      int rt = t & 3, r1 = (u - 1) & 3;
      int b1 = nuc_ring[r1];
      const double *emit_1 = emit_ring[r1];
      const double *trans_1 = trans_ring[r1];
      const double *vs_t = v_start_ring[rt];
      const uint8_t *vsp_t = vsp_ring[rt];
      for (int ch = 0; ch < 6; ch++) {
        int s0 = cds_state_at(t, ch);
        int s1 = next_cds_state(s0);
        int s2 = next_cds_state(s1);
        double step12 = emit_1[s1] + trans_1[pair_to_edge[s0][s1]]
                      + emit_u[s2] + trans_u[pair_to_edge[s1][s2]];
        double best_v0 = -INFINITY;
        uint8_t best_ph = 0;
        for (int idx0 = 0; idx0 < hidden_counts[s0]; idx0++) {
          int h0 = hidden_offsets[s0] + idx0;
          double v0 = vs_t[h0];
          if (v0 == -INFINITY) continue;
          int h1 = grammar_next(h0, s1, b1);
          if (h1 < 0) continue;
          int h2 = grammar_next(h1, s2, b0);
          if (h2 < 0) continue;
          if (v0 > best_v0) {
            best_v0 = v0;
            best_ph = vsp_t[h0];
          }
        }
        X_exon[(size_t)t * 6 + ch] = best_v0 > -INFINITY ? best_v0 + step12 - S_exon[ch] : -INFINITY;
        bp_exon_prev_h[(size_t)t * 6 + ch] = best_ph;
      }
    }

    int ir_u = u % intron_ring;
    for (int ic = 0; ic < INTRON_CHANNELS; ic++) {
      int iq = intron_ch_to_iq[ic];
      double best_v0 = -INFINITY;
      uint8_t best_ph = 0;
      for (int h0 = 10; h0 < HIDDEN_STATES; h0++) {
        if (intron_h_to_ch(h0) != ic) continue;
        double v0 = vs_u[h0];
        if (v0 > best_v0) {
          best_v0 = v0;
          best_ph = vs_prev_h[h0];
        }
      }
      X_intron_ring[ir_u * INTRON_CHANNELS + ic] = best_v0 > -INFINITY ? best_v0 - S_intron[iq] : -INFINITY;
      bp_intron_prev_h[ir_u * INTRON_CHANNELS + ic] = best_ph;
    }

    double next_v_end[HIDDEN_STATES];
    uint8_t next_bp_ph[HIDDEN_STATES];
    int32_t next_bp_t[HIDDEN_STATES];
    for (int h = 0; h < HIDDEN_STATES; h++) {
      next_v_end[h] = -INFINITY;
      next_bp_ph[h] = vs_prev_h[h];
      next_bp_t[h] = u;
    }

    next_v_end[0] = vs_u[0];

    for (int h0 = 1; h0 < HIDDEN_STATES; h0++) {
      int s0 = hidden_state[h0];
      if (!is_cds_state(s0)) continue;
      double v0 = vs_u[h0];
      if (v0 == -INFINITY) continue;
      double cand = v0 + phi_ex1;
      if (cand > next_v_end[h0]) {
        next_v_end[h0] = cand;
        next_bp_t[h0] = u;
        next_bp_ph[h0] = vs_prev_h[h0];
      }
    }

    if (u >= 1) {
      int r1 = (u - 1) & 3;
      const double *vs_prev = v_start_ring[r1];
      const uint8_t *vsp_prev = vsp_ring[r1];
      for (int h0 = 1; h0 < HIDDEN_STATES; h0++) {
        int s0 = hidden_state[h0];
        if (!is_cds_state(s0)) continue;
        double v0 = vs_prev[h0];
        if (v0 == -INFINITY) continue;
        int s1 = next_cds_state(s0);
        int h1 = grammar_next(h0, s1, b0);
        if (h1 < 0) continue;
        double cand = v0 + emit_u[s1] + trans_u[pair_to_edge[s0][s1]] + phi_ex2;
        if (cand > next_v_end[h1]) {
          next_v_end[h1] = cand;
          next_bp_t[h1] = u - 1;
          next_bp_ph[h1] = vsp_prev[h0];
        }
      }
    }

    if (u >= 2) {
      int b2 = nuc_ring[(u - 2) & 3], b1 = nuc_ring[(u - 1) & 3];
      for (int ch = 0; ch < 6; ch++) {
        int s0 = cds_state_at(u - 2, ch);
        int s1 = next_cds_state(s0);
        int s2 = next_cds_state(s1);
        int h_end = -1;
        for (int idx0 = 0; idx0 < hidden_counts[s0]; idx0++) {
          int h0 = hidden_offsets[s0] + idx0;
          if (!initial_states[b2][h0]) continue;
          int h1 = grammar_next(h0, s1, b1);
          if (h1 < 0) continue;
          int h2 = grammar_next(h1, s2, b0);
          if (h2 >= 0) { h_end = h2; break; }
        }
        if (h_end >= 0) {
          int t_min = min_t[ch];
          double s_u = S_exon[ch];
          double best_val = next_v_end[h_end];
          int32_t best_t = next_bp_t[h_end];
          uint8_t best_ph = next_bp_ph[h_end];
          for (int t = t_min; t <= u - 2; t++) {
            double x_t = X_exon[(size_t)t * 6 + ch];
            if (x_t == -INFINITY) continue;
            double cand = s_u + x_t + phi_exon(w, u - t + 1);
            if (cand > best_val) {
              best_val = cand;
              best_t = t;
              best_ph = bp_exon_prev_h[(size_t)t * 6 + ch];
            }
          }
          next_v_end[h_end] = best_val;
          next_bp_t[h_end] = best_t;
          next_bp_ph[h_end] = best_ph;
        }
        if (b2 < 4 && b1 < 4 && b0 < 4) {
          int codon = (b2 << 4) | (b1 << 2) | b0;
          if (s2 == 3 && (codon == 48 || codon == 50 || codon == 56))
            min_t[ch] = u - 1;
          else if (s2 == 7 && (codon == 60 || codon == 28 || codon == 52))
            min_t[ch] = u - 1;
        }
      }
    }

    int max_d = u + 1 < CRF_INTRON_LEN ? u + 1 : CRF_INTRON_LEN;
    for (int ic = 0; ic < INTRON_CHANNELS; ic++) {
      int iq = intron_ch_to_iq[ic];
      int h = intron_ch_rep[ic];

      if (u >= CRF_INTRON_LEN) {
        double cand_cont = val_tail[ic] > -INFINITY ? val_tail[ic] + step_int[iq] + w[OFF_INTRON_TAIL] : -INFINITY;
        int t_ent = u - CRF_INTRON_LEN;
        int ir_ent = t_ent % intron_ring;
        double x_ent = X_intron_ring[ir_ent * INTRON_CHANNELS + ic];
        double cand_enter = x_ent > -INFINITY
          ? S_intron[iq] + x_ent + w[OFF_INTRON_LEN + CRF_INTRON_LEN - 1] + w[OFF_INTRON_TAIL]
          : -INFINITY;
        if (cand_enter >= cand_cont && cand_enter > -INFINITY) {
          val_tail[ic] = cand_enter;
          tail_t[ic] = t_ent;
          tail_prev_h[ic] = bp_intron_prev_h[ir_ent * INTRON_CHANNELS + ic];
        } else {
          val_tail[ic] = cand_cont;
        }
      }

      double best_val = val_tail[ic];
      int32_t best_t = tail_t[ic];
      uint8_t best_ph = tail_prev_h[ic];
      double s_u = S_intron[iq];
      for (int d = 1; d <= max_d; d++) {
        int t = u - d + 1;
        int ir_t = t % intron_ring;
        double x_t = X_intron_ring[ir_t * INTRON_CHANNELS + ic];
        if (x_t == -INFINITY) continue;
        double cand = s_u + x_t + w[OFF_INTRON_LEN + d - 1];
        if (cand > best_val) {
          best_val = cand;
          best_t = t;
          best_ph = bp_intron_prev_h[ir_t * INTRON_CHANNELS + ic];
        }
      }
      next_v_end[h] = best_val;
      next_bp_t[h] = best_t;
      next_bp_ph[h] = best_ph;
    }

    memcpy(v_end, next_v_end, sizeof(v_end));
    memcpy(bt_u, next_bp_t, sizeof(next_bp_t));
    memcpy(bp_ph_u, next_bp_ph, sizeof(next_bp_ph));
  }

  int best_h = 0;
  for (int h = 1; h < HIDDEN_STATES; h++)
    if (terminal_allowed(h) && v_end[h] > v_end[best_h]) best_h = h;
  if (!isfinite(v_end[best_h])) anno_fail("No complete prediction path");

  int u = length - 1;
  int h_end = best_h;
  while (u >= 0) {
    int t = bp_t[(size_t)u * HIDDEN_STATES + h_end];
    int prev_h = bp_prev_h[(size_t)u * HIDDEN_STATES + h_end];
    int s_end = hidden_state[h_end];
    if (s_end == 0) {
      path[u] = 0;
    } else if (is_cds_state(s_end)) {
      int ch = s_end <= 3 ? (u - (s_end - 1) + 3) % 3 : 3 + (u + (s_end - 7)) % 3;
      for (int p = t; p <= u; p++)
        path[p] = (uint8_t)cds_state_at(p, ch);
    } else {
      for (int p = t; p <= u; p++)
        path[p] = (uint8_t)s_end;
    }
    h_end = prev_h;
    u = t - 1;
  }

  free(X_exon);
  free(bp_exon_prev_h);
  free(X_intron_ring);
  free(bp_intron_prev_h);
  free(bp_t);
  free(bp_prev_h);
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
    for (int offset = 0; offset < BOUNDARY_KMER_COUNT; offset++) {
      int disp = boundary_kmer_offsets[offset];
      for (int k = 0; k < K2; k++) {
        int idx = OFF_SPLICE_START_STOP + (kind * BOUNDARY_KMER_COUNT + offset) * K2 + k;
        format_kmer(k, 2, kbuf);
        fprintf(file, "%d\tboundary_k2\t%s\toffset=%+d,k2=%s\t%.9g\n",
                idx, kind_names[kind], disp, kbuf, (float)w[idx]);
      }
    }
  }
  for (int offset = 0; offset < BOUNDARY_CWT_COUNT; offset++) {
    int disp = boundary_cwt_offsets[offset];
    int base = offset * WAVE_COUNT;
    for (int sc = 0; sc < WAVE_COUNT; sc++) {
      int idx = OFF_BP + base + sc;
      fprintf(file, "%d\tboundary_cwt_power\t%s\toffset=%+d,scale=%d\t%.9g\n",
              idx, kind_names[1], disp, wave_sizes[sc], (float)w[idx]);
    }
  }
  for (int slot = 0; slot < DERIV_COUNT; slot++) {
    for (int offset = 0; offset < BOUNDARY_CWT_COUNT; offset++) {
      int disp = boundary_cwt_offsets[offset];
      for (int sc = 0; sc < WAVE_COUNT; sc++) {
        int idx = derivative_index(slot, offset, sc);
        fprintf(file, "%d\tboundary_cwt_deriv%d\t%s\toffset=%+d,scale=%d\t%.9g\n",
                idx, derivative_orders[slot], kind_names[1], disp, wave_sizes[sc], (float)w[idx]);
      }
    }
  }
  for (int d = 1; d <= CRF_EXON_LEN; d++) {
    int idx = OFF_EXON_LEN + d - 1;
    fprintf(file, "%d\texon_length\tCDS\tlen=%d\t%.9g\n", idx, d, (float)w[idx]);
  }
  fprintf(file, "%d\texon_length_tail\tCDS\tslope_gt_%d\t%.9g\n",
          OFF_EXON_TAIL, CRF_EXON_LEN, (float)w[OFF_EXON_TAIL]);
  for (int d = 1; d <= CRF_INTRON_LEN; d++) {
    int idx = OFF_INTRON_LEN + d - 1;
    fprintf(file, "%d\tintron_length\tintron\tlen=%d\t%.9g\n", idx, d, (float)w[idx]);
  }
  fprintf(file, "%d\tintron_length_tail\tintron\tslope_gt_%d\t%.9g\n",
          OFF_INTRON_TAIL, CRF_INTRON_LEN, (float)w[OFF_INTRON_TAIL]);

  if (fclose(file)) anno_fail("Cannot finish writing weights text: %s", path);
}

#endif
