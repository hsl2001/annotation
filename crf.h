#ifndef CRF_H
#define CRF_H
#include "anno.h"
#include "lbfgs.h"
#include <stdint.h>

#define CRF_L2_C2 1.0
#define CRF_MIN_FREQ 1

/* Linear-chain CRF decoded jointly for both strands with strand-shared emissions.
   0 intergenic; 1-3 CDS(+) codon position; 4-6 intron(+) carrying the codon position of the
   next exon base; 7-9 CDS(-) codon position; 10-12 intron(-). */
static const int boundary_kmer_offsets[] = {-8, -7, -6, -5, -4, -3, -2, -1,
                                           0, 1, 2, 3, 4, 5, 6};
static const int boundary_cwt_offsets[] = {-27, -26, -25, -24, -2, -1, 0, 1, 2};
enum {
  CRF_WINDOW_SIZE = 65536,
  DERIV_COUNT = 3, DERIV_RADIUS = 2,
  STATES = 13, EMIT_STATES = 7, HIDDEN_STATES = 39,
  K2 = 17, K3 = 65, KNUC = 5, EMIT_KMER_LAYERS = 1,
  BOUNDARY_PWM_FROM = -24, BOUNDARY_PWM_COUNT = 33,
  BOUNDARY_KMER_COUNT = sizeof(boundary_kmer_offsets) / sizeof(boundary_kmer_offsets[0]),
  BOUNDARY_CWT_COUNT = sizeof(boundary_cwt_offsets) / sizeof(boundary_cwt_offsets[0]),
  OFF_E3 = 0,
  OFF_SPLICE_START_STOP = OFF_E3 + EMIT_STATES * K3,
  OFF_BP = OFF_SPLICE_START_STOP + 4 * BOUNDARY_KMER_COUNT * K2,
  OFF_BD = OFF_BP + BOUNDARY_CWT_COUNT * WAVE_COUNT,
  OFF_PWM = OFF_BD + DERIV_COUNT * BOUNDARY_CWT_COUNT * WAVE_COUNT,
  WEIGHTS = OFF_PWM + 4 * BOUNDARY_PWM_COUNT * KNUC
};

void crf_init(void);
void crf_train(Contig *contigs, int count, double *weights, int iterations, int threads);
void crf_decode(const double *weights, const Contig *contig, uint8_t *path);
void crf_decode_threads(const double *weights, const Contig *contig, uint8_t *path, int threads);
void crf_save_weights_text(const char *path, const double *weights);

#endif

#ifdef CRF_IMPLEMENTATION

#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const int emit_kmer_sizes[EMIT_KMER_LAYERS] = {3};
static const int emit_kmer_offsets[EMIT_KMER_LAYERS] = {OFF_E3};
static const int emit_kmer_counts[EMIT_KMER_LAYERS] = {K3};

typedef struct {
  uint32_t kmer_emit[EMIT_KMER_LAYERS][2];
  uint8_t boundary_nuc[2][BOUNDARY_PWM_COUNT];
  uint8_t nuc, boundary_k2[2][BOUNDARY_KMER_COUNT];
} Site;

static const int hidden_offsets[STATES] = {0, 1, 4, 8, 10, 11, 14, 18, 20, 24, 29, 33, 38};
static const int hidden_counts[STATES] = {1, 3, 4, 2, 1, 3, 4, 2, 4, 5, 4, 5, 1};
static int trans_a[STATES * STATES], trans_b[STATES * STATES], ntrans;
static int edge_rev[STATES * STATES], edge_kind[STATES * STATES];
static int pair_to_edge[STATES][STATES];
static uint8_t hidden_state[HIDDEN_STATES], hidden_context[HIDDEN_STATES];
typedef struct { uint8_t previous, current; uint16_t transition; } HiddenEdge;
static HiddenEdge hidden_edges[KNUC][HIDDEN_STATES * 4];
static int hidden_edge_count[KNUC];
static uint8_t rev_k2[K2];

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

static void site_context(const Contig *c, int position, Site *site) {
  for (int rev = 0; rev < 2; rev++) {
    for (int layer = 0; layer < EMIT_KMER_LAYERS; layer++) {
      int size = emit_kmer_sizes[layer];
      site->kmer_emit[layer][rev] = (uint32_t)kmer(
        c->seq, c->genome_length, position - (rev ? size / 2 : (size - 1) / 2), size, rev);
    }
    for (int offset = 0; offset < BOUNDARY_PWM_COUNT; offset++) {
      int displacement = BOUNDARY_PWM_FROM + offset;
      int g = rev ? position - displacement - 1 : position + displacement;
      int b = g >= 0 && g < c->genome_length ? base_index(c->seq[g]) : -1;
      site->boundary_nuc[rev][offset] = b < 0 ? 4 : rev ? 3 - b : b;
    }
  }
}

static Site site_at(const Contig *c, int i) {
  int position = c->offset + i;
  int base = base_index(c->seq[position]);
  Site s = {.nuc = base < 0 ? 4 : base};
  for (int offset = 0; offset < BOUNDARY_KMER_COUNT; offset++) {
    int displacement = boundary_kmer_offsets[offset];
    s.boundary_k2[0][offset] = kmer(c->seq, c->genome_length, position + displacement, 2, 0);
    s.boundary_k2[1][offset] = kmer(c->seq, c->genome_length, position - displacement - 2, 2, 1);
  }
  site_context(c, position, &s);
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

static const double derivative_kernels[DERIV_COUNT][2 * DERIV_RADIUS + 1] = {
  {0, -0.5, 0, 0.5, 0},
  {0, 1, -2, 1, 0},
  {-0.5, 1, 0, -1, 0.5}
};

static int boundary_kmer_index(const Site *site, int reverse, int kind, int offset) {
  return OFF_SPLICE_START_STOP + (kind * BOUNDARY_KMER_COUNT + offset) * K2 + site->boundary_k2[reverse][offset];
}

static int boundary_pwm_index(const Site *site, int reverse, int kind, int offset) {
  return OFF_PWM + (kind * BOUNDARY_PWM_COUNT + offset) * KNUC +
         site->boundary_nuc[reverse][offset];
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
      derivative[slot][scale] = reverse && !(slot % 2) ? -value : value;
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
    if (base == 4) return -1;
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
  if (base == 4) return -1;
  if (current == 9)
    return hidden_offsets[9] + (base == 3 ? 1 : base == 1 ? 2 : 0);
  if (current == 8) {
    if (context == 3)
      return base == 3 || base == 1 ? hidden_offsets[8] + 3 : -1;
    if (context == 4) return base == 3 ? hidden_offsets[8] + 3 : -1;
    int next_context = (context == 1 && (base == 3 || base == 1)) || (context == 2 && base == 3) ? 1 :
                       context == 2 && base == 0 ? 2 : 0;
    return hidden_offsets[8] + next_context;
  }
  if (context == 3) return base == 0 ? hidden_offsets[7] : -1;
  if (context == 1 && base == 0) return -1;
  return hidden_offsets[7] + (context == 2 && base == 3);
}

static int initial_allowed(int base, int hidden) {
  return grammar_next(0, hidden_state[hidden], base) == hidden;
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
  for (int base = 0; base < KNUC; base++)
    for (int previous = 0; previous < HIDDEN_STATES; previous++)
      for (int edge = 0; edge < ntrans; edge++) {
        if (trans_a[edge] != hidden_state[previous]) continue;
        int current = grammar_next(previous, trans_b[edge], base);
        if (current < 0) continue;
        if (hidden_edge_count[base] >= HIDDEN_STATES * 4) anno_fail("Too many grammar transitions");
        hidden_edges[base][hidden_edge_count[base]++] = (HiddenEdge){(uint8_t)previous, (uint8_t)current, (uint16_t)edge};
      }
  for (int k = 0; k < 16; k++)
    rev_k2[k] = (uint8_t)(((3 - (k & 3)) << 2) | (3 - ((k >> 2) & 3)));
  rev_k2[16] = 16;
}

static int emit_state(int state) {
  if (state < 7) return state;
  return state <= 9 ? state - 6 : 4 + (state - 9) % 3;
}

static int emission_index(const Site *site, int layer, int state) {
  return emit_kmer_offsets[layer] + emit_state(state) * emit_kmer_counts[layer] +
         (int)site->kmer_emit[layer][state >= 7];
}

static double emission_score(const double *w, const Site *site, int state) {
  double score = 0;
  for (int layer = 0; layer < EMIT_KMER_LAYERS; layer++)
    score += w[emission_index(site, layer, state)];
  return score;
}

static void emission_accumulate(double *values, const Site *site, int state, double value) {
  for (int layer = 0; layer < EMIT_KMER_LAYERS; layer++)
    values[emission_index(site, layer, state)] += value;
}

typedef struct {
  Site *sites;
  double *psi, *log_norm, *log_alpha;
  int capacity;
} Work;

typedef struct {
  float power[WAVE_COUNT];
  double deriv[DERIV_COUNT][WAVE_COUNT];
} CwtCached;

typedef struct {
  const Contig *contig;
  uint8_t nuc_ring[64];
  uint8_t k2_ring[64];
  CwtCached cwt_ring[64];
} SiteStream;

static inline void site_stream_fill_cwt(SiteStream *s, int g) {
  CwtCached *ce = &s->cwt_ring[g & 63];
  if (s->nuc_ring[g & 63] == 4) {
    memset(ce, 0, sizeof(*ce));
    return;
  }
  const float *p = s->contig->cwt_power + (size_t)g * WAVE_COUNT;
  for (int sc = 0; sc < WAVE_COUNT; sc++) ce->power[sc] = p[sc];
  int ok_r1 = s->nuc_ring[(g - 1) & 63] < 4 && s->nuc_ring[(g + 1) & 63] < 4;
  if (!ok_r1) {
    memset(ce->deriv, 0, sizeof(ce->deriv));
    return;
  }
  const float *pm1 = p - WAVE_COUNT, *pp1 = p + WAVE_COUNT;
  for (int sc = 0; sc < WAVE_COUNT; sc++) {
    double v_m1 = pm1[sc], v_0 = p[sc], v_p1 = pp1[sc];
    ce->deriv[0][sc] = -0.5 * v_m1 + 0.5 * v_p1;
    ce->deriv[1][sc] = v_m1 - 2.0 * v_0 + v_p1;
  }
  int ok_r2 = s->nuc_ring[(g - 2) & 63] < 4 && s->nuc_ring[(g + 2) & 63] < 4;
  if (ok_r2) {
    const float *pm2 = p - 2 * WAVE_COUNT, *pp2 = p + 2 * WAVE_COUNT;
    for (int sc = 0; sc < WAVE_COUNT; sc++)
      ce->deriv[2][sc] = -0.5 * pm2[sc] + pm1[sc] - pp1[sc] + 0.5 * pp2[sc];
  } else {
    for (int sc = 0; sc < WAVE_COUNT; sc++) ce->deriv[2][sc] = 0.0;
  }
}

static void site_stream_init(SiteStream *s, const Contig *c, int start_position) {
  s->contig = c;
  int g_start = c->offset + start_position;
  for (int g = g_start - 29; g <= g_start + 27; g++) {
    int b = (g >= 0 && g < c->genome_length) ? base_index(c->seq[g]) : -1;
    s->nuc_ring[g & 63] = (uint8_t)(b < 0 ? 4 : b);
  }
  for (int g = g_start - 8; g <= g_start + 5; g++) {
    uint8_t b0 = s->nuc_ring[g & 63], b1 = s->nuc_ring[(g + 1) & 63];
    s->k2_ring[g & 63] = (b0 < 4 && b1 < 4) ? (uint8_t)((b0 << 2) | b1) : 16;
  }
  for (int g = g_start - 27; g <= g_start + 25; g++)
    site_stream_fill_cwt(s, g);
}

static inline Site site_scores_stream(const double *w, SiteStream *s, int position,
                                      double *emissions, double *psi) {
  const Contig *c = s->contig;
  int g0 = c->offset + position;
  int g_new = g0 + 28;
  int b_new = (g_new >= 0 && g_new < c->genome_length) ? base_index(c->seq[g_new]) : -1;
  s->nuc_ring[g_new & 63] = (uint8_t)(b_new < 0 ? 4 : b_new);

  int g_k2 = g0 + 6;
  uint8_t b0 = s->nuc_ring[g_k2 & 63], b1 = s->nuc_ring[(g_k2 + 1) & 63];
  s->k2_ring[g_k2 & 63] = (b0 < 4 && b1 < 4) ? (uint8_t)((b0 << 2) | b1) : 16;

  site_stream_fill_cwt(s, g0 + 26);

  Site site;
  site.nuc = s->nuc_ring[g0 & 63];
  site_context(c, g0, &site);
  for (int offset = 0; offset < BOUNDARY_KMER_COUNT; offset++) {
    uint8_t k2f = s->k2_ring[(g0 + boundary_kmer_offsets[offset]) & 63];
    site.boundary_k2[0][offset] = k2f;
    site.boundary_k2[1][BOUNDARY_KMER_COUNT - 1 - offset] = rev_k2[k2f];
  }

  for (int state = 0; state < STATES; state++) {
    double score = emission_score(w, &site, state);
    if (!isfinite(score))
      anno_fail("Non-finite CRF emission at position %d", position + 1);
    emissions[state] = score;
  }

  double b_score[2][4];
  for (int rev = 0; rev < 2; rev++) {
    const uint8_t *k2 = site.boundary_k2[rev];
    for (int kind = 0; kind < 4; kind++) {
      const double *wk = w + OFF_SPLICE_START_STOP + kind * (BOUNDARY_KMER_COUNT * K2);
      double score = 0;
      for (int offset = 0; offset < BOUNDARY_KMER_COUNT; offset++)
        score += wk[offset * K2 + k2[offset]];
      for (int offset = 0; offset < BOUNDARY_PWM_COUNT; offset++)
        score += w[boundary_pwm_index(&site, rev, kind, offset)];
      b_score[rev][kind] = score;
    }
    double cw = 0;
    for (int offset = 0; offset < BOUNDARY_CWT_COUNT; offset++) {
      int displacement = boundary_cwt_offsets[offset];
      int g = rev ? g0 - displacement - 1 : g0 + displacement;
      const CwtCached *ce = &s->cwt_ring[g & 63];
      const double *wbp = w + OFF_BP + offset * WAVE_COUNT;
      const double *wd0 = w + OFF_BD + (0 * BOUNDARY_CWT_COUNT + offset) * WAVE_COUNT;
      const double *wd1 = w + OFF_BD + (1 * BOUNDARY_CWT_COUNT + offset) * WAVE_COUNT;
      const double *wd2 = w + OFF_BD + (2 * BOUNDARY_CWT_COUNT + offset) * WAVE_COUNT;
      for (int scale = 0; scale < WAVE_COUNT; scale++) {
        cw += wbp[scale] * ce->power[scale];
        double d0 = ce->deriv[0][scale];
        double d1 = ce->deriv[1][scale];
        double d2 = ce->deriv[2][scale];
        if (rev) { d0 = -d0; d2 = -d2; }
        cw += wd0[scale] * d0;
        cw += wd1[scale] * d1;
        cw += wd2[scale] * d2;
      }
    }
    b_score[rev][1] += cw;
  }

  for (int edge = 0; edge < ntrans; edge++) {
    int reverse = edge_rev[edge], kind = edge_kind[edge];
    double score = kind >= 0 ? b_score[reverse][kind] : 0;
    double p = score + emissions[trans_b[edge]];
    if (!isfinite(p))
      anno_fail("Non-finite CRF transition at position %d", position + 1);
    psi[edge] = p;
  }
  return site;
}

static int gold_path_valid(const Contig *c) {
  int hidden = 0;
  for (int position = 0; position < c->length; position++) {
    int base = base_index(c->seq[c->offset + position]);
    hidden = grammar_next(hidden, c->label[c->offset + position], base < 0 ? 4 : base);
    if (hidden < 0) return 0;
  }
  return terminal_allowed(hidden);
}

static void work_ensure(Work *k, int length) {
  if (length <= k->capacity) return;
  free(k->sites), free(k->psi), free(k->log_norm), free(k->log_alpha);
  k->sites = anno_alloc((size_t)length, sizeof(Site));
  k->psi = anno_alloc((size_t)length, (size_t)ntrans * sizeof(double));
  k->log_norm = anno_alloc((size_t)length, sizeof(double));
  k->log_alpha = anno_alloc((size_t)length, HIDDEN_STATES * sizeof(double));
  k->capacity = length;
}

static void work_free(Work *k) {
  free(k->sites), free(k->psi), free(k->log_norm), free(k->log_alpha);
}

static inline int forward_edge_allowed(const double *prev_log_alpha, const HiddenEdge *step,
                                       int last) {
  return prev_log_alpha[step->previous] > -INFINITY &&
         (!last || terminal_allowed(step->current));
}

/* Normalizes unnormalized log-forward scores in place and returns their log-sum-exp. */
static double normalize_log_alpha(double *log_alpha, int position) {
  double max_score = -INFINITY;
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
    if (log_alpha[hidden] > max_score) max_score = log_alpha[hidden];
  if (!isfinite(max_score)) anno_fail("Invalid CRF partition at position %d", position + 1);
  double sum = 0;
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
    if (log_alpha[hidden] > -INFINITY) sum += exp(log_alpha[hidden] - max_score);
  double log_norm = max_score + log(sum);
  if (!isfinite(log_norm)) anno_fail("Invalid CRF partition at position %d", position + 1);
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) log_alpha[hidden] -= log_norm;
  return log_norm;
}

/*
 * Accumulates the exact observed-minus-expected gradient for one training window.
 * Forward scores are kept in log space so no reachable state underflows to zero, and the
 * backward pass uses P(prev, cur) = P(cur) * exp(log_alpha[prev] + psi - log_norm - log_alpha[cur]),
 * whose exponent is never positive, so posteriors stay in [0, 1] for arbitrarily sharp weights.
 */
static double contig_gradient(const double *w, const Contig *c,
                              Work *k, double *gradient) {
  work_ensure(k, c->length);
  SiteStream stream;
  site_stream_init(&stream, c, 0);
  double logz = 0, gold_score = 0;
  for (int position = 0; position < c->length; position++) {
    double emissions[STATES];
    double *psi = k->psi + (size_t)position * ntrans;
    Site site = site_scores_stream(w, &stream, position, emissions, psi);
    k->sites[position] = site;
    int gold_state = c->label[c->offset + position];
    gold_score += emissions[gold_state];

    double *L = k->log_alpha + (size_t)position * HIDDEN_STATES;
    if (!position) {
      for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
        L[hidden] = initial_allowed(site.nuc, hidden) &&
                    (c->length > 1 || terminal_allowed(hidden))
          ? emissions[hidden_state[hidden]] : -INFINITY;
    } else {
      int prev_gold = c->label[c->offset + position - 1];
      int gold_edge = pair_to_edge[prev_gold][gold_state];
      if (gold_edge < 0) anno_fail("Invalid transition in validated CRF training path");
      gold_score += psi[gold_edge] - emissions[gold_state];

      const double *prev_L = L - HIDDEN_STATES;
      const HiddenEdge *edges = hidden_edges[site.nuc];
      int edge_count = hidden_edge_count[site.nuc];
      int last = position == c->length - 1;
      double incoming_sum[HIDDEN_STATES] = {0};
      for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) L[hidden] = -INFINITY;
      for (int edge = 0; edge < edge_count; edge++) {
        const HiddenEdge *step = edges + edge;
        if (!forward_edge_allowed(prev_L, step, last)) continue;
        double score = prev_L[step->previous] + psi[step->transition];
        if (score > L[step->current]) L[step->current] = score;
      }
      for (int edge = 0; edge < edge_count; edge++) {
        const HiddenEdge *step = edges + edge;
        if (!forward_edge_allowed(prev_L, step, last)) continue;
        incoming_sum[step->current] +=
          exp(prev_L[step->previous] + psi[step->transition] - L[step->current]);
      }
      for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
        if (L[hidden] > -INFINITY) L[hidden] += log(incoming_sum[hidden]);
    }
    k->log_norm[position] = normalize_log_alpha(L, position);
    logz += k->log_norm[position];
  }

  double likelihood = gold_score - logz;
  double B[HIDDEN_STATES];
  const double *last_L = k->log_alpha + (size_t)(c->length - 1) * HIDDEN_STATES;
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
    B[hidden] = exp(last_L[hidden]);

  for (int position = c->length - 1; position >= 0; position--) {
    const Site *site = k->sites + position;
    const double *L = k->log_alpha + (size_t)position * HIDDEN_STATES;
    int state = c->label[c->offset + position];

    double state_marginals[STATES] = {0};
    for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
      state_marginals[hidden_state[hidden]] += B[hidden];
    state_marginals[state] -= 1.0;

    for (int st = 0; st < STATES; st++)
      emission_accumulate(gradient, site, st, -state_marginals[st]);

    if (position) {
      const double *prev_L = L - HIDDEN_STATES;
      const double *psi = k->psi + (size_t)position * ntrans;
      double log_norm = k->log_norm[position];
      int last = position == c->length - 1;
      int previous_state = c->label[c->offset + position - 1];
      int gold_edge = pair_to_edge[previous_state][state];

      double edge_marginals[STATES * STATES] = {0};
      double B_prev[HIDDEN_STATES] = {0};
      const HiddenEdge *edges = hidden_edges[site->nuc];
      int edge_count = hidden_edge_count[site->nuc];
      for (int edge = 0; edge < edge_count; edge++) {
        const HiddenEdge *step = edges + edge;
        if (!(B[step->current] > 0.0) || !forward_edge_allowed(prev_L, step, last)) continue;
        double edge_prob = B[step->current] *
          exp(prev_L[step->previous] + psi[step->transition] - log_norm - L[step->current]);
        B_prev[step->previous] += edge_prob;
        edge_marginals[step->transition] += edge_prob;
      }
      double sum_b = 0;
      for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) sum_b += B_prev[hidden];
      if (!(sum_b > 0.0) || !isfinite(sum_b))
        anno_fail("Invalid CRF posterior at position %d", position + 1);
      for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) B[hidden] = B_prev[hidden] / sum_b;

      double b_coeff[4][2] = {{0}};
      int has_boundary = 0;
      edge_marginals[gold_edge] -= 1.0;

      for (int edge = 0; edge < ntrans; edge++) {
        double m = edge_marginals[edge];
        if (m != 0) {
          int reverse = edge_rev[edge], kind = edge_kind[edge];
          if (kind >= 0) {
            b_coeff[kind][reverse] -= m;
            has_boundary = 1;
          }
        }
      }
      if (has_boundary) {
        int genome_position = c->offset + position;
        for (int rev = 0; rev < 2; rev++) {
          for (int kind = 0; kind < 4; kind++) {
            double val = b_coeff[kind][rev];
            if (val == 0) continue;
            for (int offset = 0; offset < BOUNDARY_KMER_COUNT; offset++)
              gradient[boundary_kmer_index(site, rev, kind, offset)] += val;
            for (int offset = 0; offset < BOUNDARY_PWM_COUNT; offset++)
              gradient[boundary_pwm_index(site, rev, kind, offset)] += val;
          }
          double val = b_coeff[1][rev];
          if (val == 0) continue;
          for (int offset = 0; offset < BOUNDARY_CWT_COUNT; offset++) {
            float power[WAVE_COUNT];
            double derivative[DERIV_COUNT][WAVE_COUNT];
            boundary_cwt_at(c, genome_position, rev, offset, power, derivative);
            double *gbp = gradient + OFF_BP + offset * WAVE_COUNT;
            for (int scale = 0; scale < WAVE_COUNT; scale++) {
              gbp[scale] += val * power[scale];
              for (int slot = 0; slot < DERIV_COUNT; slot++)
                gradient[derivative_index(slot, offset, scale)] +=
                  val * derivative[slot][scale];
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
} CrfWindow;

static int next_window_end(const Contig *contig, int start, int limit) {
  int end = limit - start > CRF_WINDOW_SIZE ? start + CRF_WINDOW_SIZE : limit;
  if (end == limit || !contig->label[end - 1] || !contig->label[end]) return end;
  int gene_start = end - 1, gene_end = end;
  while (gene_start > start && contig->label[gene_start - 1]) gene_start--;
  while (gene_end < limit && contig->label[gene_end]) gene_end++;
  return gene_start > start && end - gene_start <= gene_end - end ? gene_start : gene_end;
}

static size_t make_windows(Contig *contigs, int count, CrfWindow **windows, size_t *skipped) {
  size_t window_count = 0, capacity = 0;
  *skipped = 0;
  for (int i = 0; i < count; i++) {
    Contig *contig = contigs + i;
    int limit = contig->offset + contig->length;
    for (int start = contig->offset; start < limit;) {
      int end = next_window_end(contig, start, limit);
      if (end <= start) anno_fail("Cannot split contig into CRF windows: %s", contig->name);
      Contig instance = *contig;
      instance.offset = start;
      instance.length = end - start;
      if (!gold_path_valid(&instance)) {
        (*skipped)++;
        start = end;
        continue;
      }
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
  size_t count, skipped;
  int threads;
} CrfTrain;

typedef struct {
  CrfTrain *train;
  int thread_id;
  const double *w;
  double likelihood;
  size_t positions;
  double *gradient;
  Work work;
} CrfWorker;

static void *crf_worker(void *arg) {
  CrfWorker *worker = arg;
  CrfTrain *train = worker->train;
  worker->likelihood = 0;
  worker->positions = 0;
  memset(worker->gradient, 0, WEIGHTS * sizeof(*worker->gradient));
  for (size_t slot = (size_t)worker->thread_id; slot < train->count;
       slot += (size_t)train->threads) {
    const CrfWindow *window = train->windows + slot;
    Contig instance = *window->contig;
    instance.offset = window->start;
    instance.length = window->length;
    double score = contig_gradient(worker->w, &instance, &worker->work, worker->gradient);
    if (!isfinite(score)) anno_fail("Non-finite training likelihood");
    worker->likelihood += score;
    worker->positions += (size_t)instance.length;
  }
  return NULL;
}

typedef struct {
  CrfTrain train;
  CrfWorker *workers;
  pthread_t *threads;
  double likelihood;
  size_t positions;
  double *scale, *weights;
} CrfObjective;

static double crf_objective(const double *w, double *gradient, void *context) {
  CrfObjective *objective = context;
  CrfTrain *train = &objective->train;
  double likelihood = 0;
  size_t positions = 0;

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
    for (int f = 0; f < WEIGHTS; f++) gradient[f] += objective->workers[t].gradient[f];
  }
  if (!isfinite(likelihood)) anno_fail("Non-finite training likelihood");

  objective->likelihood = likelihood;
  objective->positions = positions;
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
          iteration, evaluation->train.count, evaluation->train.skipped,
          evaluation->likelihood / evaluation->positions, objective, gradient_norm,
          step, line_search_trials);
}

/* L-BFGS optimizes v with w = scale * v, so frequent dense and rare sparse features share curvature. */
static double crf_scaled_objective(const double *v, double *gradient, void *context) {
  CrfObjective *objective = context;
  for (int f = 0; f < WEIGHTS; f++)
    objective->weights[f] = objective->scale[f] * v[f];
  double loss = crf_objective(objective->weights, gradient, context);
  for (int f = 0; f < WEIGHTS; f++)
    gradient[f] *= objective->scale[f];
  return loss;
}

static inline void feature_observe(double *scale, double *freq, int f, double value) {
  scale[f] += value * value;
  if (freq && value != 0) freq[f] += 1;
}

static void feature_counts(const CrfWindow *windows, size_t count, double *scale, double *freq) {
  memset(scale, 0, WEIGHTS * sizeof(*scale));
  if (freq) memset(freq, 0, WEIGHTS * sizeof(*freq));
  for (size_t w = 0; w < count; w++) {
    Contig c = *windows[w].contig;
    c.offset = windows[w].start;
    c.length = windows[w].length;
    for (int i = 0; i < c.length; i++) {
      Site site = site_at(&c, i);
      int state = c.label[c.offset + i];
      for (int layer = 0; layer < EMIT_KMER_LAYERS; layer++)
        feature_observe(scale, freq, emission_index(&site, layer, state), 1);
      if (!i) continue;
      int previous = c.label[c.offset + i - 1];
      int reverse = previous >= 7 || state >= 7, kind = boundary_kind(previous, state);
      if (kind >= 0) {
        for (int offset = 0; offset < BOUNDARY_KMER_COUNT; offset++)
          feature_observe(scale, freq, boundary_kmer_index(&site, reverse, kind, offset), 1);
        for (int offset = 0; offset < BOUNDARY_PWM_COUNT; offset++)
          feature_observe(scale, freq, boundary_pwm_index(&site, reverse, kind, offset), 1);
        if (kind == 1) {
          int genome_position = c.offset + i;
          for (int offset = 0; offset < BOUNDARY_CWT_COUNT; offset++) {
            float power[WAVE_COUNT];
            double derivative[DERIV_COUNT][WAVE_COUNT];
            boundary_cwt_at(&c, genome_position, reverse, offset, power, derivative);
            int base = offset * WAVE_COUNT;
            for (int k = 0; k < WAVE_COUNT; k++) {
              feature_observe(scale, freq, OFF_BP + base + k, power[k]);
              for (int slot = 0; slot < DERIV_COUNT; slot++)
                feature_observe(scale, freq, derivative_index(slot, offset, k),
                                derivative[slot][k]);
            }
          }
        }
      }
    }
  }
}

static void counts_to_scales(double *scale) {
  for (int f = 0; f < WEIGHTS; f++) scale[f] = 1 / sqrt(fmax(scale[f], 1));
}

static int feature_active(double freq) {
  return freq >= CRF_MIN_FREQ;
}

static const struct {
  const char *name;
  int from, to;
} feature_groups[] = {
  {"emit_k3", OFF_E3, OFF_SPLICE_START_STOP},
  {"boundary_k2", OFF_SPLICE_START_STOP, OFF_BP},
  {"boundary_cwt_power", OFF_BP, OFF_BD},
  {"boundary_cwt_deriv", OFF_BD, OFF_PWM},
  {"boundary_pwm", OFF_PWM, WEIGHTS}
};
enum { FEATURE_GROUPS = sizeof(feature_groups) / sizeof(feature_groups[0]) };

typedef struct {
  CrfObjective *inner;
  int *index, count;
  double *v, *gradient;
} CrfActive;

static double crf_active_objective(const double *x, double *gradient, void *context) {
  CrfActive *active = context;
  for (int i = 0; i < active->count; i++) active->v[active->index[i]] = x[i];
  double loss = crf_scaled_objective(active->v, active->gradient, active->inner);
  for (int i = 0; i < active->count; i++) gradient[i] = active->gradient[active->index[i]];
  return loss;
}

static void crf_active_progress(void *context, int iteration, double objective, double step,
                                int line_search_trials, size_t n, const double *gradient) {
  crf_progress(((CrfActive *)context)->inner, iteration, objective, step, line_search_trials,
               n, gradient);
}

void crf_train(Contig *contigs, int count, double *w, int iterations, int threads) {
  if (count < 1 || iterations < 1)
    anno_fail("CRF training requires instances and positive iterations");
  if (threads < 1) anno_fail("CRF training requires at least one worker thread");
  CrfWindow *windows = NULL;
  size_t skipped = 0;
  size_t window_count = make_windows(contigs, count, &windows, &skipped);
  if (!window_count) anno_fail("No valid CRF training windows (%zu windows skipped)", skipped);
  CrfObjective objective = {
    .train = {.windows = windows, .count = window_count, .skipped = skipped, .threads = threads}
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
  double *freq = anno_alloc(WEIGHTS, sizeof(*freq));
  feature_counts(windows, window_count, objective.scale, freq);
  CrfActive active = {
    .inner = &objective,
    .index = anno_alloc(WEIGHTS, sizeof(*active.index)),
    .v = anno_alloc(WEIGHTS, sizeof(*active.v)),
    .gradient = anno_alloc(WEIGHTS, sizeof(*active.gradient))
  };
  fprintf(stderr, "Active features (min frequency %d):", CRF_MIN_FREQ);
  for (int g = 0; g < FEATURE_GROUPS; g++) {
    int group_active = 0;
    for (int f = feature_groups[g].from; f < feature_groups[g].to; f++)
      if (feature_active(freq[f])) active.index[active.count++] = f, group_active++;
    fprintf(stderr, " %s %d/%d", feature_groups[g].name, group_active,
            feature_groups[g].to - feature_groups[g].from);
  }
  fprintf(stderr, ", total %d/%d\n", active.count, WEIGHTS);
  free(freq);
  counts_to_scales(objective.scale);
  double *x = anno_alloc((size_t)active.count, sizeof(*x));
  for (int i = 0; i < active.count; i++)
    x[i] = w[active.index[i]] / objective.scale[active.index[i]];
  LbfgsStatus status = lbfgs_minimize((size_t)active.count, x, crf_active_objective,
                                      crf_active_progress, &active, iterations, &completed);
  memset(w, 0, WEIGHTS * sizeof(*w));
  for (int i = 0; i < active.count; i++)
    w[active.index[i]] = x[i] * objective.scale[active.index[i]];
  free(x), free(active.index), free(active.v), free(active.gradient);
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

static inline void viterbi_relax(uint8_t nuc, int position, const double *emissions,
                                 const double *transitions, const double *previous,
                                 double *current, uint8_t *back) {
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) {
    current[hidden] = position == 0 && initial_allowed(nuc, hidden) ? emissions[hidden_state[hidden]] : -INFINITY;
    back[hidden] = 0;
  }
  if (position) {
    const HiddenEdge *edges = hidden_edges[nuc];
    int edge_count = hidden_edge_count[nuc];
    for (int edge = 0; edge < edge_count; edge++) {
      const HiddenEdge *step = edges + edge;
      double prev_score = previous[step->previous];
      if (prev_score == -INFINITY) continue;
      double score = prev_score + transitions[step->transition];
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

typedef struct {
  const double *weights;
  const Contig *contig;
  uint8_t *back;
  int chunk_id;
  int c_start, c_end;
  double start_vec[HIDDEN_STATES];
  double end_vec[HIDDEN_STATES];
} ViterbiChunk;

static void *viterbi_chunk_worker(void *argument) {
  ViterbiChunk *chunk = argument;
  const double *w = chunk->weights;
  const Contig *c = chunk->contig;
  double previous[HIDDEN_STATES], current[HIDDEN_STATES];
  double emissions[STATES], transitions[STATES * STATES];
  SiteStream stream;

  if (chunk->chunk_id == 0) {
    memset(previous, 0, sizeof(previous));
    site_stream_init(&stream, c, 0);
    for (int pos = 0; pos < chunk->c_end; pos++) {
      Site site = site_scores_stream(w, &stream, pos, emissions, transitions);
      viterbi_relax(site.nuc, pos, emissions, transitions, previous, current,
                    chunk->back + (size_t)pos * HIDDEN_STATES);
      memcpy(previous, current, sizeof(previous));
    }
    memcpy(chunk->end_vec, previous, sizeof(previous));
    return NULL;
  }

  int w_start = chunk->c_start - CRF_WINDOW_SIZE;
  if (w_start < 0) w_start = 0;
  site_stream_init(&stream, c, w_start);
  if (w_start == 0) {
    memset(previous, 0, sizeof(previous));
  } else {
    previous[0] = 0.0;
    for (int h = 1; h < HIDDEN_STATES; h++) previous[h] = -INFINITY;
  }
  uint8_t dummy_back[HIDDEN_STATES];
  for (int pos = w_start; pos < chunk->c_start; pos++) {
    Site site = site_scores_stream(w, &stream, pos, emissions, transitions);
    viterbi_relax(site.nuc, pos, emissions, transitions, previous, current, dummy_back);
    memcpy(previous, current, sizeof(previous));
  }
  memcpy(chunk->start_vec, previous, sizeof(previous));
  for (int pos = chunk->c_start; pos < chunk->c_end; pos++) {
    Site site = site_scores_stream(w, &stream, pos, emissions, transitions);
    viterbi_relax(site.nuc, pos, emissions, transitions, previous, current,
                  chunk->back + (size_t)pos * HIDDEN_STATES);
    memcpy(previous, current, sizeof(previous));
  }
  memcpy(chunk->end_vec, previous, sizeof(previous));
  return NULL;
}

/* Exact whole-contig MAP decoding (optionally multi-threaded across chunks). */
void crf_decode_threads(const double *w, const Contig *c, uint8_t *path, int threads) {
  int length = c->length;
  uint8_t *back = anno_alloc((size_t)length, HIDDEN_STATES * sizeof(*back));
  int max_chunks = length / CRF_WINDOW_SIZE;
  int n_chunks = threads > 1 ? (threads < max_chunks ? threads : max_chunks) : 1;
  if (n_chunks < 1) n_chunks = 1;

  double previous[HIDDEN_STATES] = {0}, current[HIDDEN_STATES];
  if (n_chunks == 1) {
    SiteStream stream;
    site_stream_init(&stream, c, 0);
    double emissions[STATES], transitions[STATES * STATES];
    for (int position = 0; position < length; position++) {
      Site site = site_scores_stream(w, &stream, position, emissions, transitions);
      viterbi_relax(site.nuc, position, emissions, transitions, previous, current,
                    back + (size_t)position * HIDDEN_STATES);
      memcpy(previous, current, sizeof(previous));
    }
  } else {
    ViterbiChunk *chunks = anno_alloc((size_t)n_chunks, sizeof(*chunks));
    pthread_t *tids = anno_alloc((size_t)n_chunks, sizeof(*tids));
    for (int k = 0; k < n_chunks; k++) {
      chunks[k] = (ViterbiChunk){
        .weights = w,
        .contig = c,
        .back = back,
        .chunk_id = k,
        .c_start = (int)((long long)k * length / n_chunks),
        .c_end = (int)((long long)(k + 1) * length / n_chunks)
      };
      if (pthread_create(tids + k, NULL, viterbi_chunk_worker, chunks + k))
        anno_fail("pthread_create failed");
    }
    for (int k = 0; k < n_chunks; k++)
      if (pthread_join(tids[k], NULL)) anno_fail("pthread_join failed");

    double emissions[STATES], transitions[STATES * STATES];
    uint8_t dummy_back[HIDDEN_STATES];
    for (int k = 1; k < n_chunks; k++) {
      double true_prev[HIDDEN_STATES], approx_prev[HIDDEN_STATES];
      double true_cur[HIDDEN_STATES], approx_cur[HIDDEN_STATES];
      memcpy(true_prev, chunks[k - 1].end_vec, sizeof(true_prev));
      memcpy(approx_prev, chunks[k].start_vec, sizeof(approx_prev));
      if (memcmp(true_prev, approx_prev, sizeof(true_prev)) == 0) continue;

      SiteStream stream;
      site_stream_init(&stream, c, chunks[k].c_start);
      int coalesced = 0;
      for (int pos = chunks[k].c_start; pos < chunks[k].c_end; pos++) {
        Site site = site_scores_stream(w, &stream, pos, emissions, transitions);
        viterbi_relax(site.nuc, pos, emissions, transitions, approx_prev, approx_cur, dummy_back);
        viterbi_relax(site.nuc, pos, emissions, transitions, true_prev, true_cur,
                      back + (size_t)pos * HIDDEN_STATES);
        memcpy(approx_prev, approx_cur, sizeof(approx_prev));
        memcpy(true_prev, true_cur, sizeof(true_prev));
        if (memcmp(true_prev, approx_prev, sizeof(true_prev)) == 0) {
          coalesced = 1;
          break;
        }
      }
      if (!coalesced)
        memcpy(chunks[k].end_vec, true_prev, sizeof(true_prev));
    }
    memcpy(previous, chunks[n_chunks - 1].end_vec, sizeof(previous));
    free(chunks), free(tids);
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

void crf_decode(const double *w, const Contig *c, uint8_t *path) {
  crf_decode_threads(w, c, path, 1);
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
  static const char *const state_names[EMIT_STATES] = {
    "intergenic",
    "CDS:0", "CDS:1", "CDS:2",
    "intron:0", "intron:1", "intron:2"
  };
  static const char *const kind_names[4] = {
    "5ss", "3ss", "start_codon", "stop_codon"
  };

  FILE *file = fopen(path, "w");
  if (!file) anno_fail("Cannot write weights text: %s", path);
  fprintf(file, "# index\tgroup\ttarget\tfeature\tweight\n");

  char kbuf[4];
  for (int layer = 0; layer < EMIT_KMER_LAYERS; layer++) {
    int size = emit_kmer_sizes[layer];
    for (int s = 0; s < EMIT_STATES; s++)
      for (int k = 0; k < emit_kmer_counts[layer]; k++) {
        int idx = emit_kmer_offsets[layer] + s * emit_kmer_counts[layer] + k;
        if (!(float)w[idx]) continue;
        format_kmer(k, size, kbuf);
        fprintf(file, "%d\temit_k%d\t%s\tk%d=%s\t%.9g\n",
                idx, size, state_names[s], size, kbuf, (float)w[idx]);
      }
  }
  for (int kind = 0; kind < 4; kind++)
    for (int offset = 0; offset < BOUNDARY_KMER_COUNT; offset++)
      for (int k = 0; k < K2; k++) {
        int idx = OFF_SPLICE_START_STOP + (kind * BOUNDARY_KMER_COUNT + offset) * K2 + k;
        if (!(float)w[idx]) continue;
        format_kmer(k, 2, kbuf);
        fprintf(file, "%d\tboundary_k2\t%s\toffset=%+d,k2=%s\t%.9g\n",
                idx, kind_names[kind], boundary_kmer_offsets[offset], kbuf, (float)w[idx]);
      }
  for (int offset = 0; offset < BOUNDARY_CWT_COUNT; offset++)
    for (int sc = 0; sc < WAVE_COUNT; sc++) {
      int idx = OFF_BP + offset * WAVE_COUNT + sc;
      if (!(float)w[idx]) continue;
      fprintf(file, "%d\tboundary_cwt_power\t%s\toffset=%+d,scale=%d\t%.9g\n",
              idx, kind_names[1], boundary_cwt_offsets[offset], wave_sizes[sc], (float)w[idx]);
    }
  for (int slot = 0; slot < DERIV_COUNT; slot++)
    for (int offset = 0; offset < BOUNDARY_CWT_COUNT; offset++)
      for (int sc = 0; sc < WAVE_COUNT; sc++) {
        int idx = derivative_index(slot, offset, sc);
        if (!(float)w[idx]) continue;
        fprintf(file, "%d\tboundary_cwt_deriv%d\t%s\toffset=%+d,scale=%d\t%.9g\n",
                idx, slot + 1, kind_names[1], boundary_cwt_offsets[offset], wave_sizes[sc],
                (float)w[idx]);
      }
  for (int kind = 0; kind < 4; kind++)
    for (int offset = 0; offset < BOUNDARY_PWM_COUNT; offset++)
      for (int b = 0; b < KNUC; b++) {
        int idx = OFF_PWM + (kind * BOUNDARY_PWM_COUNT + offset) * KNUC + b;
        if (!(float)w[idx]) continue;
        fprintf(file, "%d\tboundary_pwm\t%s\toffset=%+d,base=%c\t%.9g\n",
                idx, kind_names[kind], BOUNDARY_PWM_FROM + offset, "ACGTN"[b], (float)w[idx]);
      }

  if (fclose(file)) anno_fail("Cannot finish writing weights text: %s", path);
}

#endif
