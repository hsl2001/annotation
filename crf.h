#ifndef CRF_H
#define CRF_H

#include "anno.h"
#include "lbfgs.h"
#include <stdint.h>

#ifndef CRF_L2_C2
#define CRF_L2_C2 2.0
#endif
#ifndef CRF_WINDOW_SIZE
#define CRF_WINDOW_SIZE 65536
#endif
#ifndef CRF_DECODE_CHUNK
#define CRF_DECODE_CHUNK 262144
#endif
#ifndef CRF_GPU_GANGS
#define CRF_GPU_GANGS 26
#endif
#ifndef CRF_GPU_VEC
#define CRF_GPU_VEC 32
#endif
#ifndef CRF_GPU_THREADS
#define CRF_GPU_THREADS (CRF_GPU_GANGS * CRF_GPU_VEC)
#endif
#ifndef CRF_GPU_BATCH_POS
#define CRF_GPU_BATCH_POS 3200000
#endif

/* Linear-chain CRF decoded jointly for both strands.
   0 intergenic; 1-3 CDS(+) codon position; 4-6 intron(+) carrying the codon position of the
   next exon base; 7-9 CDS(-) codon position; 10-12 intron(-). */
enum { STATES = 13, HIDDEN_STATES = 47 };
enum { K2 = 17, K4 = 257, KNUC = 5, RADIUS = 6,
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
#include <stdatomic.h>
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
static HiddenEdge fwd_edges[KNUC][HIDDEN_STATES * 4];
static int fwd_offset[KNUC][HIDDEN_STATES + 1];
static HiddenEdge bwd_edges[KNUC][HIDDEN_STATES * 4];
static int bwd_offset[KNUC][HIDDEN_STATES + 1];
static HiddenEdge tr_edges[KNUC][HIDDEN_STATES * 4];
static int tr_offset[KNUC][STATES * STATES + 1];
static uint8_t initial_states[KNUC][HIDDEN_STATES];
static uint8_t full_init_states[KNUC][HIDDEN_STATES];
static uint8_t term_states[HIDDEN_STATES];

#ifdef _OPENACC
#pragma acc routine seq
#endif
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

#ifdef _OPENACC
#pragma acc routine seq
#endif
static Site site_at_raw(const char *seq, const float *cwt_power, int genome_length, int position) {
  int base = base_index(seq[position]);
  Site s = {
    .k4_emit = {kmer(seq, genome_length, position - 2, 4, 0),
                kmer(seq, genome_length, position - 1, 4, 1)},
    .k4_trans = {kmer(seq, genome_length, position - 2, 4, 0),
                 kmer(seq, genome_length, position - 2, 4, 1)},
    .nuc = base < 0 ? 4 : base
  };
  for (int offset = 0; offset < FLANKS; offset++) {
    int displacement = -RADIUS + 2 * offset;
    s.flank[0][offset] = kmer(seq, genome_length, position + displacement, 2, 0);
    s.flank[1][offset] = kmer(seq, genome_length, position - displacement - 2, 2, 1);
  }
  for (int scale = 0; scale < WAVE_COUNT; scale++)
    s.cwt_power[scale] = cwt_power[(size_t)position * WAVE_COUNT + scale];
  if (position > 0 && position + 1 < genome_length &&
      base_index(seq[position - 1]) >= 0 && base >= 0 &&
      base_index(seq[position + 1]) >= 0) {
    for (int scale = 0; scale < WAVE_COUNT; scale++) {
      double left = cwt_power[(size_t)(position - 1) * WAVE_COUNT + scale];
      double right = cwt_power[(size_t)(position + 1) * WAVE_COUNT + scale];
      double slope = (right - left) / (right + left + 1e-12);
      s.derivative[scale] = (int16_t)lrint(slope * INT16_MAX);
    }
  }
  return s;
}

static Site site_at(const Contig *c, int i) {
  return site_at_raw(c->seq, c->cwt_power, c->genome_length, c->offset + i);
}

#ifdef _OPENACC
#pragma acc routine seq
#endif
static int boundary_kind(int previous, int current) {
  if ((previous >= 1 && previous <= 3 && current >= 4 && current <= 6) ||
      (previous >= 10 && current >= 7 && current <= 9)) return 0;
  if ((previous >= 4 && previous <= 6 && current >= 1 && current <= 3) ||
      (previous >= 7 && previous <= 9 && current >= 10)) return 1;
  if ((previous == 0 && current == 1) || (previous == 7 && current == 0)) return 2;
  if ((previous == 3 && current == 0) || (previous == 0 && current == 9)) return 3;
  return -1;
}

#ifdef _OPENACC
#pragma acc routine seq
#endif
static void boundary_cwt_raw(const char *seq, const float *cwt_power, int genome_length,
                             int position, int reverse, int offset,
                             float *power, double *derivative) {
  int displacement = -RADIUS + 2 * offset;
  int g = reverse ? position - displacement - 2 : position + displacement;
  if (g < 0 || g >= genome_length || base_index(seq[g]) < 0) {
    for (int scale = 0; scale < WAVE_COUNT; scale++) power[scale] = 0, derivative[scale] = 0;
    return;
  }
  for (int scale = 0; scale < WAVE_COUNT; scale++)
    power[scale] = cwt_power[(size_t)g * WAVE_COUNT + scale];
  if (g > 0 && g + 1 < genome_length &&
      base_index(seq[g - 1]) >= 0 && base_index(seq[g + 1]) >= 0) {
    for (int scale = 0; scale < WAVE_COUNT; scale++) {
      double left = cwt_power[(size_t)(g - 1) * WAVE_COUNT + scale];
      double right = cwt_power[(size_t)(g + 1) * WAVE_COUNT + scale];
      double slope = (right - left) / (right + left + 1e-12);
      derivative[scale] = reverse ? -slope : slope;
    }
  } else {
    for (int scale = 0; scale < WAVE_COUNT; scale++) derivative[scale] = 0;
  }
}

static void boundary_cwt_at(const Contig *c, int position, int reverse, int offset,
                            float *power, double *derivative) {
  boundary_cwt_raw(c->seq, c->cwt_power, c->genome_length, position, reverse, offset,
                   power, derivative);
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
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
    term_states[hidden] = (uint8_t)terminal_allowed(hidden);
  for (int base = 0; base < KNUC; base++) {
    for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
      full_init_states[base][hidden] = (uint8_t)initial_allowed(base, hidden, 1);
    int fc = 0;
    for (int curr = 0; curr < HIDDEN_STATES; curr++) {
      fwd_offset[base][curr] = fc;
      for (int e = 0; e < hidden_edge_count[base]; e++)
        if (hidden_edges[base][e].current == curr)
          fwd_edges[base][fc++] = hidden_edges[base][e];
    }
    fwd_offset[base][HIDDEN_STATES] = fc;

    int bc = 0;
    for (int prev = 0; prev < HIDDEN_STATES; prev++) {
      bwd_offset[base][prev] = bc;
      for (int e = 0; e < hidden_edge_count[base]; e++)
        if (hidden_edges[base][e].previous == prev)
          bwd_edges[base][bc++] = hidden_edges[base][e];
    }
    bwd_offset[base][HIDDEN_STATES] = bc;

    int tc = 0;
    for (int tr = 0; tr < ntrans; tr++) {
      tr_offset[base][tr] = tc;
      for (int e = 0; e < hidden_edge_count[base]; e++)
        if (hidden_edges[base][e].transition == tr)
          tr_edges[base][tc++] = hidden_edges[base][e];
    }
    tr_offset[base][ntrans] = tc;
  }
#ifdef _OPENACC
  #pragma acc enter data copyin(hidden_offsets, hidden_counts, trans_a, trans_b, ntrans, \
                                edge_pair, edge_rev, edge_kind, pair_to_edge, \
                                hidden_state, hidden_context, hidden_edges, hidden_edge_count, \
                                initial_states, full_init_states, term_states, \
                                fwd_edges, fwd_offset, bwd_edges, bwd_offset, tr_edges, tr_offset)
#endif
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

#ifndef _OPENACC
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
#endif

static double log_add(double left, double right) {
  if (left == -INFINITY) return right;
  if (right == -INFINITY) return left;
  double maximum = fmax(left, right);
  return maximum + log1p(exp(-fabs(left - right)));
}

#ifdef _OPENACC
#pragma acc routine seq
#endif
static int site_scores_raw(
  const double *w, const char *seq, const float *cwt_power, int genome_length,
  int genome_position, const Site *site, int ntrans_val,
  const int *epair, const int *erev, const int *ekind, const int *tb,
  double *emissions, double *psi
) {
  float deriv_pos[WAVE_COUNT];
  for (int scale = 0; scale < WAVE_COUNT; scale++)
    deriv_pos[scale] = site->derivative[scale] / (float)INT16_MAX;

  for (int state = 0; state < STATES; state++) {
    int reverse = state >= 7;
    double score = w[OFF_ENUC + state * KNUC + site->nuc] +
                   w[OFF_E4 + state * K4 + site->k4_emit[reverse]];
    const double *wp = w + OFF_EP + state * WAVE_COUNT;
    const double *wd = w + OFF_ED + state * WAVE_COUNT;
    for (int scale = 0; scale < WAVE_COUNT; scale++) {
      score += wp[scale] * site->cwt_power[scale];
      score += wd[scale] * (reverse ? -deriv_pos[scale] : deriv_pos[scale]);
    }
    if (!isfinite(score)) return -1;
    emissions[state] = score;
  }

  if (!psi) return 0;

  double b_splice[2][4], b_cwt[2][4];
  for (int rev = 0; rev < 2; rev++) {
    float power[FLANKS][WAVE_COUNT];
    double derivative[FLANKS][WAVE_COUNT];
    for (int offset = 0; offset < FLANKS; offset++)
      boundary_cwt_raw(seq, cwt_power, genome_length, genome_position, rev, offset,
                       power[offset], derivative[offset]);
    const uint8_t *fl = site->flank[rev];
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

  for (int edge = 0; edge < ntrans_val; edge++) {
    int pair = epair[edge], reverse = erev[edge], kind = ekind[edge];
    double score = w[OFF_TNUC + pair * KNUC + site->nuc] +
                   w[OFF_T4 + pair * K4 + site->k4_trans[reverse]];
    if (kind >= 0) score += b_splice[reverse][kind];
    const double *wp = w + OFF_TP + pair * WAVE_COUNT;
    const double *wd = w + OFF_TD + pair * WAVE_COUNT;
    for (int scale = 0; scale < WAVE_COUNT; scale++) {
      score += wp[scale] * site->cwt_power[scale];
      score += wd[scale] * (reverse ? -deriv_pos[scale] : deriv_pos[scale]);
    }
    if (kind >= 0) score += b_cwt[reverse][kind];
    double p = score + emissions[tb[edge]];
    if (!isfinite(p)) return -2;
    psi[edge] = p;
  }
  return 0;
}

#ifndef _OPENACC
static Site site_scores(const double *w, const Contig *c, int position,
                        double *emissions, double *psi) {
  Site site = site_at(c, position);
  int status = site_scores_raw(w, c->seq, c->cwt_power, c->genome_length,
                               c->offset + position, &site, ntrans,
                               edge_pair, edge_rev, edge_kind, trans_b,
                               emissions, psi);
  if (status == -1) anno_fail("Non-finite CRF emission at position %d", position + 1);
  if (status == -2) anno_fail("Non-finite CRF transition at position %d", position + 1);
  return site;
}
#endif

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

#ifndef _OPENACC
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
static double contig_gradient(const double *w, const Contig *c, Work *k, double *gradient) {
  double mass = gold_grammar_mass(c, 0, c->length);
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
#endif

typedef struct {
  Contig *contig;
  int contig_index, start, length;
  double mass;
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
      Contig instance = *contig;
      instance.offset = start;
      instance.length = end - start;
      double mass = gold_grammar_mass(&instance, 0, instance.length);
      (*windows)[window_count++] = (CrfWindow){contig, i, start, end - start, mass};
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

#ifdef _OPENACC
typedef struct {
  size_t total_seq_len, valid_wcount, skipped, positions;
  size_t max_batch_pos, max_batch_win, batch_count;
  size_t *batch_win_offset, *batch_win_count, *batch_pos_count;
  size_t *contig_seq_offset;
  char *flat_seq;
  float *flat_cwt;
  uint8_t *flat_label;
  int *win_buf_start, *win_seq_offset, *win_genome_len, *win_contig_start, *win_len;
  double *win_mass, *win_ll;
  uint32_t *pos_win;
  Site *sites;
  double *exp_psi, *alpha, *beta_buf, *inv_scale, *max_psi_arr, *gold_arr;
  double *g_priv, *g_gpu, *w_dev;
  int *err_flag;
} CrfGpu;

static void crf_gpu_prepare_batch(CrfGpu *gpu, const CrfWindow *windows, size_t batch_idx) {
  size_t w0 = gpu->batch_win_offset[batch_idx];
  size_t wc = gpu->batch_win_count[batch_idx];
  int cursor = 0;
  size_t valid_seen = 0, b_slot = 0;
  for (size_t i = 0; b_slot < wc; i++) {
    if (windows[i].mass == -INFINITY) continue;
    if (valid_seen >= w0) {
      int ci = windows[i].contig_index;
      int len = windows[i].length;
      gpu->win_buf_start[b_slot] = cursor;
      gpu->win_seq_offset[b_slot] = (int)gpu->contig_seq_offset[ci];
      gpu->win_genome_len[b_slot] = windows[i].contig->genome_length;
      gpu->win_contig_start[b_slot] = windows[i].start;
      gpu->win_len[b_slot] = len;
      gpu->win_mass[b_slot] = windows[i].mass;
      for (int p = 0; p < len; p++) gpu->pos_win[cursor + p] = (uint32_t)b_slot;
      cursor += len;
      b_slot++;
    }
    valid_seen++;
  }
  size_t b_pos = (size_t)cursor;
  int *win_buf_start = gpu->win_buf_start;
  int *win_seq_offset = gpu->win_seq_offset;
  int *win_genome_len = gpu->win_genome_len;
  int *win_contig_start = gpu->win_contig_start;
  int *win_len = gpu->win_len;
  double *win_mass = gpu->win_mass;
  uint32_t *pos_win = gpu->pos_win;
  const char *flat_seq = gpu->flat_seq;
  const float *flat_cwt = gpu->flat_cwt;
  Site *sites = gpu->sites;
  size_t total_seq = gpu->total_seq_len;
  size_t total_cwt = total_seq * WAVE_COUNT;
  (void)win_len, (void)win_mass;

  #pragma acc update device(win_buf_start[0:wc], win_seq_offset[0:wc], win_genome_len[0:wc], \
                            win_contig_start[0:wc], win_len[0:wc], win_mass[0:wc], pos_win[0:b_pos])
  #pragma acc parallel loop gang vector present(flat_seq[0:total_seq], flat_cwt[0:total_cwt], \
                                                win_buf_start[0:wc], win_seq_offset[0:wc], \
                                                win_genome_len[0:wc], win_contig_start[0:wc], \
                                                pos_win[0:b_pos], sites[0:b_pos])
  for (size_t idx = 0; idx < b_pos; idx++) {
    uint32_t s = pos_win[idx];
    int p = (int)idx - win_buf_start[s];
    int seq_off = win_seq_offset[s];
    int g_len = win_genome_len[s];
    int g_pos = win_contig_start[s] + p;
    sites[idx] = site_at_raw(flat_seq + seq_off, flat_cwt + (size_t)seq_off * WAVE_COUNT,
                             g_len, g_pos);
  }
}

static void crf_gpu_init(CrfGpu *gpu, const Contig *contigs, int count,
                         const CrfWindow *windows, size_t window_count) {
  memset(gpu, 0, sizeof(*gpu));
  gpu->contig_seq_offset = anno_alloc((size_t)count, sizeof(*gpu->contig_seq_offset));
  size_t total_seq = 0;
  for (int c = 0; c < count; c++) {
    gpu->contig_seq_offset[c] = total_seq;
    total_seq += (size_t)contigs[c].length;
  }
  gpu->total_seq_len = total_seq;
  gpu->flat_seq = anno_alloc(total_seq, 1);
  gpu->flat_cwt = anno_alloc(total_seq * WAVE_COUNT, sizeof(float));
  gpu->flat_label = anno_alloc(total_seq, 1);
  for (int c = 0; c < count; c++) {
    size_t off = gpu->contig_seq_offset[c], len = (size_t)contigs[c].length;
    memcpy(gpu->flat_seq + off, contigs[c].seq, len);
    memcpy(gpu->flat_cwt + off * WAVE_COUNT, contigs[c].cwt_power, len * WAVE_COUNT * sizeof(float));
    memcpy(gpu->flat_label + off, contigs[c].label, len);
  }

  size_t valid_wcount = 0, skipped = 0, positions = 0;
  for (size_t i = 0; i < window_count; i++) {
    if (windows[i].mass == -INFINITY) skipped++;
    else {
      valid_wcount++;
      positions += (size_t)windows[i].length;
    }
  }
  if (!positions) anno_fail("No valid CRF training windows (%zu windows skipped)", skipped);
  gpu->valid_wcount = valid_wcount;
  gpu->skipped = skipped;
  gpu->positions = positions;

  gpu->batch_win_offset = anno_alloc(valid_wcount, sizeof(size_t));
  gpu->batch_win_count = anno_alloc(valid_wcount, sizeof(size_t));
  gpu->batch_pos_count = anno_alloc(valid_wcount, sizeof(size_t));

  size_t batch_count = 0, cur_w0 = 0, cur_wc = 0, cur_pos = 0;
  size_t max_batch_pos = 0, max_batch_win = 0;
  for (size_t i = 0; i < window_count; i++) {
    if (windows[i].mass == -INFINITY) continue;
    size_t wlen = (size_t)windows[i].length;
    if (cur_wc > 0 && cur_pos + wlen > CRF_GPU_BATCH_POS) {
      gpu->batch_win_offset[batch_count] = cur_w0;
      gpu->batch_win_count[batch_count] = cur_wc;
      gpu->batch_pos_count[batch_count] = cur_pos;
      if (cur_pos > max_batch_pos) max_batch_pos = cur_pos;
      if (cur_wc > max_batch_win) max_batch_win = cur_wc;
      batch_count++;
      cur_w0 += cur_wc;
      cur_wc = 0;
      cur_pos = 0;
    }
    cur_wc++;
    cur_pos += wlen;
  }
  if (cur_wc > 0) {
    gpu->batch_win_offset[batch_count] = cur_w0;
    gpu->batch_win_count[batch_count] = cur_wc;
    gpu->batch_pos_count[batch_count] = cur_pos;
    if (cur_pos > max_batch_pos) max_batch_pos = cur_pos;
    if (cur_wc > max_batch_win) max_batch_win = cur_wc;
    batch_count++;
  }
  gpu->batch_count = batch_count;
  gpu->max_batch_pos = max_batch_pos;
  gpu->max_batch_win = max_batch_win;

  gpu->win_buf_start = anno_alloc(max_batch_win, sizeof(int));
  gpu->win_seq_offset = anno_alloc(max_batch_win, sizeof(int));
  gpu->win_genome_len = anno_alloc(max_batch_win, sizeof(int));
  gpu->win_contig_start = anno_alloc(max_batch_win, sizeof(int));
  gpu->win_len = anno_alloc(max_batch_win, sizeof(int));
  gpu->win_mass = anno_alloc(max_batch_win, sizeof(double));
  gpu->win_ll = anno_alloc(max_batch_win, sizeof(double));
  gpu->pos_win = anno_alloc(max_batch_pos, sizeof(uint32_t));
  gpu->sites = anno_alloc(max_batch_pos, sizeof(Site));
  gpu->exp_psi = anno_alloc(max_batch_pos, (size_t)ntrans * sizeof(double));
  gpu->alpha = anno_alloc(max_batch_pos, HIDDEN_STATES * sizeof(double));
  gpu->beta_buf = anno_alloc(max_batch_win * 2, HIDDEN_STATES * sizeof(double));
  gpu->inv_scale = anno_alloc(max_batch_pos, sizeof(double));
  gpu->max_psi_arr = anno_alloc(max_batch_pos, sizeof(double));
  gpu->gold_arr = anno_alloc(max_batch_pos, sizeof(double));
  gpu->g_priv = anno_alloc((size_t)WEIGHTS * CRF_GPU_THREADS, sizeof(double));
  gpu->g_gpu = anno_alloc(WEIGHTS, sizeof(double));
  gpu->w_dev = anno_alloc(WEIGHTS, sizeof(double));
  gpu->err_flag = anno_alloc(1, sizeof(int));

  char *flat_seq = gpu->flat_seq;
  float *flat_cwt = gpu->flat_cwt;
  uint8_t *flat_label = gpu->flat_label;
  int *win_buf_start = gpu->win_buf_start;
  int *win_seq_offset = gpu->win_seq_offset;
  int *win_genome_len = gpu->win_genome_len;
  int *win_contig_start = gpu->win_contig_start;
  int *win_len = gpu->win_len;
  double *win_mass = gpu->win_mass;
  double *win_ll = gpu->win_ll;
  uint32_t *pos_win = gpu->pos_win;
  Site *sites = gpu->sites;
  double *exp_psi = gpu->exp_psi;
  double *alpha = gpu->alpha;
  double *beta_buf = gpu->beta_buf;
  double *inv_scale = gpu->inv_scale;
  double *max_psi_arr = gpu->max_psi_arr;
  double *gold_arr = gpu->gold_arr;
  double *g_priv = gpu->g_priv;
  double *g_gpu = gpu->g_gpu;
  double *w_dev = gpu->w_dev;
  int *err_flag = gpu->err_flag;
  size_t total_cwt = total_seq * WAVE_COUNT;
  size_t ep_len = max_batch_pos * (size_t)ntrans;
  size_t al_len = max_batch_pos * HIDDEN_STATES;
  size_t bb_len = max_batch_win * 2 * HIDDEN_STATES;
  size_t gp_len = (size_t)WEIGHTS * CRF_GPU_THREADS;

  #pragma acc enter data copyin(flat_seq[0:total_seq], flat_cwt[0:total_cwt], flat_label[0:total_seq]) \
                         create(win_buf_start[0:max_batch_win], win_seq_offset[0:max_batch_win], \
                                win_genome_len[0:max_batch_win], win_contig_start[0:max_batch_win], \
                                win_len[0:max_batch_win], win_mass[0:max_batch_win], win_ll[0:max_batch_win], \
                                pos_win[0:max_batch_pos], sites[0:max_batch_pos], exp_psi[0:ep_len], \
                                alpha[0:al_len], beta_buf[0:bb_len], inv_scale[0:max_batch_pos], \
                                max_psi_arr[0:max_batch_pos], gold_arr[0:max_batch_pos], \
                                g_priv[0:gp_len], g_gpu[0:WEIGHTS], w_dev[0:WEIGHTS], err_flag[0:1])

  (void)flat_seq, (void)flat_cwt, (void)flat_label, (void)win_buf_start, (void)win_seq_offset;
  (void)win_genome_len, (void)win_contig_start, (void)win_len, (void)win_mass, (void)win_ll;
  (void)pos_win, (void)sites, (void)exp_psi, (void)alpha, (void)beta_buf, (void)inv_scale;
  (void)max_psi_arr, (void)gold_arr, (void)g_priv, (void)g_gpu, (void)w_dev, (void)err_flag;
  (void)total_cwt, (void)ep_len, (void)al_len, (void)bb_len, (void)gp_len;

  if (batch_count == 1) crf_gpu_prepare_batch(gpu, windows, 0);
}

static void crf_gpu_free(CrfGpu *gpu) {
  char *flat_seq = gpu->flat_seq;
  float *flat_cwt = gpu->flat_cwt;
  uint8_t *flat_label = gpu->flat_label;
  int *win_buf_start = gpu->win_buf_start;
  int *win_seq_offset = gpu->win_seq_offset;
  int *win_genome_len = gpu->win_genome_len;
  int *win_contig_start = gpu->win_contig_start;
  int *win_len = gpu->win_len;
  double *win_mass = gpu->win_mass;
  double *win_ll = gpu->win_ll;
  uint32_t *pos_win = gpu->pos_win;
  Site *sites = gpu->sites;
  double *exp_psi = gpu->exp_psi;
  double *alpha = gpu->alpha;
  double *beta_buf = gpu->beta_buf;
  double *inv_scale = gpu->inv_scale;
  double *max_psi_arr = gpu->max_psi_arr;
  double *gold_arr = gpu->gold_arr;
  double *g_priv = gpu->g_priv;
  double *g_gpu = gpu->g_gpu;
  double *w_dev = gpu->w_dev;
  int *err_flag = gpu->err_flag;
  size_t total_seq = gpu->total_seq_len;
  size_t total_cwt = total_seq * WAVE_COUNT;
  size_t max_batch_win = gpu->max_batch_win;
  size_t max_batch_pos = gpu->max_batch_pos;
  size_t ep_len = max_batch_pos * (size_t)ntrans;
  size_t al_len = max_batch_pos * HIDDEN_STATES;
  size_t bb_len = max_batch_win * 2 * HIDDEN_STATES;
  size_t gp_len = (size_t)WEIGHTS * CRF_GPU_THREADS;

  #pragma acc exit data delete(flat_seq[0:total_seq], flat_cwt[0:total_cwt], flat_label[0:total_seq], \
                               win_buf_start[0:max_batch_win], win_seq_offset[0:max_batch_win], \
                               win_genome_len[0:max_batch_win], win_contig_start[0:max_batch_win], \
                               win_len[0:max_batch_win], win_mass[0:max_batch_win], win_ll[0:max_batch_win], \
                               pos_win[0:max_batch_pos], sites[0:max_batch_pos], exp_psi[0:ep_len], \
                               alpha[0:al_len], beta_buf[0:bb_len], inv_scale[0:max_batch_pos], \
                               max_psi_arr[0:max_batch_pos], gold_arr[0:max_batch_pos], \
                               g_priv[0:gp_len], g_gpu[0:WEIGHTS], w_dev[0:WEIGHTS], err_flag[0:1])
  (void)flat_seq, (void)flat_cwt, (void)flat_label, (void)win_buf_start, (void)win_seq_offset;
  (void)win_genome_len, (void)win_contig_start, (void)win_len, (void)win_mass, (void)win_ll;
  (void)pos_win, (void)sites, (void)exp_psi, (void)alpha, (void)beta_buf, (void)inv_scale;
  (void)max_psi_arr, (void)gold_arr, (void)g_priv, (void)g_gpu, (void)w_dev, (void)err_flag;
  (void)total_cwt, (void)ep_len, (void)al_len, (void)bb_len, (void)gp_len;
  free(gpu->contig_seq_offset);
  free(gpu->batch_win_offset), free(gpu->batch_win_count), free(gpu->batch_pos_count);
  free(gpu->flat_seq), free(gpu->flat_cwt), free(gpu->flat_label);
  free(gpu->win_buf_start), free(gpu->win_seq_offset), free(gpu->win_genome_len);
  free(gpu->win_contig_start), free(gpu->win_len), free(gpu->win_mass), free(gpu->win_ll);
  free(gpu->pos_win), free(gpu->sites), free(gpu->exp_psi), free(gpu->alpha);
  free(gpu->beta_buf), free(gpu->inv_scale), free(gpu->max_psi_arr), free(gpu->gold_arr);
  free(gpu->g_priv), free(gpu->g_gpu), free(gpu->w_dev), free(gpu->err_flag);
}
#else
typedef struct {
  CrfTrain *train;
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
  size_t slot = atomic_fetch_add_explicit(&train->next_instance, 1, memory_order_relaxed);
  while (slot < train->count) {
    const CrfWindow *window = train->windows + slot;
    Contig instance = *window->contig;
    instance.offset = window->start;
    instance.length = window->length;
    double score = contig_gradient(worker->w, &instance, &worker->work, worker->gradient);
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
#endif

typedef struct {
  CrfTrain train;
#ifdef _OPENACC
  CrfGpu gpu;
#else
  CrfWorker *workers;
  pthread_t *threads;
#endif
  double likelihood;
  size_t positions, skipped;
  double *scale, *weights;
} CrfObjective;

static double crf_objective(const double *w, double *gradient, void *context) {
  CrfObjective *objective = context;
  CrfTrain *train = &objective->train;
  double likelihood = 0;
  size_t positions = 0, skipped = 0;

#ifdef _OPENACC
  CrfGpu *gpu = &objective->gpu;
  memcpy(gpu->w_dev, w, WEIGHTS * sizeof(double));
  gpu->err_flag[0] = 0;

  const char *flat_seq = gpu->flat_seq;
  const float *flat_cwt = gpu->flat_cwt;
  const uint8_t *flat_label = gpu->flat_label;
  int *win_buf_start = gpu->win_buf_start;
  int *win_seq_offset = gpu->win_seq_offset;
  int *win_genome_len = gpu->win_genome_len;
  int *win_contig_start = gpu->win_contig_start;
  int *win_len = gpu->win_len;
  double *win_mass = gpu->win_mass;
  double *win_ll = gpu->win_ll;
  uint32_t *pos_win = gpu->pos_win;
  Site *sites = gpu->sites;
  double *exp_psi = gpu->exp_psi;
  double *alpha = gpu->alpha;
  double *beta_buf = gpu->beta_buf;
  double *inv_scale = gpu->inv_scale;
  double *max_psi_arr = gpu->max_psi_arr;
  double *gold_arr = gpu->gold_arr;
  double *g_priv = gpu->g_priv;
  double *g_gpu = gpu->g_gpu;
  double *w_dev = gpu->w_dev;
  int *err_flag = gpu->err_flag;
  size_t total_seq = gpu->total_seq_len;
  size_t total_cwt = total_seq * WAVE_COUNT;
  size_t gp_len = (size_t)WEIGHTS * CRF_GPU_THREADS;

  #pragma acc update device(w_dev[0:WEIGHTS], err_flag[0:1])
  #pragma acc parallel loop gang vector present(g_priv[0:gp_len], g_gpu[0:WEIGHTS])
  for (size_t i = 0; i < gp_len; i++) {
    g_priv[i] = 0.0;
    if (i < WEIGHTS) g_gpu[i] = 0.0;
  }

  for (size_t b = 0; b < gpu->batch_count; b++) {
    if (gpu->batch_count > 1) crf_gpu_prepare_batch(gpu, train->windows, b);
    size_t wc = gpu->batch_win_count[b];
    size_t b_pos = gpu->batch_pos_count[b];
    size_t ep_len = b_pos * (size_t)ntrans;
    size_t al_len = b_pos * HIDDEN_STATES;
    size_t bb_len = wc * 2 * HIDDEN_STATES;

    /* Kernel 1: Site scores across all positions in batch */
    #pragma acc parallel loop gang vector present(trans_b, ntrans, edge_pair, edge_rev, edge_kind, pair_to_edge, \
                                                  w_dev[0:WEIGHTS], flat_seq[0:total_seq], flat_cwt[0:total_cwt], \
                                                  flat_label[0:total_seq], win_buf_start[0:wc], win_seq_offset[0:wc], \
                                                  win_genome_len[0:wc], win_contig_start[0:wc], pos_win[0:b_pos], \
                                                  sites[0:b_pos], exp_psi[0:ep_len], alpha[0:al_len], \
                                                  max_psi_arr[0:b_pos], gold_arr[0:b_pos], err_flag[0:1])
    for (size_t idx = 0; idx < b_pos; idx++) {
      uint32_t s = pos_win[idx];
      int p = (int)idx - win_buf_start[s];
      int seq_off = win_seq_offset[s];
      int g_len = win_genome_len[s];
      int g_pos = win_contig_start[s] + p;
      Site site = sites[idx];
      double emissions[STATES], psi[STATES * STATES];
      int status = site_scores_raw(w_dev, flat_seq + seq_off, flat_cwt + (size_t)seq_off * WAVE_COUNT,
                                   g_len, g_pos, &site, ntrans,
                                   edge_pair, edge_rev, edge_kind, trans_b,
                                   emissions, p ? psi : NULL);
      if (status < 0) {
        err_flag[0] = status;
        continue;
      }
      int gold_state = flat_label[seq_off + g_pos];
      double gs = emissions[gold_state];
      if (!p) {
        double *A = alpha + idx * HIDDEN_STATES;
        for (int st = 0; st < STATES; st++) A[st] = emissions[st];
        max_psi_arr[idx] = 0.0;
        gold_arr[idx] = gs;
      } else {
        double max_psi = psi[0];
        for (int edge = 1; edge < ntrans; edge++)
          if (psi[edge] > max_psi) max_psi = psi[edge];
        max_psi_arr[idx] = max_psi;
        double *ep = exp_psi + idx * (size_t)ntrans;
        for (int edge = 0; edge < ntrans; edge++) ep[edge] = exp(psi[edge] - max_psi);
        int prev_gold = flat_label[seq_off + g_pos - 1];
        int gold_edge = pair_to_edge[prev_gold][gold_state];
        if (gold_edge >= 0) gs += psi[gold_edge] - emissions[gold_state];
        gold_arr[idx] = gs;
      }
    }

    /* Kernel 2: Forward-Backward and marginals per window */
    #pragma acc parallel loop gang present(hidden_offsets, hidden_counts, fwd_edges, fwd_offset, \
                                           bwd_edges, bwd_offset, tr_edges, tr_offset, hidden_state, \
                                           initial_states, full_init_states, term_states, ntrans, pair_to_edge, \
                                           flat_label[0:total_seq], win_buf_start[0:wc], win_seq_offset[0:wc], \
                                           win_genome_len[0:wc], win_contig_start[0:wc], win_len[0:wc], \
                                           win_mass[0:wc], win_ll[0:wc], sites[0:b_pos], exp_psi[0:ep_len], \
                                           alpha[0:al_len], beta_buf[0:bb_len], inv_scale[0:b_pos], \
                                           max_psi_arr[0:b_pos], gold_arr[0:b_pos])
    for (size_t s = 0; s < wc; s++) {
      int buf_start = win_buf_start[s], len = win_len[s];
      int seq_off = win_seq_offset[s], c_start = win_contig_start[s], g_len = win_genome_len[s];
      int reaches_end = (c_start + len == g_len);
      double logz = 0, gold_score = gold_arr[buf_start];

      double *A0 = alpha + (size_t)buf_start * HIDDEN_STATES;
      double *B = beta_buf + s * 2 * HIDDEN_STATES;
      double *B_prev = B + HIDDEN_STATES;
      #pragma acc loop vector
      for (int st = 0; st < STATES; st++) B[st] = A0[st];

      int b0 = sites[buf_start].nuc;
      double max_e = -INFINITY;
      #pragma acc loop vector reduction(max:max_e)
      for (int h = 0; h < HIDDEN_STATES; h++) {
        int ok = (c_start == 0) ? full_init_states[b0][h] : initial_states[b0][h];
        if (ok && B[hidden_state[h]] > max_e) max_e = B[hidden_state[h]];
      }
      double sum0 = 0;
      #pragma acc loop vector reduction(+:sum0)
      for (int h = 0; h < HIDDEN_STATES; h++) {
        int ok = (c_start == 0) ? full_init_states[b0][h] : initial_states[b0][h];
        double v = ok ? exp(B[hidden_state[h]] - max_e) : 0.0;
        if (len == 1 && reaches_end && !term_states[h]) v = 0.0;
        A0[h] = v;
        sum0 += v;
      }
      double inv0 = 1.0 / sum0;
      #pragma acc loop vector
      for (int h = 0; h < HIDDEN_STATES; h++) A0[h] *= inv0;
      inv_scale[buf_start] = inv0;
      logz += max_e + log(sum0);

      #pragma acc loop seq
      for (int p = 1; p < len; p++) {
        size_t idx = (size_t)(buf_start + p);
        int nb = sites[idx].nuc;
        double *A = alpha + idx * HIDDEN_STATES;
        const double *prev_A = A - HIDDEN_STATES;
        const double *ep = exp_psi + idx * (size_t)ntrans;
        double sum = 0;
        #pragma acc loop vector reduction(+:sum)
        for (int h = 0; h < HIDDEN_STATES; h++) {
          double v = 0;
          int e0 = fwd_offset[nb][h], e1 = fwd_offset[nb][h + 1];
          #pragma acc loop seq
          for (int e = e0; e < e1; e++)
            v += prev_A[fwd_edges[nb][e].previous] * ep[fwd_edges[nb][e].transition];
          if (p == len - 1 && reaches_end && !term_states[h]) v = 0.0;
          A[h] = v;
          sum += v;
        }
        double inv = 1.0 / sum;
        #pragma acc loop vector
        for (int h = 0; h < HIDDEN_STATES; h++) A[h] *= inv;
        inv_scale[idx] = inv;
        logz += max_psi_arr[idx] + log(sum);
        gold_score += gold_arr[idx];
      }

      #pragma acc loop vector
      for (int h = 0; h < HIDDEN_STATES; h++)
        B[h] = (!reaches_end || term_states[h]) ? 1.0 : 0.0;

      #pragma acc loop seq
      for (int p = len - 1; p >= 0; p--) {
        size_t idx = (size_t)(buf_start + p);
        int nb = sites[idx].nuc;
        double *A = alpha + idx * HIDDEN_STATES;
        int state = flat_label[seq_off + c_start + p];

        if (p > 0) {
          const double *prev_A = A - HIDDEN_STATES;
          const double *ep = exp_psi + idx * (size_t)ntrans;
          double inv = inv_scale[idx];
          int prev_state = flat_label[seq_off + c_start + p - 1];
          int gold_edge = pair_to_edge[prev_state][state];

          #pragma acc loop vector
          for (int h = 0; h < HIDDEN_STATES; h++) {
            double v = 0;
            int e0 = bwd_offset[nb][h], e1 = bwd_offset[nb][h + 1];
            #pragma acc loop seq
            for (int e = e0; e < e1; e++)
              v += ep[bwd_edges[nb][e].transition] * B[bwd_edges[nb][e].current];
            B_prev[h] = v * inv;
          }

          double *em_out = exp_psi + idx * (size_t)ntrans;
          #pragma acc loop vector
          for (int tr = 0; tr < ntrans; tr++) {
            double em = 0;
            int e0 = tr_offset[nb][tr], e1 = tr_offset[nb][tr + 1];
            #pragma acc loop seq
            for (int e = e0; e < e1; e++)
              em += prev_A[tr_edges[nb][e].previous] * B[tr_edges[nb][e].current];
            double val = em * ep[tr] * inv;
            if (tr == gold_edge) val -= 1.0;
            em_out[tr] = val;
          }
        }

        #pragma acc loop vector
        for (int st = 0; st < STATES; st++) {
          double sm = 0;
          int h0 = hidden_offsets[st], h1 = h0 + hidden_counts[st];
          #pragma acc loop seq
          for (int h = h0; h < h1; h++) sm += A[h] * B[h];
          if (st == state) sm -= 1.0;
          A[st] = sm;
        }

        if (p > 0) {
          #pragma acc loop vector
          for (int h = 0; h < HIDDEN_STATES; h++) B[h] = B_prev[h];
        }
      }
      win_ll[s] = win_mass[s] - logz + gold_score;
    }

    /* Kernel 3a: Discrete feature gradients (zero-atomic strided accumulation) */
    #pragma acc parallel loop gang num_gangs(CRF_GPU_GANGS) vector_length(CRF_GPU_VEC) \
                              present(ntrans, edge_pair, edge_rev, edge_kind, win_buf_start[0:wc], \
                                      pos_win[0:b_pos], sites[0:b_pos], alpha[0:al_len], \
                                      exp_psi[0:ep_len], g_priv[0:gp_len])
    for (int g = 0; g < CRF_GPU_GANGS; g++) {
      #pragma acc loop vector
      for (int v = 0; v < CRF_GPU_VEC; v++) {
        int tid = g * CRF_GPU_VEC + v;
        double *gp = g_priv + tid;
        #pragma acc loop seq
        for (size_t idx = (size_t)tid; idx < b_pos; idx += CRF_GPU_THREADS) {
          uint32_t s = pos_win[idx];
          int is_start = ((int)idx == win_buf_start[s]);
          Site site = sites[idx];
          double *A = alpha + idx * HIDDEN_STATES;
          #pragma acc loop seq
          for (int st = 0; st < STATES; st++) {
            double val = -A[st];
            gp[(OFF_ENUC + st * KNUC + site.nuc) * CRF_GPU_THREADS] += val;
            gp[(OFF_E4 + st * K4 + site.k4_emit[st >= 7]) * CRF_GPU_THREADS] += val;
          }
          if (!is_start) {
            const double *em = exp_psi + idx * (size_t)ntrans;
            double b_coeff[4][2] = {{0}};
            #pragma acc loop seq
            for (int edge = 0; edge < ntrans; edge++) {
              double val = -em[edge];
              int pair = edge_pair[edge], reverse = edge_rev[edge], kind = edge_kind[edge];
              gp[(OFF_TNUC + pair * KNUC + site.nuc) * CRF_GPU_THREADS] += val;
              gp[(OFF_T4 + pair * K4 + site.k4_trans[reverse]) * CRF_GPU_THREADS] += val;
              if (kind >= 0) b_coeff[kind][reverse] += val;
            }
            #pragma acc loop seq
            for (int kind = 0; kind < 4; kind++) {
              #pragma acc loop seq
              for (int rev = 0; rev < 2; rev++) {
                double val = b_coeff[kind][rev];
                A[STATES + kind * 2 + rev] = val;
                if (val != 0) {
                  #pragma acc loop seq
                  for (int offset = 0; offset < FLANKS; offset++)
                    gp[(OFF_SPLICE + (kind * FLANKS + offset) * K2 + site.flank[rev][offset]) * CRF_GPU_THREADS] += val;
                }
              }
            }
          } else {
            #pragma acc loop seq
            for (int i = 0; i < 8; i++) A[STATES + i] = 0.0;
          }
        }
      }
    }

    /* Kernel 3b: Continuous CWT feature gradients via register reductions */
    #pragma acc parallel loop gang collapse(2) present(sites[0:b_pos], alpha[0:al_len], g_gpu[0:WEIGHTS])
    for (int st = 0; st < STATES; st++) {
      for (int scale = 0; scale < WAVE_COUNT; scale++) {
        int reverse = st >= 7;
        double sum_p = 0, sum_d = 0;
        #pragma acc loop vector reduction(+:sum_p, sum_d)
        for (size_t idx = 0; idx < b_pos; idx++) {
          double val = -alpha[idx * HIDDEN_STATES + st];
          float d = sites[idx].derivative[scale] / (float)INT16_MAX;
          sum_p += val * sites[idx].cwt_power[scale];
          sum_d += val * (reverse ? -d : d);
        }
        g_gpu[OFF_EP + st * WAVE_COUNT + scale] += sum_p;
        g_gpu[OFF_ED + st * WAVE_COUNT + scale] += sum_d;
      }
    }

    #pragma acc parallel loop gang collapse(2) present(ntrans, edge_pair, edge_rev, win_buf_start[0:wc], \
                                                       pos_win[0:b_pos], sites[0:b_pos], exp_psi[0:ep_len], g_gpu[0:WEIGHTS])
    for (int edge = 0; edge < ntrans; edge++) {
      for (int scale = 0; scale < WAVE_COUNT; scale++) {
        int pair = edge_pair[edge], reverse = edge_rev[edge];
        double sum_p = 0, sum_d = 0;
        #pragma acc loop vector reduction(+:sum_p, sum_d)
        for (size_t idx = 0; idx < b_pos; idx++) {
          uint32_t s = pos_win[idx];
          if ((int)idx != win_buf_start[s]) {
            double val = -exp_psi[idx * (size_t)ntrans + edge];
            float d = sites[idx].derivative[scale] / (float)INT16_MAX;
            sum_p += val * sites[idx].cwt_power[scale];
            sum_d += val * (reverse ? -d : d);
          }
        }
        g_gpu[OFF_TP + pair * WAVE_COUNT + scale] += sum_p;
        g_gpu[OFF_TD + pair * WAVE_COUNT + scale] += sum_d;
      }
    }

    #pragma acc parallel loop gang collapse(3) present(flat_seq[0:total_seq], flat_cwt[0:total_cwt], \
                                                       win_buf_start[0:wc], win_seq_offset[0:wc], \
                                                       win_genome_len[0:wc], win_contig_start[0:wc], \
                                                       pos_win[0:b_pos], alpha[0:al_len], g_gpu[0:WEIGHTS])
    for (int kind = 0; kind < 4; kind++) {
      for (int offset = 0; offset < FLANKS; offset++) {
        for (int scale = 0; scale < WAVE_COUNT; scale++) {
          int displacement = -RADIUS + 2 * offset;
          int base = (kind * FLANKS + offset) * WAVE_COUNT;
          double sum_p = 0, sum_d = 0;
          #pragma acc loop vector reduction(+:sum_p, sum_d)
          for (size_t idx = 0; idx < b_pos; idx++) {
            const double *A = alpha + idx * HIDDEN_STATES;
            uint32_t s = pos_win[idx];
            int p = (int)idx - win_buf_start[s];
            int seq_off = win_seq_offset[s];
            int g_len = win_genome_len[s];
            int g_pos = win_contig_start[s] + p;
            const char *c_seq = flat_seq + seq_off;
            const float *c_cwt = flat_cwt + (size_t)seq_off * WAVE_COUNT;
            #pragma acc loop seq
            for (int rev = 0; rev < 2; rev++) {
              double val = A[STATES + kind * 2 + rev];
              if (val != 0) {
                int g = rev ? g_pos - displacement - 2 : g_pos + displacement;
                if (g >= 0 && g < g_len && base_index(c_seq[g]) >= 0) {
                  sum_p += val * c_cwt[(size_t)g * WAVE_COUNT + scale];
                  if (g > 0 && g + 1 < g_len &&
                      base_index(c_seq[g - 1]) >= 0 && base_index(c_seq[g + 1]) >= 0) {
                    double left = c_cwt[(size_t)(g - 1) * WAVE_COUNT + scale];
                    double right = c_cwt[(size_t)(g + 1) * WAVE_COUNT + scale];
                    double slope = (right - left) / (right + left + 1e-12);
                    sum_d += val * (rev ? -slope : slope);
                  }
                }
              }
            }
          }
          g_gpu[OFF_BP + base + scale] += sum_p;
          g_gpu[OFF_BD + base + scale] += sum_d;
        }
      }
    }

    #pragma acc update self(win_ll[0:wc], err_flag[0:1])
    if (err_flag[0] == -1) anno_fail("Non-finite CRF emission");
    if (err_flag[0] == -2) anno_fail("Non-finite CRF transition");
    for (size_t s = 0; s < wc; s++) {
      if (!isfinite(win_ll[s])) anno_fail("Non-finite training likelihood");
      likelihood += win_ll[s];
    }
  }

  #pragma acc parallel loop gang present(g_priv[0:gp_len], g_gpu[0:WEIGHTS])
  for (int f = 0; f < WEIGHTS; f++) {
    double s = 0;
    #pragma acc loop vector reduction(+:s)
    for (int t = 0; t < CRF_GPU_THREADS; t++) s += g_priv[(size_t)f * CRF_GPU_THREADS + t];
    g_gpu[f] += s;
  }
  #pragma acc update self(g_gpu[0:WEIGHTS])
  memcpy(gradient, g_gpu, WEIGHTS * sizeof(double));
  positions = gpu->positions;
  skipped = gpu->skipped;
#else
  atomic_store_explicit(&train->next_instance, 0, memory_order_relaxed);
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
#endif

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
  atomic_init(&objective.train.next_instance, 0);
#ifdef _OPENACC
  fprintf(stderr, "Using OpenACC GPU (%d gangs x %d vectors) for %zu training windows\n",
          CRF_GPU_GANGS, CRF_GPU_VEC, objective.train.count);
  crf_gpu_init(&objective.gpu, contigs, count, windows, window_count);
#else
  fprintf(stderr, "Using %d worker threads for %zu training windows\n",
          objective.train.threads, objective.train.count);
  objective.workers = anno_alloc((size_t)objective.train.threads, sizeof(*objective.workers));
  objective.threads = anno_alloc((size_t)objective.train.threads, sizeof(*objective.threads));
  for (int t = 0; t < objective.train.threads; t++)
    objective.workers[t] = (CrfWorker){
      .train = &objective.train,
      .gradient = anno_alloc(WEIGHTS, sizeof(*objective.workers->gradient))
    };
#endif

  int completed = 0;
  objective.scale = anno_alloc(WEIGHTS, sizeof(*objective.scale));
  objective.weights = anno_alloc(WEIGHTS, sizeof(*objective.weights));
  feature_scales(contigs, count, objective.scale);
  for (int f = 0; f < WEIGHTS; f++) w[f] /= objective.scale[f];
  LbfgsStatus status = lbfgs_minimize(WEIGHTS, w, crf_scaled_objective, crf_progress, &objective,
                                      iterations, &completed);
  for (int f = 0; f < WEIGHTS; f++) w[f] *= objective.scale[f];
  free(objective.scale), free(objective.weights);
#ifdef _OPENACC
  crf_gpu_free(&objective.gpu);
#else
  for (int t = 0; t < objective.train.threads; t++) {
    work_free(&objective.workers[t].work);
    free(objective.workers[t].gradient);
  }
  free(objective.workers), free(objective.threads);
#endif
  free(windows);
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

/* Exact whole-contig MAP decoding. */
void crf_decode(const double *w, const Contig *c, uint8_t *path) {
  int length = c->length;
  uint8_t *back = anno_alloc((size_t)length, HIDDEN_STATES * sizeof(*back));
  const char *seq = c->seq;
  const float *cpow = c->cwt_power;
  int chunk_cap = length < CRF_DECODE_CHUNK ? length : CRF_DECODE_CHUNK;
  double *chunk_psi = anno_alloc((size_t)chunk_cap, (size_t)ntrans * sizeof(double));
  uint8_t *chunk_nuc = anno_alloc((size_t)chunk_cap, sizeof(uint8_t));
  int err_flag[1] = {0};

  Site site0 = site_at(c, 0);
  double init_emit[STATES];
  if (site_scores_raw(w, seq, cpow, length, 0, &site0, ntrans,
                      edge_pair, edge_rev, edge_kind, trans_b, init_emit, NULL) < 0)
    anno_fail("Non-finite prediction emission");

  double previous[HIDDEN_STATES] = {0}, current[HIDDEN_STATES];
#ifdef _OPENACC
  size_t total_cwt = (size_t)length * WAVE_COUNT;
  size_t cp_len = (size_t)chunk_cap * (size_t)ntrans;
  #pragma acc data copyin(w[0:WEIGHTS], seq[0:length], cpow[0:total_cwt], err_flag[0:1]) \
                   create(chunk_psi[0:cp_len], chunk_nuc[0:chunk_cap])
#endif
  {
    for (int chunk_start = 0; chunk_start < length; chunk_start += CRF_DECODE_CHUNK) {
      int chunk_len = length - chunk_start > CRF_DECODE_CHUNK ? CRF_DECODE_CHUNK : length - chunk_start;
#ifdef _OPENACC
      #pragma acc parallel loop gang vector present(trans_b, ntrans, edge_pair, edge_rev, edge_kind, \
                                                    w[0:WEIGHTS], seq[0:length], cpow[0:total_cwt], \
                                                    chunk_psi[0:cp_len], chunk_nuc[0:chunk_cap], err_flag[0:1])
#endif
      for (int local = 0; local < chunk_len; local++) {
        int pos = chunk_start + local;
        Site site = site_at_raw(seq, cpow, length, pos);
        chunk_nuc[local] = site.nuc;
        if (!pos) continue;
        double emissions[STATES];
        int status = site_scores_raw(w, seq, cpow, length, pos, &site, ntrans,
                                     edge_pair, edge_rev, edge_kind, trans_b,
                                     emissions, chunk_psi + (size_t)local * (size_t)ntrans);
        if (status < 0) err_flag[0] = status;
      }
#ifdef _OPENACC
      #pragma acc update self(chunk_psi[0:(size_t)chunk_len * (size_t)ntrans], chunk_nuc[0:chunk_len], err_flag[0:1])
#endif
      if (err_flag[0] == -1) anno_fail("Non-finite prediction emission");
      if (err_flag[0] == -2) anno_fail("Non-finite prediction transition");

      for (int local = 0; local < chunk_len; local++) {
        int pos = chunk_start + local;
        int b = chunk_nuc[local];
        uint8_t *bt = back + (size_t)pos * HIDDEN_STATES;
        double maximum = -INFINITY;
        if (!pos) {
          for (int h = 0; h < HIDDEN_STATES; h++) {
            current[h] = full_init_states[b][h] ? init_emit[hidden_state[h]] : -INFINITY;
            bt[h] = 0;
            if (current[h] > maximum) maximum = current[h];
          }
        } else {
          const double *psi = chunk_psi + (size_t)local * (size_t)ntrans;
          for (int h = 0; h < HIDDEN_STATES; h++) {
            double best_s = -INFINITY;
            uint8_t best_p = 0;
            int e0 = fwd_offset[b][h], e1 = fwd_offset[b][h + 1];
            for (int e = e0; e < e1; e++) {
              double sc = previous[fwd_edges[b][e].previous] + psi[fwd_edges[b][e].transition];
              if (sc > best_s) {
                best_s = sc;
                best_p = fwd_edges[b][e].previous;
              }
            }
            current[h] = best_s;
            bt[h] = best_p;
            if (best_s > maximum) maximum = best_s;
          }
        }
        if (!isfinite(maximum)) anno_fail("No finite prediction path");
        for (int h = 0; h < HIDDEN_STATES; h++) current[h] -= maximum;
        memcpy(previous, current, sizeof(previous));
      }
    }
  }

  int best = 0;
  for (int hidden = 1; hidden < HIDDEN_STATES; hidden++)
    if (terminal_allowed(hidden) && previous[hidden] > previous[best]) best = hidden;
  if (!isfinite(previous[best])) anno_fail("No complete prediction path");
  for (int position = length - 1; position >= 0; position--) {
    path[position] = hidden_state[best];
    best = back[(size_t)position * HIDDEN_STATES + best];
  }
  free(chunk_psi);
  free(chunk_nuc);
  free(back);
}

#endif
