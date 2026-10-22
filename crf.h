#ifndef CRF_H
#define CRF_H
#include "anno.h"
#include <stdint.h>

#define CRF_RESIDUAL_WEIGHT 3.0
#define HMM_PSEUDOCOUNT 1.0

enum {
  CRF_WINDOW_SIZE = 65536,
  CRF_SELF_HALO = GFF_CWT_RADIUS + 1,
  STATES = 13, EMIT_STATES = 3, HIDDEN_STATES = 86,
  TRANSITION_FEATURES = 3 + GFF_BOUND_KINDS,
  K3 = 65, KNUC = 5,
  BOUNDARY_PWM_FROM = -24, BOUNDARY_PWM_COUNT = 33,
  OFF_E3 = 0,
  OFF_PWM = OFF_E3 + EMIT_STATES * K3,
  OFF_TRANSITION = OFF_PWM + GFF_BOUND_KINDS * BOUNDARY_PWM_COUNT * KNUC,
  WEIGHTS = OFF_TRANSITION + TRANSITION_FEATURES
};

void crf_init(void);
void crf_set_gff_cwt(const GffCwtModel *model);
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
  uint8_t nuc;
  uint8_t k3[2];
  uint8_t boundary_nuc[2][BOUNDARY_PWM_COUNT];
} Site;

static const int hidden_offsets[STATES] = {0, 1, 4, 8, 10, 15, 25, 40, 42, 46, 51, 66, 81};
static const int hidden_counts[STATES] = {1, 3, 4, 2, 5, 10, 15, 2, 4, 5, 15, 15, 5};
static int trans_a[STATES * STATES], trans_b[STATES * STATES], ntrans;
static int edge_rev[STATES * STATES], edge_kind[STATES * STATES], edge_trans[STATES * STATES];
static uint8_t hidden_state[HIDDEN_STATES], hidden_context[HIDDEN_STATES], hidden_sub[HIDDEN_STATES];
typedef struct { uint8_t previous, current; uint16_t transition; } HiddenEdge;
static HiddenEdge hidden_edges[KNUC][HIDDEN_STATES * 4];
static int hidden_edge_count[KNUC];
static GffCwtModel crf_gff_cwt;
static double crf_bound_energy[GFF_BOUND_KINDS];
static int crf_gff_cwt_set;

void crf_set_gff_cwt(const GffCwtModel *model) {
  if (model) {
    crf_gff_cwt = *model;
    for (int k = 0; k < GFF_BOUND_KINDS; k++)
      crf_bound_energy[k] = gff_cwt_bound_energy(model, k);
    crf_gff_cwt_set = 1;
  } else {
    memset(&crf_gff_cwt, 0, sizeof(crf_gff_cwt));
    memset(crf_bound_energy, 0, sizeof(crf_bound_energy));
    crf_gff_cwt_set = 0;
  }
}

static int kmer(const char *seq, int length, int from, int size, int reverse) {
  int code = 0;
  for (int offset = 0; offset < size; offset++) {
    int position = from + (reverse ? size - 1 - offset : offset);
    int base = position >= 0 && position < length ? base_index(seq[position]) : -1;
    if (base < 0) return 1 << (2 * size);
    code = code * 4 + (reverse ? 3 - base : base);
  }
  return code;
}

static Site site_at(const Contig *c, int position) {
  int base = position >= 0 && position < c->length ? base_index(c->seq[position]) : -1;
  Site s = {.nuc = (uint8_t)(base < 0 ? 4 : base)};
  for (int rev = 0; rev < 2; rev++) {
    s.k3[rev] = (uint8_t)kmer(c->seq, c->length, position - 1, 3, rev);
    for (int offset = 0; offset < BOUNDARY_PWM_COUNT; offset++) {
      int disp = BOUNDARY_PWM_FROM + offset;
      int g = rev ? position - disp - 1 : position + disp;
      int b = g >= 0 && g < c->length ? base_index(c->seq[g]) : -1;
      s.boundary_nuc[rev][offset] = (uint8_t)(b < 0 ? 4 : rev ? 3 - b : b);
    }
  }
  return s;
}

static inline int emission_index(const Site *site, int state) {
  return OFF_E3 + emit_region(state) * K3 + site->k3[state >= 7];
}

static inline int boundary_pwm_index(const Site *site, int reverse, int kind, int offset) {
  return OFF_PWM + (kind * BOUNDARY_PWM_COUNT + offset) * KNUC +
         site->boundary_nuc[reverse][offset];
}

/* Returns 1 if displacement disp is already deterministically fixed by the grammar. */
static inline int pwm_offset_locked(int kind, int disp) {
  if (kind == 0) return disp == 0 || disp == 1;
  if (kind == 1) return disp == -2 || disp == -1;
  if (kind == 2) return disp >= 0 && disp <= 2;
  return disp >= -3 && disp <= -1;
}

static double boundary_score(const double *w, const Site *site, int rev, int kind) {
  double score = 0.0;
  for (int offset = 0; offset < BOUNDARY_PWM_COUNT; offset++)
    score += w[boundary_pwm_index(site, rev, kind, offset)];
  return score;
}

static inline int state_transition_index(int a, int b) {
  int kind = boundary_kind(a, b);
  if (kind >= 0) return OFF_TRANSITION + 3 + kind;
  return OFF_TRANSITION + emit_region(a);
}

static int allowed(int a, int b) {
  if (a == 0) return b == 0 || b == 1 || b == 9;
  if (a <= 3) return b == 1 + a % 3 || b == 4 + a % 3 || (a == 3 && b == 0);
  if (a <= 6) return b == a || b == a - 3;
  if (a <= 9) return b == 7 + (a + 1) % 3 || b == 10 + (a + 1) % 3 || (a == 7 && b == 0);
  return b == a || b == a - 3;
}

static int initial_allowed(int base, int hidden) {
  return !hidden ||
         (base == 0 && hidden == hidden_offsets[1] + 2) ||
         (base == 3 && hidden == hidden_offsets[9] + 3) ||
         (base == 1 && hidden == hidden_offsets[9] + 4);
}

static int terminal_allowed(int hidden) {
  int state = hidden_state[hidden];
  return !state || ((state == 3 || state == 7) && hidden_context[hidden] == 1);
}

static int find_or_add_edge(int a, int b) {
  int trans = state_transition_index(a, b);
  for (int e = 0; e < ntrans; e++)
    if (trans_a[e] == a && trans_b[e] == b && edge_trans[e] == trans)
      return e;
  if (ntrans >= STATES * STATES) anno_fail("Too many CRF transition edges");
  int e = ntrans++;
  trans_a[e] = a;
  trans_b[e] = b;
  edge_rev[e] = boundary_rev(a, b);
  edge_kind[e] = boundary_kind(a, b);
  edge_trans[e] = trans;
  return e;
}

static void add_hidden_edge(int base, int previous, int current) {
  if (current < 0) return;
  if (hidden_edge_count[base] >= HIDDEN_STATES * 4) anno_fail("Too many grammar transitions");
  int a = hidden_state[previous], b = hidden_state[current];
  if (!allowed(a, b)) anno_fail("Invalid state transition: %d -> %d", a, b);
  int edge = find_or_add_edge(a, b);
  hidden_edges[base][hidden_edge_count[base]++] =
    (HiddenEdge){(uint8_t)previous, (uint8_t)current, (uint16_t)edge};
}

static void add_coding_edges_plus(int base, int previous, int target_state, int context) {
  if (base == 4) return;
  if (target_state == 1) {
    add_hidden_edge(base, previous, hidden_offsets[1] + (base == 3 ? 1 : 0));
    return;
  }
  if (target_state == 2) {
    if (context == 0)
      add_hidden_edge(base, previous, hidden_offsets[2] + 0);
    else if (context == 1)
      add_hidden_edge(base, previous, hidden_offsets[2] + (base == 0 ? 1 : base == 2 ? 2 : 0));
    else if (context == 2 && base == 3)
      add_hidden_edge(base, previous, hidden_offsets[2] + 3);
    return;
  }
  if (context == 0)
    add_hidden_edge(base, previous, hidden_offsets[3] + 0);
  else if (context == 1)
    add_hidden_edge(base, previous, hidden_offsets[3] + (base == 0 || base == 2 ? 1 : 0));
  else if (context == 2)
    add_hidden_edge(base, previous, hidden_offsets[3] + (base == 0 ? 1 : 0));
  else if (context == 3 && base == 2)
    add_hidden_edge(base, previous, hidden_offsets[3] + 0);
}

static void add_coding_edges_minus(int base, int previous, int target_state, int context) {
  if (base == 4) return;
  if (target_state == 9) {
    add_hidden_edge(base, previous, hidden_offsets[9] + (base == 3 ? 1 : base == 1 ? 2 : 0));
    return;
  }
  if (target_state == 8) {
    if (context == 0)
      add_hidden_edge(base, previous, hidden_offsets[8] + 0);
    else if (context == 1)
      add_hidden_edge(base, previous, hidden_offsets[8] + (base == 1 || base == 3 ? 1 : 0));
    else if (context == 2)
      add_hidden_edge(base, previous, hidden_offsets[8] + (base == 3 ? 1 : base == 0 ? 2 : 0));
    else if (context == 3 && (base == 1 || base == 3))
      add_hidden_edge(base, previous, hidden_offsets[8] + 3);
    else if (context == 4 && base == 3)
      add_hidden_edge(base, previous, hidden_offsets[8] + 3);
    return;
  }
  if (context == 0)
    add_hidden_edge(base, previous, hidden_offsets[7] + 0);
  else if (context == 1 && base != 0)
    add_hidden_edge(base, previous, hidden_offsets[7] + 0);
  else if (context == 2)
    add_hidden_edge(base, previous, hidden_offsets[7] + (base == 3 ? 1 : 0));
  else if (context == 3 && base == 0)
    add_hidden_edge(base, previous, hidden_offsets[7] + 0);
}

/* Deterministic ORF and GT-AG splice grammar over 13 observed states. */
static void add_grammar_edges(int base, int previous) {
  int state = hidden_state[previous], context = hidden_context[previous];
  int sub = hidden_sub[previous], minus = state >= 7;
  if (!state) {
    add_hidden_edge(base, previous, 0);
    if (base == 0) add_hidden_edge(base, previous, hidden_offsets[1] + 2);
    if (base == 3) add_hidden_edge(base, previous, hidden_offsets[9] + 3);
    if (base == 1) add_hidden_edge(base, previous, hidden_offsets[9] + 4);
    return;
  }
  if (is_cds(state)) {
    if (state == 3 && context == 1) {
      add_hidden_edge(base, previous, 0);
      return;
    }
    if (state == 7 && context == 1)
      add_hidden_edge(base, previous, 0);
    int next_c = minus ? 7 + (state + 1) % 3 : 1 + state % 3;
    int next_i = next_c + 3, reset = next_c == (minus ? 9 : 1);
    if (minus) add_coding_edges_minus(base, previous, next_c, reset ? 0 : context);
    else add_coding_edges_plus(base, previous, next_c, reset ? 0 : context);
    int nctx_i = hidden_counts[next_i] / 5;
    int i_ctx = reset ? 0 : context;
    if (i_ctx < nctx_i && base == (minus ? 1 : 2))
      add_hidden_edge(base, previous, hidden_offsets[next_i] + i_ctx);
    return;
  }
  int nctx = hidden_counts[state] / 5, off = hidden_offsets[state] + context;
  if (sub == 0) {
    if (base == 3) add_hidden_edge(base, previous, off + 1 * nctx);
  } else if (sub == 1) {
    if (base < 4) add_hidden_edge(base, previous, off + 2 * nctx);
  } else if (sub == 2) {
    if (base < 4) add_hidden_edge(base, previous, off + 2 * nctx);
    if (base == 0) add_hidden_edge(base, previous, off + 3 * nctx);
  } else if (sub == 3) {
    if (base == (minus ? 1 : 2)) add_hidden_edge(base, previous, off + 4 * nctx);
  } else if (minus) {
    add_coding_edges_minus(base, previous, state - 3, context);
  } else {
    add_coding_edges_plus(base, previous, state - 3, context);
  }
}

void crf_init(void) {
  ntrans = 0;
  memset(hidden_state, 0, sizeof(hidden_state));
  memset(hidden_context, 0, sizeof(hidden_context));
  memset(hidden_sub, 0, sizeof(hidden_sub));

  for (int state = 1; state < STATES; state++) {
    int off = hidden_offsets[state], cnt = hidden_counts[state];
    int nsub = is_intron(state) ? 5 : 1, nctx = cnt / nsub;
    for (int sub = 0; sub < nsub; sub++)
      for (int c = 0; c < nctx; c++) {
        int h = off + sub * nctx + c;
        hidden_state[h] = (uint8_t)state;
        hidden_context[h] = (uint8_t)c;
        hidden_sub[h] = (uint8_t)sub;
      }
  }

  memset(hidden_edge_count, 0, sizeof(hidden_edge_count));
  for (int base = 0; base < KNUC; base++)
    for (int previous = 0; previous < HIDDEN_STATES; previous++)
      add_grammar_edges(base, previous);
}

static Site site_scores(const double *w, const Contig *c, int position,
                        const double *cwt_bonus, double *emissions, double *psi) {
  Site site = site_at(c, position);
  for (int state = 0; state < STATES; state++) {
    double score = w[emission_index(&site, state)];
    if (!isfinite(score))
      anno_fail("Non-finite emission at position %d", position + 1);
    emissions[state] = score;
  }
  double b_score[2][GFF_BOUND_KINDS] = {{0}};
  for (int rev = 0; rev < 2; rev++)
    for (int kind = 0; kind < GFF_BOUND_KINDS; kind++) {
      b_score[rev][kind] = boundary_score(w, &site, rev, kind);
      if (cwt_bonus)
        b_score[rev][kind] += cwt_bonus[rev * GFF_BOUND_KINDS + kind];
    }
  for (int edge = 0; edge < ntrans; edge++) {
    int rev = edge_rev[edge], kind = edge_kind[edge];
    double score = w[edge_trans[edge]] + (kind >= 0 ? b_score[rev][kind] : 0.0);
    double p = score + emissions[trans_b[edge]];
    if (!isfinite(p))
      anno_fail("Non-finite transition at position %d", position + 1);
    psi[edge] = p;
  }
  return site;
}

typedef struct {
  Site *sites;
  double *psi, *log_norm, *log_alpha;
  double *bound_marg, *pred_cwt, *res_cwt, *grad_bound, *cwt_bonus;
  int capacity;
} Work;

static void work_free(Work *k) {
  free(k->sites), free(k->psi), free(k->log_norm), free(k->log_alpha);
  free(k->bound_marg), free(k->pred_cwt), free(k->res_cwt), free(k->grad_bound), free(k->cwt_bonus);
  memset(k, 0, sizeof(*k));
}

static void work_ensure(Work *k, int length) {
  if (length <= k->capacity) return;
  work_free(k);
  k->sites = anno_alloc((size_t)length, sizeof(Site));
  k->psi = anno_alloc((size_t)length, (size_t)ntrans * sizeof(double));
  k->log_norm = anno_alloc((size_t)length, sizeof(double));
  k->log_alpha = anno_alloc((size_t)length, HIDDEN_STATES * sizeof(double));
  k->bound_marg = anno_alloc((size_t)length, (2 * GFF_BOUND_KINDS) * sizeof(double));
  k->pred_cwt = anno_alloc((size_t)length, CWT_CHANNELS * sizeof(double));
  k->res_cwt = anno_alloc((size_t)length, CWT_CHANNELS * sizeof(double));
  k->grad_bound = anno_alloc((size_t)length, (2 * GFF_BOUND_KINDS) * sizeof(double));
  k->cwt_bonus = anno_alloc((size_t)length, (2 * GFF_BOUND_KINDS) * sizeof(double));
  k->capacity = length;
}

static double normalize_log_alpha(double *log_alpha, int position) {
  double max_score = -INFINITY;
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
    if (log_alpha[hidden] > max_score) max_score = log_alpha[hidden];
  if (!isfinite(max_score)) anno_fail("Invalid partition at position %d", position + 1);
  double sum = 0;
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
    if (log_alpha[hidden] > -INFINITY) sum += exp(log_alpha[hidden] - max_score);
  double log_norm = max_score + log(sum);
  if (!isfinite(log_norm)) anno_fail("Invalid partition at position %d", position + 1);
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) log_alpha[hidden] -= log_norm;
  return log_norm;
}

static inline double edge_posterior(const double *prev_L, const double *L,
                                    const double *psi, double log_norm,
                                    const HiddenEdge *step) {
  return (prev_L[step->previous] > -INFINITY && L[step->current] > -INFINITY)
    ? exp(prev_L[step->previous] + psi[step->transition] - log_norm - L[step->current])
    : 0.0;
}

static void backward_step(const Work *k, int position, double *B, double *edge_accum) {
  const Site *site = k->sites + position;
  const double *L = k->log_alpha + (size_t)position * HIDDEN_STATES;
  const double *prev_L = L - HIDDEN_STATES;
  const double *psi = k->psi + (size_t)position * ntrans;
  double log_norm = k->log_norm[position];
  double B_prev[HIDDEN_STATES] = {0};
  const HiddenEdge *edges = hidden_edges[site->nuc];
  int edge_count = hidden_edge_count[site->nuc];
  for (int edge = 0; edge < edge_count; edge++) {
    const HiddenEdge *step = edges + edge;
    int u = step->previous, v = step->current, tr = step->transition;
    if (!(B[v] > 0.0)) continue;
    double q = edge_posterior(prev_L, L, psi, log_norm, step);
    if (q == 0.0) continue;
    double edge_prob = B[v] * q;
    B_prev[u] += edge_prob;
    if (edge_accum) edge_accum[tr] += edge_prob;
  }
  double sum_b = 0;
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) sum_b += B_prev[hidden];
  if (!(sum_b > 0.0) || !isfinite(sum_b))
    anno_fail("Invalid posterior at position %d", position + 1);
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) B[hidden] = B_prev[hidden] / sum_b;
}

typedef struct {
  double emit_counts[EMIT_STATES][K3];
  double pwm_counts[GFF_BOUND_KINDS][BOUNDARY_PWM_COUNT][KNUC];
  double trans_counts[TRANSITION_FEATURES];
  double nuc_counts[4];
} HmmStats;

static void forward_backward_window(const double *w, const Contig *c, Work *k,
                                    int start, int length, const double *cwt_bonus) {
  for (int position = 0; position < length; position++) {
    double emissions[STATES];
    double *psi = k->psi + (size_t)position * ntrans;
    const double *bonus = cwt_bonus ? cwt_bonus + (size_t)position * (2 * GFF_BOUND_KINDS) : NULL;
    Site site = site_scores(w, c, start + position, bonus, emissions, psi);
    k->sites[position] = site;

    double *L = k->log_alpha + (size_t)position * HIDDEN_STATES;
    if (!position) {
      for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
        L[hidden] = emissions[hidden_state[hidden]];
    } else {
      const double *prev_L = L - HIDDEN_STATES;
      const HiddenEdge *edges = hidden_edges[site.nuc];
      int edge_count = hidden_edge_count[site.nuc];
      double incoming_sum[HIDDEN_STATES] = {0};
      for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) L[hidden] = -INFINITY;
      for (int edge = 0; edge < edge_count; edge++) {
        const HiddenEdge *step = edges + edge;
        if (!(prev_L[step->previous] > -INFINITY)) continue;
        double score = prev_L[step->previous] + psi[step->transition];
        if (score > L[step->current]) L[step->current] = score;
      }
      for (int edge = 0; edge < edge_count; edge++) {
        const HiddenEdge *step = edges + edge;
        if (!(prev_L[step->previous] > -INFINITY)) continue;
        incoming_sum[step->current] +=
          exp(prev_L[step->previous] + psi[step->transition] - L[step->current]);
      }
      for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
        if (L[hidden] > -INFINITY)
          L[hidden] += log(incoming_sum[hidden]);
    }
    k->log_norm[position] = normalize_log_alpha(L, position);
  }
}

/*
 * Evaluates the expected CWT residual loss:
 *   L_CWT = 0.5 * ||W_FASTA - W_hat(xi)||^2
 *         + 0.5 * sum_{t,s,k} E_k * xi_{t,s,k} * (1 - xi_{t,s,k})
 * and converts the negative gradient -dL_CWT/d(xi_{t,s,k}) into a dimensionless
 * log-likelihood bonus for the E-step posterior tilt.
 */
static double evaluate_cwt_residual_and_bonus(const Contig *c, Work *k,
                                              int start, int length,
                                              int residual_start, int residual_length,
                                              double res_weight) {
  gff_cwt_predict_soft(k->bound_marg, length, &crf_gff_cwt, k->pred_cwt);
  double r_loss = 0.0;
  for (int position = 0; position < length; position++) {
    int g = start + position;
    double *res = k->res_cwt + (size_t)position * CWT_CHANNELS;
    if (!c->cwt || base_index(c->seq[g]) < 0) {
      memset(res, 0, CWT_CHANNELS * sizeof(*res));
      continue;
    }
    const float *obs = c->cwt + (size_t)g * CWT_CHANNELS;
    const double *pred = k->pred_cwt + (size_t)position * CWT_CHANNELS;
    for (int ch = 0; ch < CWT_CHANNELS; ch++) {
      double diff = (double)obs[ch] - pred[ch];
      res[ch] = diff;
      if (position >= residual_start && position < residual_start + residual_length)
        r_loss += 0.5 * diff * diff;
    }
    if (position && position >= residual_start && position < residual_start + residual_length) {
      const double *bm = k->bound_marg + (size_t)position * (2 * GFF_BOUND_KINDS);
      for (int rev = 0; rev < 2; rev++)
        for (int kind = 0; kind < GFF_BOUND_KINDS; kind++) {
          double xi = bm[rev * GFF_BOUND_KINDS + kind];
          r_loss += 0.5 * crf_bound_energy[kind] * xi * (1.0 - xi);
        }
    }
  }
  gff_cwt_residual_adjoint(k->res_cwt, length, &crf_gff_cwt, k->grad_bound);
  memset(k->cwt_bonus, 0, (size_t)length * (2 * GFF_BOUND_KINDS) * sizeof(*k->cwt_bonus));
  for (int position = 1; position < length; position++) {
    const double *bm = k->bound_marg + (size_t)position * (2 * GFF_BOUND_KINDS);
    const double *gb = k->grad_bound + (size_t)position * (2 * GFF_BOUND_KINDS);
    double *bonus = k->cwt_bonus + (size_t)position * (2 * GFF_BOUND_KINDS);
    for (int rev = 0; rev < 2; rev++)
      for (int kind = 0; kind < GFF_BOUND_KINDS; kind++) {
        double ek = crf_bound_energy[kind];
        if (ek <= 0.0) continue;
        double xi = bm[rev * GFF_BOUND_KINDS + kind];
        double neg_grad = gb[rev * GFF_BOUND_KINDS + kind] - 0.5 * ek * (1.0 - 2.0 * xi);
        bonus[rev * GFF_BOUND_KINDS + kind] = res_weight * (neg_grad / ek);
      }
  }
  return r_loss;
}

/*
 * Performs one Gauss-Newton / Proximal-EM step on a window:
 * 1. Forward-backward under current HMM parameters w -> boundary marginals xi(w).
 * 2. Compute CWT residual R = W_FASTA - W_hat(xi) and boundary log-bonus b_CWT.
 * 3. Forward-backward under tilted score w + b_CWT -> accumulate M-step counts.
 */
static double window_em_step(const double *w, const Contig *c, Work *k,
                             HmmStats *stats, double res_weight,
                             int start, int length,
                             int residual_start, int residual_length) {
  if (start < 0 || length <= 0 || start > c->length - length ||
      residual_start < 0 || residual_length < 0 ||
      residual_start > length - residual_length)
    anno_fail("Invalid window span");
  work_ensure(k, length);

  forward_backward_window(w, c, k, start, length, NULL);
  const double *last_L = k->log_alpha + (size_t)(length - 1) * HIDDEN_STATES;
  double B[HIDDEN_STATES];
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) B[hidden] = exp(last_L[hidden]);
  memset(k->bound_marg, 0, (size_t)length * (2 * GFF_BOUND_KINDS) * sizeof(*k->bound_marg));
  for (int position = length - 1; position >= 1; position--) {
    double edge_m[STATES * STATES] = {0};
    backward_step(k, position, B, edge_m);
    double *bm = k->bound_marg + (size_t)position * (2 * GFF_BOUND_KINDS);
    for (int edge = 0; edge < ntrans; edge++) {
      int kind = edge_kind[edge];
      if (kind >= 0) bm[edge_rev[edge] * GFF_BOUND_KINDS + kind] += edge_m[edge];
    }
  }

  double r_loss = evaluate_cwt_residual_and_bonus(
    c, k, start, length, residual_start, residual_length, res_weight);
  if (!stats) return r_loss;

  forward_backward_window(w, c, k, start, length, k->cwt_bonus);
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) B[hidden] = exp(last_L[hidden]);
  for (int position = length - 1; position >= 0; position--) {
    const Site *site = k->sites + position;
    int in_core = position >= residual_start && position < residual_start + residual_length;
    if (in_core) {
      if (site->nuc < 4) {
        stats->nuc_counts[site->nuc] += 0.5;
        stats->nuc_counts[3 - site->nuc] += 0.5;
      }
      for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) {
        int st = hidden_state[hidden], r = emit_region(st);
        if (r == 0) {
          stats->emit_counts[0][site->k3[0]] += 0.5 * B[hidden];
          stats->emit_counts[0][site->k3[1]] += 0.5 * B[hidden];
        } else {
          stats->emit_counts[r][site->k3[st >= 7]] += B[hidden];
        }
      }
    }

    if (position) {
      double edge_m[STATES * STATES] = {0};
      backward_step(k, position, B, edge_m);
      if (in_core) {
        for (int edge = 0; edge < ntrans; edge++) {
          double m = edge_m[edge];
          if (m == 0.0) continue;
          int tr_idx = edge_trans[edge] - OFF_TRANSITION;
          stats->trans_counts[tr_idx] += m;
          int kind = edge_kind[edge];
          if (kind >= 0) {
            int rev = edge_rev[edge];
            for (int d = 0; d < BOUNDARY_PWM_COUNT; d++)
              stats->pwm_counts[kind][d][site->boundary_nuc[rev][d]] += m;
          }
        }
      }
    }
  }
  return r_loss;
}

typedef struct {
  Contig *contig;
  int start, length, residual_start, residual_length;
} CrfWindow;

static size_t make_windows(Contig *contigs, int count, CrfWindow **windows) {
  size_t window_count = 0, capacity = 0;
  for (int i = 0; i < count; i++) {
    Contig *contig = contigs + i;
    int limit = contig->length;
    for (int core_start = 0; core_start < limit;) {
      int core_end = limit - core_start > CRF_WINDOW_SIZE
        ? core_start + CRF_WINDOW_SIZE : limit;
      int start = core_start > CRF_SELF_HALO ? core_start - CRF_SELF_HALO : 0;
      int end = limit - core_end > CRF_SELF_HALO ? core_end + CRF_SELF_HALO : limit;
      if (window_count == capacity) {
        if (capacity > SIZE_MAX / 2 / sizeof(**windows))
          anno_fail("Too many windows");
        size_t next_capacity = capacity ? capacity * 2 : 16;
        CrfWindow *grown = realloc(*windows, next_capacity * sizeof(**windows));
        if (!grown) anno_fail("Out of memory");
        *windows = grown;
        capacity = next_capacity;
      }
      (*windows)[window_count++] = (CrfWindow){
        contig, start, end - start, core_start - start, core_end - core_start
      };
      core_start = core_end;
    }
  }
  return window_count;
}

typedef struct {
  CrfWindow *windows;
  size_t count;
  atomic_size_t next_window;
  int threads;
} CrfTrain;

typedef struct {
  CrfTrain *train;
  const double *w;
  double res_loss;
  size_t positions;
  HmmStats stats;
  Work work;
} CrfWorker;

static void *crf_worker(void *arg) {
  CrfWorker *worker = arg;
  CrfTrain *train = worker->train;
  worker->res_loss = 0;
  worker->positions = 0;
  memset(&worker->stats, 0, sizeof(worker->stats));
  for (;;) {
    size_t slot = atomic_fetch_add_explicit(&train->next_window, 1, memory_order_relaxed);
    if (slot >= train->count) break;
    const CrfWindow *window = train->windows + slot;
    double res = window_em_step(
      worker->w, window->contig, &worker->work, &worker->stats, CRF_RESIDUAL_WEIGHT,
      window->start, window->length, window->residual_start, window->residual_length);
    if (!isfinite(res)) anno_fail("Non-finite training residual");
    worker->res_loss += res;
    worker->positions += (size_t)window->residual_length;
  }
  return NULL;
}

static void hmm_m_step(const HmmStats *stats, double *w) {
  double bg_nuc[4] = {0.25, 0.25, 0.25, 0.25};
  double nuc_total = 0.0;
  for (int b = 0; b < 4; b++) nuc_total += stats->nuc_counts[b] + HMM_PSEUDOCOUNT;
  if (nuc_total > 0.0)
    for (int b = 0; b < 4; b++)
      bg_nuc[b] = (stats->nuc_counts[b] + HMM_PSEUDOCOUNT) / nuc_total;

  double emit_prob[EMIT_STATES][K3];
  for (int r = 0; r < EMIT_STATES; r++) {
    double total = 0.0;
    for (int k = 0; k < 64; k++) total += stats->emit_counts[r][k] + HMM_PSEUDOCOUNT;
    for (int k = 0; k < 64; k++)
      emit_prob[r][k] = (stats->emit_counts[r][k] + HMM_PSEUDOCOUNT) / total;
    emit_prob[r][64] = 1.0;
  }
  for (int k = 0; k < K3; k++) {
    w[OFF_E3 + 0 * K3 + k] = 0.0;
    w[OFF_E3 + 1 * K3 + k] = k < 64 ? log(emit_prob[1][k] / emit_prob[0][k]) : 0.0;
    w[OFF_E3 + 2 * K3 + k] = k < 64 ? log(emit_prob[2][k] / emit_prob[0][k]) : 0.0;
  }

  for (int kind = 0; kind < GFF_BOUND_KINDS; kind++)
    for (int d = 0; d < BOUNDARY_PWM_COUNT; d++) {
      int disp = BOUNDARY_PWM_FROM + d;
      int base_idx = OFF_PWM + (kind * BOUNDARY_PWM_COUNT + d) * KNUC;
      if (pwm_offset_locked(kind, disp)) {
        for (int b = 0; b < KNUC; b++) w[base_idx + b] = 0.0;
        continue;
      }
      double total = 0.0;
      for (int b = 0; b < 4; b++) total += stats->pwm_counts[kind][d][b] + HMM_PSEUDOCOUNT;
      for (int b = 0; b < 4; b++) {
        double p = (stats->pwm_counts[kind][d][b] + HMM_PSEUDOCOUNT) / total;
        w[base_idx + b] = log(p / bg_nuc[b]);
      }
      w[base_idx + 4] = 0.0;
    }

  double stay_I = stats->trans_counts[0] + HMM_PSEUDOCOUNT;
  double exit_I_start = stats->trans_counts[3 + 2] + HMM_PSEUDOCOUNT;
  double tot_I = stay_I + exit_I_start;

  double stay_E = stats->trans_counts[1] + HMM_PSEUDOCOUNT;
  double exit_E_don = stats->trans_counts[3 + 0] + HMM_PSEUDOCOUNT;
  double exit_E_stop = stats->trans_counts[3 + 3] + HMM_PSEUDOCOUNT;
  double tot_E = stay_E + exit_E_don + exit_E_stop;

  double stay_N = stats->trans_counts[2] + HMM_PSEUDOCOUNT;
  double exit_N_acc = stats->trans_counts[3 + 1] + HMM_PSEUDOCOUNT;
  double tot_N = stay_N + exit_N_acc;

  double p_stay_I = stay_I / tot_I;
  w[OFF_TRANSITION + 0] = 0.0;
  w[OFF_TRANSITION + 1] = log((stay_E / tot_E) / p_stay_I);
  w[OFF_TRANSITION + 2] = log((stay_N / tot_N) / p_stay_I);
  w[OFF_TRANSITION + 3 + 0] = log((exit_E_don / tot_E) / p_stay_I);
  w[OFF_TRANSITION + 3 + 1] = log((exit_N_acc / tot_N) / p_stay_I);
  w[OFF_TRANSITION + 3 + 2] = log((exit_I_start / tot_I) / p_stay_I);
  w[OFF_TRANSITION + 3 + 3] = 0.0;
}

void crf_train(Contig *contigs, int count, double *w, int iterations, int threads) {
  if (count < 1 || iterations < 1)
    anno_fail("Training requires instances and positive iterations");
  if (threads < 1) anno_fail("Training requires at least one worker thread");
  if (!crf_gff_cwt_set)
    anno_fail("Self-training requires an embedded GFF->CWT model");
  for (int i = 0; i < count; i++)
    if (contigs[i].label)
      anno_fail("Self-training input must not contain GFF-derived labels");

  CrfWindow *windows = NULL;
  size_t window_count = make_windows(contigs, count, &windows);
  if (!window_count) anno_fail("No valid training windows");
  CrfTrain train = {.windows = windows, .count = window_count, .threads = threads};
  fprintf(stderr, "Using %d worker threads for %zu training windows\n",
          train.threads, train.count);
  CrfWorker *workers = anno_alloc((size_t)threads, sizeof(*workers));
  pthread_t *tids = anno_alloc((size_t)threads, sizeof(*tids));
  for (int t = 0; t < threads; t++) workers[t].train = &train;

  for (int iter = 1; iter <= iterations; iter++) {
    atomic_init(&train.next_window, 0);
    for (int t = 0; t < threads; t++) {
      workers[t].w = w;
      if (threads > 1 && pthread_create(tids + t, NULL, crf_worker, workers + t))
        anno_fail("Cannot spawn training thread");
    }
    HmmStats total_stats = {0};
    double res_loss = 0.0;
    size_t positions = 0;
    for (int t = 0; t < threads; t++) {
      if (threads > 1 && pthread_join(tids[t], NULL))
        anno_fail("Cannot join training thread");
      if (threads == 1) crf_worker(workers + t);
      res_loss += workers[t].res_loss;
      positions += workers[t].positions;
      for (int r = 0; r < EMIT_STATES; r++)
        for (int k = 0; k < K3; k++)
          total_stats.emit_counts[r][k] += workers[t].stats.emit_counts[r][k];
      for (int k = 0; k < GFF_BOUND_KINDS; k++)
        for (int d = 0; d < BOUNDARY_PWM_COUNT; d++)
          for (int b = 0; b < KNUC; b++)
            total_stats.pwm_counts[k][d][b] += workers[t].stats.pwm_counts[k][d][b];
      for (int f = 0; f < TRANSITION_FEATURES; f++)
        total_stats.trans_counts[f] += workers[t].stats.trans_counts[f];
      for (int b = 0; b < 4; b++)
        total_stats.nuc_counts[b] += workers[t].stats.nuc_counts[b];
    }

    double next_w[WEIGHTS];
    memcpy(next_w, w, sizeof(next_w));
    hmm_m_step(&total_stats, next_w);
    double delta_norm = 0.0;
    for (int f = 0; f < WEIGHTS; f++) {
      double step = 0.5 * (next_w[f] - w[f]);
      w[f] += step;
      delta_norm = hypot(delta_norm, step);
    }
    fprintf(stderr, "Iteration %d: %zu windows, cwt-residual/bp %.8g, param-step %.3g\n",
            iter, train.count, res_loss / positions, delta_norm);
    if (iter > 1 && delta_norm < 1e-3) {
      fprintf(stderr, "EM self-training converged after %d iterations\n", iter);
      break;
    }
  }

  for (int t = 0; t < threads; t++) work_free(&workers[t].work);
  free(workers), free(tids), free(windows);
}

void crf_decode(const double *w, const Contig *c, uint8_t *path) {
  int length = c->length;
  uint8_t *back = anno_alloc((size_t)length, HIDDEN_STATES * sizeof(*back));
  double previous[HIDDEN_STATES] = {0}, current[HIDDEN_STATES];
  double emissions[STATES], transitions[STATES * STATES];

  for (int pos = 0; pos < length; pos++) {
    Site site = site_scores(w, c, pos, NULL, emissions, transitions);
    uint8_t *bptr = back + (size_t)pos * HIDDEN_STATES;
    for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) {
      current[hidden] = pos == 0 && initial_allowed(site.nuc, hidden)
        ? emissions[hidden_state[hidden]] : -INFINITY;
      bptr[hidden] = 0;
    }
    if (pos) {
      const HiddenEdge *edges = hidden_edges[site.nuc];
      int edge_count = hidden_edge_count[site.nuc];
      for (int edge = 0; edge < edge_count; edge++) {
        const HiddenEdge *step = edges + edge;
        double prev_score = previous[step->previous];
        if (prev_score == -INFINITY) continue;
        double score = prev_score + transitions[step->transition];
        if (score > current[step->current]) {
          current[step->current] = score;
          bptr[step->current] = step->previous;
        }
      }
    }
    double maximum = -INFINITY;
    for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
      maximum = fmax(maximum, current[hidden]);
    if (!isfinite(maximum)) anno_fail("No finite prediction path");
    for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
      previous[hidden] = current[hidden] - maximum;
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
