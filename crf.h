#ifndef CRF_H
#define CRF_H
#include "anno.h"
#include "lbfgs.h"
#include <stdint.h>

#define CRF_L2_C2 1.0
#define CRF_RESIDUAL_WEIGHT 0.05

static const int boundary_kmer_offsets[] = {-4, -3, -2, -1, 0, 1, 2};
enum {
  CRF_WINDOW_SIZE = 65536,
  CRF_SELF_HALO = GFF_CWT_RADIUS + 1,
  STATES = 13, EMIT_STATES = 3, HIDDEN_STATES = 137, TRANSITION_FEATURES = 21,
  K2 = 17, K3 = 65, KNUC = 5,
  BOUNDARY_PWM_FROM = -24, BOUNDARY_PWM_COUNT = 33,
  BOUNDARY_KMER_COUNT = sizeof(boundary_kmer_offsets) / sizeof(boundary_kmer_offsets[0]),
  OFF_E3 = 0,
  OFF_SPLICE_START_STOP = OFF_E3 + EMIT_STATES * K3,
  OFF_PWM = OFF_SPLICE_START_STOP + 4 * BOUNDARY_KMER_COUNT * K2,
  OFF_TRANSITION = OFF_PWM + 4 * BOUNDARY_PWM_COUNT * KNUC,
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
  uint8_t boundary_k2[2][BOUNDARY_KMER_COUNT];
  uint8_t boundary_nuc[2][BOUNDARY_PWM_COUNT];
} Site;

static const int hidden_offsets[STATES] = {0, 1, 5, 11, 14, 19, 39, 69, 72, 76, 82, 102, 132};
static const int hidden_counts[STATES] = {1, 4, 6, 3, 5, 20, 30, 3, 4, 6, 20, 30, 5};
static int trans_a[STATES * STATES], trans_b[STATES * STATES], ntrans;
static int edge_rev[STATES * STATES], edge_kind[STATES * STATES], edge_trans[STATES * STATES];
static uint8_t hidden_state[HIDDEN_STATES], hidden_context[HIDDEN_STATES];
static uint8_t hidden_motif[HIDDEN_STATES], hidden_sub[HIDDEN_STATES];
typedef struct { uint8_t previous, current; uint16_t transition; } HiddenEdge;
static HiddenEdge hidden_edges[KNUC][HIDDEN_STATES * 4];
static int hidden_edge_count[KNUC];
static GffCwtModel crf_gff_cwt;
static int crf_gff_cwt_set;

void crf_set_gff_cwt(const GffCwtModel *model) {
  if (model) {
    crf_gff_cwt = *model;
    crf_gff_cwt_set = 1;
  } else {
    memset(&crf_gff_cwt, 0, sizeof(crf_gff_cwt));
    crf_gff_cwt_set = 0;
  }
}

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

static Site site_at(const Contig *c, int position) {
  int base = position >= 0 && position < c->length ? base_index(c->seq[position]) : -1;
  Site s = {.nuc = (uint8_t)(base < 0 ? 4 : base)};
  for (int rev = 0; rev < 2; rev++) {
    s.k3[rev] = (uint8_t)kmer(c->seq, c->length, position - 1, 3, rev);
    for (int offset = 0; offset < BOUNDARY_KMER_COUNT; offset++) {
      int disp = boundary_kmer_offsets[offset];
      s.boundary_k2[rev][offset] = (uint8_t)kmer(
        c->seq, c->length, rev ? position - disp - 2 : position + disp, 2, rev);
    }
    for (int offset = 0; offset < BOUNDARY_PWM_COUNT; offset++) {
      int disp = BOUNDARY_PWM_FROM + offset;
      int g = rev ? position - disp - 1 : position + disp;
      int b = g >= 0 && g < c->length ? base_index(c->seq[g]) : -1;
      s.boundary_nuc[rev][offset] = (uint8_t)(b < 0 ? 4 : rev ? 3 - b : b);
    }
  }
  return s;
}

static int motif_transition_index(int prev_motif, int cur_motif, int reverse) {
  static const int indices[7][7] = {
    { 0,  1, -1, -1, -1, -1, -1},
    {-1,  2,  3,  4,  5, -1, -1},
    {-1, -1,  6,  7,  8, -1, -1},
    {10, -1, -1,  9, 11, -1, -1},
    {-1, -1, -1, -1, 12, 13, 14},
    {-1, -1, -1, -1, -1, 15, 16},
    {-1, 18, 19, 20, -1, -1, 17}
  };
  if (reverse) {
    int swap = prev_motif;
    prev_motif = cur_motif;
    cur_motif = swap;
  }
  if (prev_motif < 0 || prev_motif >= 7 || cur_motif < 0 || cur_motif >= 7 ||
      indices[prev_motif][cur_motif] < 0)
    anno_fail("Invalid CRF transition feature: %d -> %d", prev_motif, cur_motif);
  return OFF_TRANSITION + indices[prev_motif][cur_motif];
}

static inline int emission_index(const Site *site, int state) {
  int es = state == 0 ? 0 : is_cds(state) ? 1 : 2;
  return OFF_E3 + es * K3 + site->k3[state >= 7];
}

static inline int boundary_kmer_index(const Site *site, int reverse, int kind, int offset) {
  return OFF_SPLICE_START_STOP + (kind * BOUNDARY_KMER_COUNT + offset) * K2 +
         site->boundary_k2[reverse][offset];
}

static inline int boundary_pwm_index(const Site *site, int reverse, int kind, int offset) {
  return OFF_PWM + (kind * BOUNDARY_PWM_COUNT + offset) * KNUC +
         site->boundary_nuc[reverse][offset];
}

static double boundary_score(const double *w, const Site *site, int rev, int kind) {
  double score = 0.0;
  for (int offset = 0; offset < BOUNDARY_KMER_COUNT; offset++)
    score += w[boundary_kmer_index(site, rev, kind, offset)];
  for (int offset = 0; offset < BOUNDARY_PWM_COUNT; offset++)
    score += w[boundary_pwm_index(site, rev, kind, offset)];
  return score;
}

static void boundary_accumulate(double *grad, const Site *site, int rev, int kind, double coeff) {
  for (int offset = 0; offset < BOUNDARY_KMER_COUNT; offset++)
    grad[boundary_kmer_index(site, rev, kind, offset)] += coeff;
  for (int offset = 0; offset < BOUNDARY_PWM_COUNT; offset++)
    grad[boundary_pwm_index(site, rev, kind, offset)] += coeff;
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
         (base == 0 && hidden == hidden_offsets[1] + 0) ||
         (base == 3 && hidden == hidden_offsets[9] + 0) ||
         (base == 1 && hidden == hidden_offsets[9] + 1);
}

static int terminal_allowed(int hidden) {
  int state = hidden_state[hidden];
  return !state || ((state == 3 || state == 7) && hidden_context[hidden] == 2);
}

static int find_or_add_edge(int a, int b, int trans) {
  for (int e = 0; e < ntrans; e++)
    if (trans_a[e] == a && trans_b[e] == b && edge_trans[e] == trans)
      return e;
  if (ntrans >= STATES * STATES) anno_fail("Too many CRF transition edges");
  int e = ntrans++;
  trans_a[e] = a;
  trans_b[e] = b;
  edge_rev[e] = a >= 7 || b >= 7;
  edge_kind[e] = boundary_kind(a, b);
  edge_trans[e] = trans;
  return e;
}

static void add_hidden_edge(int base, int previous, int current) {
  if (current < 0) return;
  if (hidden_edge_count[base] >= HIDDEN_STATES * 4) anno_fail("Too many grammar transitions");
  int a = hidden_state[previous], b = hidden_state[current];
  if (!allowed(a, b)) anno_fail("Invalid CRF state transition: %d -> %d", a, b);
  int reverse = a >= 7 || b >= 7;
  int trans = motif_transition_index(hidden_motif[previous], hidden_motif[current], reverse);
  int edge = find_or_add_edge(a, b, trans);
  hidden_edges[base][hidden_edge_count[base]++] =
    (HiddenEdge){(uint8_t)previous, (uint8_t)current, (uint16_t)edge};
}

static void add_coding_edges_plus(int base, int previous, int target_state, int context) {
  if (base == 4) return;
  if (target_state == 1) {
    add_hidden_edge(base, previous, hidden_offsets[1] + (base == 3 ? 2 : 1));
    if (base == 3) add_hidden_edge(base, previous, hidden_offsets[1] + 3);
    return;
  }
  if (target_state == 2) {
    if (context == 0) {
      if (base == 3) add_hidden_edge(base, previous, hidden_offsets[2] + 0);
    } else if (context == 1) {
      add_hidden_edge(base, previous, hidden_offsets[2] + 1);
    } else if (context == 2) {
      int next = base == 0 ? 2 : base == 2 ? 3 : 1;
      add_hidden_edge(base, previous, hidden_offsets[2] + next);
    } else if (context == 3) {
      if (base == 0) add_hidden_edge(base, previous, hidden_offsets[2] + 4);
      else if (base == 2) add_hidden_edge(base, previous, hidden_offsets[2] + 5);
    }
    return;
  }
  if (context == 0) {
    if (base == 2) add_hidden_edge(base, previous, hidden_offsets[3] + 0);
  } else if (context == 1) {
    add_hidden_edge(base, previous, hidden_offsets[3] + 1);
  } else if (context == 2) {
    if (base != 0 && base != 2) add_hidden_edge(base, previous, hidden_offsets[3] + 1);
  } else if (context == 3) {
    if (base != 0) add_hidden_edge(base, previous, hidden_offsets[3] + 1);
  } else if (context == 4) {
    if (base == 0 || base == 2) add_hidden_edge(base, previous, hidden_offsets[3] + 2);
  } else if (context == 5) {
    if (base == 0) add_hidden_edge(base, previous, hidden_offsets[3] + 2);
  }
}

static void add_coding_edges_minus(int base, int previous, int target_state, int context) {
  if (base == 4) return;
  if (target_state == 9) {
    int next = base == 3 ? 3 : base == 1 ? 4 : 2;
    add_hidden_edge(base, previous, hidden_offsets[9] + next);
    if (base == 1) add_hidden_edge(base, previous, hidden_offsets[9] + 5);
    return;
  }
  if (target_state == 8) {
    if (context == 0) {
      if (base == 3 || base == 1) add_hidden_edge(base, previous, hidden_offsets[8] + 0);
    } else if (context == 1) {
      if (base == 3) add_hidden_edge(base, previous, hidden_offsets[8] + 0);
    } else if (context == 2) {
      add_hidden_edge(base, previous, hidden_offsets[8] + 1);
    } else if (context == 3) {
      add_hidden_edge(base, previous, hidden_offsets[8] + (base == 3 || base == 1 ? 2 : 1));
    } else if (context == 4) {
      add_hidden_edge(base, previous, hidden_offsets[8] + (base == 3 ? 2 : 1));
    } else if (context == 5) {
      if (base == 0) add_hidden_edge(base, previous, hidden_offsets[8] + 3);
    }
    return;
  }
  if (context == 0) {
    if (base == 0) add_hidden_edge(base, previous, hidden_offsets[7] + 0);
  } else if (context == 1) {
    add_hidden_edge(base, previous, hidden_offsets[7] + 1);
  } else if (context == 2) {
    if (base != 0) add_hidden_edge(base, previous, hidden_offsets[7] + 1);
  } else if (context == 3) {
    if (base == 3) add_hidden_edge(base, previous, hidden_offsets[7] + 2);
  }
}

/* Prefix and motif classes carry start/stop/splice recognition without changing 13 observed states. */
static void add_grammar_edges(int base, int previous) {
  int state = hidden_state[previous], context = hidden_context[previous];
  int sub = hidden_sub[previous], minus = state >= 7;
  if (!state) {
    add_hidden_edge(base, previous, 0);
    if (base == 0) add_hidden_edge(base, previous, hidden_offsets[1] + 0);
    if (base == 3) add_hidden_edge(base, previous, hidden_offsets[9] + 0);
    if (base == 1) add_hidden_edge(base, previous, hidden_offsets[9] + 1);
    return;
  }
  if (is_cds(state)) {
    if ((state == 3 || state == 7) && context == 2) {
      add_hidden_edge(base, previous, 0);
      return;
    }
    int next_c = minus ? 7 + (state + 1) % 3 : 1 + state % 3;
    int next_i = next_c + 3, reset = next_c == (minus ? 9 : 1);
    if (minus) add_coding_edges_minus(base, previous, next_c, reset ? 0 : context);
    else add_coding_edges_plus(base, previous, next_c, reset ? 0 : context);
    add_hidden_edge(base, previous, hidden_offsets[next_i] + (reset ? 0 : context));
    return;
  }
  int nctx = hidden_counts[state] / 5, off = hidden_offsets[state] + context;
  if (sub == 0) {
    add_hidden_edge(base, previous, off + 1 * nctx);
    if (minus) add_hidden_edge(base, previous, off + 3 * nctx);
    add_hidden_edge(base, previous, off + 4 * nctx);
  } else if (sub == 1 || sub == 2) {
    add_hidden_edge(base, previous, off + 2 * nctx);
    add_hidden_edge(base, previous, off + 3 * nctx);
    if (!minus && sub == 1) add_hidden_edge(base, previous, off + 4 * nctx);
  } else if (sub == 3) {
    add_hidden_edge(base, previous, off + 4 * nctx);
  } else if (minus) {
    add_coding_edges_minus(base, previous, state - 3, context);
  } else {
    add_coding_edges_plus(base, previous, state - 3, context);
  }
}

void crf_init(void) {
  static const uint8_t coding_motifs[2][3][6] = {
    {{1, 2, 2, 3, 0, 0}, {1, 2, 2, 2, 3, 3}, {1, 2, 3, 0, 0, 0}},
    {{3, 2, 1, 0, 0, 0}, {3, 2, 2, 1, 0, 0}, {3, 3, 2, 2, 2, 1}}
  };
  static const uint8_t intron_motifs[2][5] = {{4, 4, 5, 6, 6}, {6, 6, 5, 4, 4}};

  ntrans = 0;
  memset(hidden_state, 0, sizeof(hidden_state));
  memset(hidden_context, 0, sizeof(hidden_context));
  memset(hidden_motif, 0, sizeof(hidden_motif));
  memset(hidden_sub, 0, sizeof(hidden_sub));

  for (int state = 1; state < STATES; state++) {
    int off = hidden_offsets[state], cnt = hidden_counts[state], minus = state >= 7;
    int nsub = is_intron(state) ? 5 : 1, nctx = cnt / nsub;
    for (int sub = 0; sub < nsub; sub++)
      for (int c = 0; c < nctx; c++) {
        int h = off + sub * nctx + c;
        hidden_state[h] = (uint8_t)state;
        hidden_context[h] = (uint8_t)c;
        hidden_sub[h] = (uint8_t)sub;
        hidden_motif[h] = is_intron(state)
          ? intron_motifs[minus][sub]
          : coding_motifs[minus][(state - 1) % 3][c];
      }
  }

  memset(hidden_edge_count, 0, sizeof(hidden_edge_count));
  for (int base = 0; base < KNUC; base++)
    for (int previous = 0; previous < HIDDEN_STATES; previous++)
      add_grammar_edges(base, previous);
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
  double b_score[2][GFF_BOUND_KINDS] = {{0}};
  for (int rev = 0; rev < 2; rev++)
    for (int kind = 0; kind < GFF_BOUND_KINDS; kind++)
      b_score[rev][kind] = boundary_score(w, &site, rev, kind);
  for (int edge = 0; edge < ntrans; edge++) {
    int rev = edge_rev[edge], kind = edge_kind[edge];
    double score = w[edge_trans[edge]] + (kind >= 0 ? b_score[rev][kind] : 0.0);
    double p = score + emissions[trans_b[edge]];
    if (!isfinite(p))
      anno_fail("Non-finite CRF transition at position %d", position + 1);
    psi[edge] = p;
  }
  return site;
}

typedef struct {
  Site *sites;
  double *psi, *log_norm, *log_alpha, *aux;
  double *region_marg, *bound_marg, *pred_cwt, *res_cwt;
  double *grad_region, *grad_bound;
  int capacity;
} Work;

static void work_free(Work *k) {
  free(k->sites), free(k->psi), free(k->log_norm), free(k->log_alpha), free(k->aux);
  free(k->region_marg), free(k->bound_marg), free(k->pred_cwt), free(k->res_cwt);
  free(k->grad_region), free(k->grad_bound);
  memset(k, 0, sizeof(*k));
}

static void work_ensure(Work *k, int length) {
  if (length <= k->capacity) return;
  work_free(k);
  k->sites = anno_alloc((size_t)length, sizeof(Site));
  k->psi = anno_alloc((size_t)length, (size_t)ntrans * sizeof(double));
  k->log_norm = anno_alloc((size_t)length, sizeof(double));
  k->log_alpha = anno_alloc((size_t)length, HIDDEN_STATES * sizeof(double));
  k->aux = anno_alloc((size_t)length, HIDDEN_STATES * sizeof(double));
  k->region_marg = anno_alloc((size_t)length, GFF_REGIONS * sizeof(double));
  k->bound_marg = anno_alloc((size_t)length, (2 * GFF_BOUND_KINDS) * sizeof(double));
  k->pred_cwt = anno_alloc((size_t)length, CWT_CHANNELS * sizeof(double));
  k->res_cwt = anno_alloc((size_t)length, CWT_CHANNELS * sizeof(double));
  k->grad_region = anno_alloc((size_t)length, GFF_REGIONS * sizeof(double));
  k->grad_bound = anno_alloc((size_t)length, (2 * GFF_BOUND_KINDS) * sizeof(double));
  k->capacity = length;
}

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

static inline double edge_posterior(const double *prev_L, const double *L,
                                    const double *psi, double log_norm,
                                    const HiddenEdge *step) {
  return (prev_L[step->previous] > -INFINITY && L[step->current] > -INFINITY)
    ? exp(prev_L[step->previous] + psi[step->transition] - log_norm - L[step->current])
    : 0.0;
}

static void backward_step(const Work *k, int position, double *B, double *D,
                          const double *prev_A, const double *grad_reg_pos,
                          const double *grad_bound_pos, double c_bar,
                          double *edge_accum) {
  const Site *site = k->sites + position;
  const double *L = k->log_alpha + (size_t)position * HIDDEN_STATES;
  const double *prev_L = L - HIDDEN_STATES;
  const double *psi = k->psi + (size_t)position * ntrans;
  double log_norm = k->log_norm[position];
  double B_prev[HIDDEN_STATES] = {0}, D_prev[HIDDEN_STATES] = {0};
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
    if (!D) {
      if (edge_accum) edge_accum[tr] += edge_prob;
    } else {
      int kind = edge_kind[tr];
      double ce = kind >= 0 ? grad_bound_pos[edge_rev[tr] * GFF_BOUND_KINDS + kind] : 0.0;
      double cn = grad_reg_pos[gff_region(hidden_state[v])];
      double r_step = ce + cn;
      D_prev[u] += q * (B[v] * r_step + D[v]);
      if (edge_accum)
        edge_accum[tr] += q * (B[v] * (prev_A[u] + r_step - c_bar) + D[v]);
    }
  }
  double sum_b = 0;
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) sum_b += B_prev[hidden];
  if (!(sum_b > 0.0) || !isfinite(sum_b))
    anno_fail("Invalid CRF posterior at position %d", position + 1);
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) B[hidden] = B_prev[hidden] / sum_b;
  if (D)
    for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) D[hidden] = D_prev[hidden] / sum_b;
}

/*
 * Computes the differentiable (FASTA->GFF->CWT) vs (FASTA->CWT) residual loss
 * and accumulates its exact gradient with respect to w over one window.
 */
static double window_residual_gradient(const double *w, const Contig *c, Work *k,
                                       double *gradient, double res_weight,
                                       int start, int length,
                                       int residual_start, int residual_length) {
  if (start < 0 || length <= 0 || start > c->length - length ||
      residual_start < 0 || residual_length < 0 ||
      residual_start > length - residual_length)
    anno_fail("Invalid CRF window span");
  work_ensure(k, length);

  for (int position = 0; position < length; position++) {
    double emissions[STATES];
    double *psi = k->psi + (size_t)position * ntrans;
    Site site = site_scores(w, c, start + position, emissions, psi);
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

  const double *last_L = k->log_alpha + (size_t)(length - 1) * HIDDEN_STATES;
  double B[HIDDEN_STATES], D[HIDDEN_STATES] = {0};
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) B[hidden] = exp(last_L[hidden]);
  memset(k->region_marg, 0, (size_t)length * GFF_REGIONS * sizeof(*k->region_marg));
  memset(k->bound_marg, 0, (size_t)length * (2 * GFF_BOUND_KINDS) * sizeof(*k->bound_marg));
  for (int position = length - 1; position >= 0; position--) {
    double *rm = k->region_marg + (size_t)position * GFF_REGIONS;
    for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
      rm[gff_region(hidden_state[hidden])] += B[hidden];
    if (position) {
      double edge_m[STATES * STATES] = {0};
      backward_step(k, position, B, NULL, NULL, NULL, NULL, 0.0, edge_m);
      double *bm = k->bound_marg + (size_t)position * (2 * GFF_BOUND_KINDS);
      for (int edge = 0; edge < ntrans; edge++) {
        int kind = edge_kind[edge];
        if (kind >= 0)
          bm[edge_rev[edge] * GFF_BOUND_KINDS + kind] += edge_m[edge];
      }
    }
  }

  gff_cwt_predict_soft(k->region_marg, k->bound_marg, length, &crf_gff_cwt, k->pred_cwt);
  double r_loss = 0.0;
  for (int position = 0; position < length; position++) {
    int g = start + position;
    double *res = k->res_cwt + (size_t)position * CWT_CHANNELS;
    const double *pred = k->pred_cwt + (size_t)position * CWT_CHANNELS;
    if (!c->cwt || base_index(c->seq[g]) < 0 ||
        position < residual_start || position >= residual_start + residual_length) {
      memset(res, 0, CWT_CHANNELS * sizeof(*res));
      continue;
    }
    const float *obs = c->cwt + (size_t)g * CWT_CHANNELS;
    for (int ch = 0; ch < CWT_CHANNELS; ch++) {
      double diff = pred[ch] - (double)obs[ch];
      r_loss += 0.5 * res_weight * diff * diff;
      res[ch] = res_weight * diff;
    }
  }
  gff_cwt_residual_adjoint(k->res_cwt, length, &crf_gff_cwt, k->grad_region, k->grad_bound);

  for (int position = 0; position < length; position++) {
    const double *L = k->log_alpha + (size_t)position * HIDDEN_STATES;
    double *A = k->aux + (size_t)position * HIDDEN_STATES;
    const double *gr = k->grad_region + (size_t)position * GFF_REGIONS;
    if (!position) {
      for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
        A[hidden] = L[hidden] > -INFINITY ? gr[gff_region(hidden_state[hidden])] : 0.0;
    } else {
      const Site *site = k->sites + position;
      const double *prev_L = L - HIDDEN_STATES;
      const double *psi = k->psi + (size_t)position * ntrans;
      double log_norm = k->log_norm[position];
      const double *prev_A = A - HIDDEN_STATES;
      const double *gb = k->grad_bound + (size_t)position * (2 * GFF_BOUND_KINDS);
      memset(A, 0, HIDDEN_STATES * sizeof(*A));
      const HiddenEdge *edges = hidden_edges[site->nuc];
      int edge_count = hidden_edge_count[site->nuc];
      for (int edge = 0; edge < edge_count; edge++) {
        const HiddenEdge *step = edges + edge;
        double q = edge_posterior(prev_L, L, psi, log_norm, step);
        if (q == 0.0) continue;
        int tr = step->transition, kind = edge_kind[tr];
        double ce = kind >= 0 ? gb[edge_rev[tr] * GFF_BOUND_KINDS + kind] : 0.0;
        A[step->current] += q * (prev_A[step->previous] + ce);
      }
      for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
        if (L[hidden] > -INFINITY)
          A[hidden] += gr[gff_region(hidden_state[hidden])];
    }
  }

  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) B[hidden] = exp(last_L[hidden]);
  double c_bar = 0.0;
  const double *last_A = k->aux + (size_t)(length - 1) * HIDDEN_STATES;
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) c_bar += B[hidden] * last_A[hidden];

  for (int position = length - 1; position >= 0; position--) {
    const Site *site = k->sites + position;
    const double *A = k->aux + (size_t)position * HIDDEN_STATES;
    double state_grad[STATES] = {0};
    for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) {
      int st = hidden_state[hidden];
      state_grad[st] += B[hidden] * (A[hidden] - c_bar) + D[hidden];
    }
    for (int st = 0; st < STATES; st++)
      gradient[emission_index(site, st)] += state_grad[st];

    if (position) {
      double edge_grad[STATES * STATES] = {0};
      const double *prev_A = A - HIDDEN_STATES;
      const double *gr = k->grad_region + (size_t)position * GFF_REGIONS;
      const double *gb = k->grad_bound + (size_t)position * (2 * GFF_BOUND_KINDS);
      backward_step(k, position, B, D, prev_A, gr, gb, c_bar, edge_grad);
      for (int edge = 0; edge < ntrans; edge++) {
        double g = edge_grad[edge];
        if (g == 0.0) continue;
        gradient[edge_trans[edge]] += g;
        int kind = edge_kind[edge];
        if (kind >= 0)
          boundary_accumulate(gradient, site, edge_rev[edge], kind, g);
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
          anno_fail("Too many CRF windows");
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
  double *gradient;
  Work work;
} CrfWorker;

static void *crf_worker(void *arg) {
  CrfWorker *worker = arg;
  CrfTrain *train = worker->train;
  worker->res_loss = 0;
  worker->positions = 0;
  memset(worker->gradient, 0, WEIGHTS * sizeof(*worker->gradient));
  for (;;) {
    size_t slot = atomic_fetch_add_explicit(&train->next_window, 1, memory_order_relaxed);
    if (slot >= train->count) break;
    const CrfWindow *window = train->windows + slot;
    double res = window_residual_gradient(
      worker->w, window->contig, &worker->work, worker->gradient, CRF_RESIDUAL_WEIGHT,
      window->start, window->length, window->residual_start, window->residual_length);
    if (!isfinite(res)) anno_fail("Non-finite training residual");
    worker->res_loss += res;
    worker->positions += (size_t)window->residual_length;
  }
  return NULL;
}

typedef struct {
  CrfTrain train;
  CrfWorker *workers;
  pthread_t *threads;
  double res_loss;
  size_t positions;
} CrfObjective;

static double crf_objective(const double *w, double *gradient, void *context) {
  CrfObjective *objective = context;
  CrfTrain *train = &objective->train;
  double res_loss = 0;
  size_t positions = 0;

  memset(gradient, 0, WEIGHTS * sizeof(*gradient));
  atomic_init(&train->next_window, 0);
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
    res_loss += objective->workers[t].res_loss;
    positions += objective->workers[t].positions;
    for (int f = 0; f < WEIGHTS; f++) gradient[f] += objective->workers[t].gradient[f];
  }
  if (!isfinite(res_loss) || !positions) anno_fail("Non-finite training residual");

  double normalization = 1.0 / (double)positions;
  objective->res_loss = res_loss;
  objective->positions = positions;
  double penalty = 0;
  for (int f = 0; f < WEIGHTS; f++) {
    penalty += w[f] * w[f];
    gradient[f] = normalization * gradient[f] + 2 * CRF_L2_C2 * w[f];
    if (!isfinite(gradient[f])) anno_fail("Non-finite training gradient");
  }
  double loss = normalization * res_loss + CRF_L2_C2 * penalty;
  if (!isfinite(loss)) anno_fail("Non-finite training objective");
  return loss;
}

static void crf_progress(void *context, int iteration, double objective, double step,
                         int line_search_trials, size_t n, const double *gradient) {
  const CrfObjective *evaluation = context;
  double gradient_norm = 0;
  for (size_t f = 0; f < n; f++) gradient_norm = hypot(gradient_norm, gradient[f]);
  fprintf(stderr,
          "Iteration %d: %zu windows, cwt-residual/bp %.8g, objective %.8g, "
          "gradient-norm %.3g, step %.3g, line-search %d\n",
          iteration, evaluation->train.count,
          evaluation->res_loss / evaluation->positions,
          objective, gradient_norm, step, line_search_trials);
}

void crf_train(Contig *contigs, int count, double *w, int iterations, int threads) {
  if (count < 1 || iterations < 1)
    anno_fail("CRF training requires instances and positive iterations");
  if (threads < 1) anno_fail("CRF training requires at least one worker thread");
  if (!crf_gff_cwt_set)
    anno_fail("Self-training requires an embedded GFF->CWT model");
  for (int i = 0; i < count; i++)
    if (contigs[i].label)
      anno_fail("Self-training input must not contain GFF-derived labels");

  CrfWindow *windows = NULL;
  size_t window_count = make_windows(contigs, count, &windows);
  if (!window_count) anno_fail("No valid CRF training windows");
  CrfObjective objective = {
    .train = {.windows = windows, .count = window_count, .threads = threads}
  };
  fprintf(stderr, "Using %d worker threads for %zu training windows\n",
          objective.train.threads, objective.train.count);
  objective.workers = anno_alloc((size_t)threads, sizeof(*objective.workers));
  objective.threads = anno_alloc((size_t)threads, sizeof(*objective.threads));
  for (int t = 0; t < threads; t++)
    objective.workers[t] = (CrfWorker){
      .train = &objective.train,
      .gradient = anno_alloc(WEIGHTS, sizeof(*objective.workers->gradient))
    };

  int completed = 0;
  LbfgsStatus status = lbfgs_minimize(WEIGHTS, w, crf_objective, crf_progress, &objective,
                                      iterations, &completed);
  for (int t = 0; t < threads; t++) {
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

void crf_decode(const double *w, const Contig *c, uint8_t *path) {
  int length = c->length;
  uint8_t *back = anno_alloc((size_t)length, HIDDEN_STATES * sizeof(*back));
  double previous[HIDDEN_STATES] = {0}, current[HIDDEN_STATES];
  double emissions[STATES], transitions[STATES * STATES];

  for (int pos = 0; pos < length; pos++) {
    Site site = site_scores(w, c, pos, emissions, transitions);
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
