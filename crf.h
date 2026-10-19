#ifndef CRF_H
#define CRF_H
#include "anno.h"
#include "lbfgs.h"
#include <stdint.h>

#define CRF_L2_C2 1.0

/* Linear-chain CRF decoded jointly for both strands with strand-shared emissions.
   0 intergenic; 1-3 CDS(+) codon position; 4-6 intron(+) carrying the codon position of the
   next exon base; 7-9 CDS(-) codon position; 10-12 intron(-). */
static const int boundary_kmer_offsets[] = {-8, -6, -4, -2, 0, 2, 4, 6};
static const int boundary_cwt_offsets[] = {-27, -26, -25, -24, -2, -1, 0, 1, 2};
enum {
  CRF_WINDOW_SIZE = 65536,
  DERIV_COUNT = 3, DERIV_RADIUS = 2,
  STATES = 13, EMIT_STATES = 7, HIDDEN_STATES = 39,
  K2 = 17, K3 = 65, KNUC = 5,
  BOUNDARY_KMER_COUNT = sizeof(boundary_kmer_offsets) / sizeof(boundary_kmer_offsets[0]),
  BOUNDARY_CWT_COUNT = sizeof(boundary_cwt_offsets) / sizeof(boundary_cwt_offsets[0]),
  OFF_E3 = 0, OFF_SPLICE_START_STOP = OFF_E3 + EMIT_STATES * K3,
  OFF_BP = OFF_SPLICE_START_STOP + 4 * BOUNDARY_KMER_COUNT * K2,
  OFF_BD = OFF_BP + BOUNDARY_CWT_COUNT * WAVE_COUNT,
  WEIGHTS = OFF_BD + DERIV_COUNT * BOUNDARY_CWT_COUNT * WAVE_COUNT
};

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
  uint8_t k3_emit[2], nuc, boundary_k2[2][BOUNDARY_KMER_COUNT];
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

static const double derivative_kernels[DERIV_COUNT][2 * DERIV_RADIUS + 1] = {
  {0, -0.5, 0, 0.5, 0},
  {0, 1, -2, 1, 0},
  {-0.5, 1, 0, -1, 0.5}
};

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
}

static int emit_state(int state) {
  if (state < 7) return state;
  return state <= 9 ? state - 6 : 4 + (state - 9) % 3;
}

static int emission_index(const Site *site, int state) {
  return OFF_E3 + emit_state(state) * K3 + site->k3_emit[state >= 7];
}

typedef struct {
  Site *sites;
  double *exp_psi, *alpha;
  int capacity;
} Work;

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
  double b_score[2][4];
  for (int rev = 0; rev < 2; rev++) {
    for (int kind = 0; kind < 4; kind++) {
      double sp = 0;
      for (int offset = 0; offset < BOUNDARY_KMER_COUNT; offset++)
        sp += w[boundary_kmer_index(&site, rev, kind, offset)];
      b_score[rev][kind] = sp;
    }
    double cw = 0;
    for (int offset = 0; offset < BOUNDARY_CWT_COUNT; offset++) {
      float power[WAVE_COUNT];
      double derivative[DERIV_COUNT][WAVE_COUNT];
      boundary_cwt_at(c, genome_position, rev, offset, power, derivative);
      const double *wbp = w + OFF_BP + offset * WAVE_COUNT;
      for (int scale = 0; scale < WAVE_COUNT; scale++) {
        cw += wbp[scale] * power[scale];
        for (int slot = 0; slot < DERIV_COUNT; slot++)
          cw += w[derivative_index(slot, offset, scale)] * derivative[slot][scale];
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
  free(k->sites), free(k->exp_psi), free(k->alpha);
  k->sites = anno_alloc((size_t)length, sizeof(Site));
  k->exp_psi = anno_alloc((size_t)length, (size_t)ntrans * sizeof(double));
  k->alpha = anno_alloc((size_t)length, HIDDEN_STATES * sizeof(double));
  k->capacity = length;
}

static void work_free(Work *k) {
  free(k->sites), free(k->exp_psi), free(k->alpha);
}

/* Accumulates the exact observed-minus-expected gradient for one training window. */
static double contig_gradient(const double *w, const Contig *c,
                              Work *k, double *gradient) {
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
        if (initial_allowed(site.nuc, hidden) &&
            (c->length > 1 || terminal_allowed(hidden))) {
          double e = emissions[hidden_state[hidden]];
          if (e > max_e) max_e = e;
        }
      if (!isfinite(max_e)) anno_fail("Invalid CRF partition at position %d", position + 1);
      double sum = 0;
      for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) {
        double v = initial_allowed(site.nuc, hidden) &&
                   (c->length > 1 || terminal_allowed(hidden))
          ? exp(emissions[hidden_state[hidden]] - max_e) : 0.0;
        A[hidden] = v;
        sum += v;
      }
      if (!(sum > 0) || !isfinite(sum)) anno_fail("Invalid CRF partition at position %d", position + 1);
      double inv = 1.0 / sum;
      for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) A[hidden] *= inv;
      logz += max_e + log(sum);
    } else {
      int prev_gold = c->label[c->offset + position - 1];
      int gold_edge = pair_to_edge[prev_gold][gold_state];
      if (gold_edge < 0) anno_fail("Invalid transition in validated CRF training path");
      gold_score += psi[gold_edge] - emissions[gold_state];

      const double *prev_A = A - HIDDEN_STATES;
      const HiddenEdge *edges = hidden_edges[site.nuc];
      int edge_count = hidden_edge_count[site.nuc];
      int last = position == c->length - 1;
      double max_psi = -INFINITY;
      for (int edge = 0; edge < edge_count; edge++) {
        const HiddenEdge *step = edges + edge;
        if (prev_A[step->previous] > 0.0 && (!last || terminal_allowed(step->current))) {
          double p = psi[step->transition];
          if (p > max_psi) max_psi = p;
        }
      }
      if (!isfinite(max_psi))
        anno_fail("Invalid CRF partition at position %d", position + 1);

      double *ep = k->exp_psi + (size_t)position * ntrans;
      for (int edge = 0; edge < ntrans; edge++)
        ep[edge] = psi[edge] > max_psi ? 0.0 : exp(psi[edge] - max_psi);

      for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) A[hidden] = 0.0;
      for (int edge = 0; edge < edge_count; edge++) {
        const HiddenEdge *step = edges + edge;
        if (prev_A[step->previous] > 0.0 && (!last || terminal_allowed(step->current)))
          A[step->current] += prev_A[step->previous] * ep[step->transition];
      }

      double sum = 0;
      for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) sum += A[hidden];
      if (!(sum > 0) || !isfinite(sum))
        anno_fail("Invalid CRF partition at position %d", position + 1);
      double inv = 1.0 / sum;
      for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) A[hidden] *= inv;
      for (int edge = 0; edge < ntrans; edge++) ep[edge] *= inv;
      logz += max_psi + log(sum);
    }
  }

  double likelihood = gold_score - logz;
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

    for (int st = 0; st < STATES; st++)
      gradient[emission_index(site, st)] -= state_marginals[st];

    if (position) {
      const double *prev_A = A - HIDDEN_STATES;
      const double *ep = k->exp_psi + (size_t)position * ntrans;
      int previous_state = c->label[c->offset + position - 1];
      int gold_edge = pair_to_edge[previous_state][state];

      double ratio[HIDDEN_STATES];
      for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
        ratio[hidden] = A[hidden] > 0.0 && B[hidden] > 0.0 ? B[hidden] / A[hidden] : 0.0;

      double edge_marginals[STATES * STATES] = {0};
      double B_prev[HIDDEN_STATES] = {0};
      const HiddenEdge *edges = hidden_edges[site->nuc];
      int edge_count = hidden_edge_count[site->nuc];
      for (int edge = 0; edge < edge_count; edge++) {
        const HiddenEdge *step = edges + edge;
        double r = ratio[step->current];
        if (r > 0.0 && prev_A[step->previous] > 0.0) {
          double edge_prob = prev_A[step->previous] * ep[step->transition] * r;
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
  for (int f = 0; f < WEIGHTS; f++) objective->weights[f] = objective->scale[f] * v[f];
  double loss = crf_objective(objective->weights, gradient, context);
  for (int f = 0; f < WEIGHTS; f++) gradient[f] *= objective->scale[f];
  return loss;
}

/* scale_f = 1 / sqrt(max(1, sum of squared gold-path values of feature f)). */
static void feature_scales(const CrfWindow *windows, size_t count, double *scale) {
  memset(scale, 0, WEIGHTS * sizeof(*scale));
  for (size_t w = 0; w < count; w++) {
    Contig c = *windows[w].contig;
    c.offset = windows[w].start;
    c.length = windows[w].length;
    for (int i = 0; i < c.length; i++) {
      Site site = site_at(&c, i);
      int state = c.label[c.offset + i];
      scale[emission_index(&site, state)] += 1;
      if (!i) continue;
      int previous = c.label[c.offset + i - 1];
      int reverse = previous >= 7 || state >= 7, kind = boundary_kind(previous, state);
      if (kind >= 0) {
        for (int offset = 0; offset < BOUNDARY_KMER_COUNT; offset++)
          scale[boundary_kmer_index(&site, reverse, kind, offset)] += 1;
        if (kind == 1) {
          int genome_position = c.offset + i;
          for (int offset = 0; offset < BOUNDARY_CWT_COUNT; offset++) {
            float power[WAVE_COUNT];
            double derivative[DERIV_COUNT][WAVE_COUNT];
            boundary_cwt_at(&c, genome_position, reverse, offset, power, derivative);
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
    }
  }
  for (int f = 0; f < WEIGHTS; f++) scale[f] = 1 / sqrt(fmax(scale[f], 1));
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
  feature_scales(windows, window_count, objective.scale);
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
    current[hidden] = position == 0 && initial_allowed(site.nuc, hidden) ? emissions[hidden_state[hidden]] : -INFINITY;
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
  for (int s = 0; s < EMIT_STATES; s++) {
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
                idx, slot + 1, kind_names[1], disp, wave_sizes[sc], (float)w[idx]);
      }
    }
  }

  if (fclose(file)) anno_fail("Cannot finish writing weights text: %s", path);
}

#endif
