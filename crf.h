#ifndef CRF_H
#define CRF_H
#include "anno.h"
#include "lbfgs.h"
#include <stdint.h>

#define CRF_L2_C2 1.0
#define CRF_F1_BASE_WEIGHT 0.03
#define CRF_F1_TARGET_LOSS 0.25
#define CRF_MIN_FREQ 1

static const int boundary_kmer_offsets[] = {-4, -3, -2, -1, 0, 1, 2};
static const int boundary_cwt_offsets[] = {-26, -24, -6, -4, -2, -1, 0, 2, 4};
enum {
  CRF_WINDOW_SIZE = 65536,
  GC_RADIUS = 2048,
  DERIV_COUNT = 3, DERIV_RADIUS = 2,
  STATES = 13, EMIT_STATES = 7, HIDDEN_STATES = 137, TRANSITION_FEATURES = 21,
  K2 = 17, K3 = 65, KNUC = 5, EMIT_KMER_LAYERS = 1,
  BOUNDARY_PWM_FROM = -24, BOUNDARY_PWM_COUNT = 33,
  BOUNDARY_KMER_COUNT = sizeof(boundary_kmer_offsets) / sizeof(boundary_kmer_offsets[0]),
  BOUNDARY_CWT_COUNT = sizeof(boundary_cwt_offsets) / sizeof(boundary_cwt_offsets[0]),
  OFF_E3 = 0,
  OFF_SPLICE_START_STOP = OFF_E3 + EMIT_STATES * K3,
  OFF_BP = OFF_SPLICE_START_STOP + 4 * BOUNDARY_KMER_COUNT * K2,
  OFF_BD = OFF_BP + 4 * BOUNDARY_CWT_COUNT * WAVE_COUNT,
  OFF_PWM = OFF_BD + 4 * DERIV_COUNT * BOUNDARY_CWT_COUNT * WAVE_COUNT,
  OFF_BCROSS_D1 = OFF_PWM + 4 * BOUNDARY_PWM_COUNT * KNUC,
  OFF_BCROSS_D2 = OFF_BCROSS_D1 + 4 * BOUNDARY_KMER_COUNT * K2,
  OFF_BPWM_CWT = OFF_BCROSS_D2 + 4 * BOUNDARY_KMER_COUNT * K2,
  OFF_GC_DIFF = OFF_BPWM_CWT + 4 * BOUNDARY_CWT_COUNT * KNUC,
  OFF_GC_LVL = OFF_GC_DIFF + EMIT_STATES,
  OFF_E3_GC = OFF_GC_LVL + EMIT_STATES,
  OFF_TRANS_GC = OFF_E3_GC + EMIT_STATES * K3,
  OFF_TRANSITION = OFF_TRANS_GC + TRANSITION_FEATURES,
  OFF_F1_WEIGHT = OFF_TRANSITION + TRANSITION_FEATURES,
  WEIGHTS = OFF_F1_WEIGHT + 1
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
#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const int emit_kmer_sizes[EMIT_KMER_LAYERS] = {3};
static const int emit_kmer_offsets[EMIT_KMER_LAYERS] = {OFF_E3};
static const int emit_kmer_counts[EMIT_KMER_LAYERS] = {K3};

typedef struct {
  uint32_t kmer_emit[EMIT_KMER_LAYERS][2];
  uint8_t boundary_nuc[2][BOUNDARY_PWM_COUNT];
  uint8_t boundary_cwt_nuc[2][BOUNDARY_CWT_COUNT];
  uint8_t nuc, boundary_k2[2][BOUNDARY_KMER_COUNT];
  float gc_dev, gc_contrast;
} Site;

static const int hidden_offsets[STATES] = {0, 1, 5, 11, 14, 19, 39, 69, 72, 76, 82, 102, 132};
static const int hidden_counts[STATES] = {1, 4, 6, 3, 5, 20, 30, 3, 4, 6, 20, 30, 5};
static int trans_a[STATES * STATES], trans_b[STATES * STATES], ntrans;
static int edge_rev[STATES * STATES], edge_kind[STATES * STATES], edge_trans[STATES * STATES];
static int pair_to_edge[STATES][STATES];
static uint8_t hidden_state[HIDDEN_STATES], hidden_context[HIDDEN_STATES];
static uint8_t hidden_motif[HIDDEN_STATES], hidden_sub[HIDDEN_STATES];
typedef struct { uint8_t previous, current; uint16_t transition; } HiddenEdge;
static HiddenEdge hidden_edges[KNUC][HIDDEN_STATES * 4];
static int hidden_edge_count[KNUC];
static uint8_t rev_k2[K2];
static uint8_t k3_gc_count[K3];

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

static double local_gc_at(const Contig *c, int position) {
  if (c->genome_length < 128) return 0.40;
  int lo = position - GC_RADIUS < 0 ? 0 : position - GC_RADIUS;
  int hi = position + GC_RADIUS >= c->genome_length ? c->genome_length - 1 : position + GC_RADIUS;
  int gc = 0, atgc = 0;
  for (int g = lo; g <= hi; g++) {
    int b = base_index(c->seq[g]);
    if (b >= 0) {
      atgc++;
      if (b == 1 || b == 2) gc++;
    }
  }
  return atgc > 0 ? (double)gc / (double)atgc : 0.40;
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
    for (int offset = 0; offset < BOUNDARY_CWT_COUNT; offset++) {
      int displacement = boundary_cwt_offsets[offset];
      int g = rev ? position - displacement - 1 : position + displacement;
      int b = g >= 0 && g < c->genome_length ? base_index(c->seq[g]) : -1;
      site->boundary_cwt_nuc[rev][offset] = b < 0 ? 4 : rev ? 3 - b : b;
    }
  }
  if (c->genome_length >= 128) {
    double lgc = local_gc_at(c, position);
    uint32_t k3 = site->kmer_emit[0][0];
    site->gc_dev = (float)(lgc - 0.40);
    site->gc_contrast = k3 < 64 ? (float)((double)k3_gc_count[k3] / 3.0 - lgc) : 0.0f;
  } else {
    site->gc_dev = 0.0f;
    site->gc_contrast = 0.0f;
  }
}

static Site __attribute__((unused)) site_at(const Contig *c, int i) {
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

static int motif_transition_try(int prev_motif, int cur_motif, int reverse) {
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
  if (prev_motif < 0 || prev_motif >= 7 || cur_motif < 0 || cur_motif >= 7) return -1;
  int index = indices[prev_motif][cur_motif];
  return index < 0 ? -1 : OFF_TRANSITION + index;
}

static int motif_transition_index(int prev_motif, int cur_motif, int reverse) {
  int trans = motif_transition_try(prev_motif, cur_motif, reverse);
  if (trans < 0)
    anno_fail("Invalid CRF transition feature: %d -> %d", prev_motif, cur_motif);
  return trans;
}

static int motif_at(const uint8_t *y, int length, int i) {
  int st = y[i], minus = st >= 7;
  if (!st) return 0;
  if (is_intron(st)) {
    int dl = (i == 0 || y[i - 1] != st) ? 0 : (i == 1 || y[i - 2] != st) ? 1 : 2;
    int dr = (i + 1 == length || y[i + 1] != st) ? 0 :
             (i + 2 == length || y[i + 2] != st) ? 1 : 2;
    int d5 = minus ? dr : dl, d3 = minus ? dl : dr;
    return d5 == 0 ? 4 : d3 == 0 ? 6 : d5 == 1 ? 4 : d3 == 1 ? 6 : 5;
  }
  int p = i, q = i, left = minus ? 9 : 1, right = minus ? 7 : 3;
  while (p > 0 && y[p] != left) p--;
  if (p == 0 || y[p - 1] == 0) return minus ? 3 : 1;
  while (q + 1 < length && y[q] != right) q++;
  if (q + 1 == length || y[q + 1] == 0) return minus ? 1 : 3;
  return 2;
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

static int bp_index(int kind, int offset, int scale) {
  return OFF_BP + (kind * BOUNDARY_CWT_COUNT + offset) * WAVE_COUNT + scale;
}

static int derivative_index(int kind, int slot, int offset, int scale) {
  return OFF_BD + ((kind * DERIV_COUNT + slot) * BOUNDARY_CWT_COUNT + offset) * WAVE_COUNT + scale;
}

static int bcross_d1_index(const Site *site, int reverse, int kind, int offset) {
  return OFF_BCROSS_D1 + (kind * BOUNDARY_KMER_COUNT + offset) * K2 +
         site->boundary_k2[reverse][offset];
}

static int bcross_d2_index(const Site *site, int reverse, int kind, int offset) {
  return OFF_BCROSS_D2 + (kind * BOUNDARY_KMER_COUNT + offset) * K2 +
         site->boundary_k2[reverse][offset];
}

static int bpwm_cwt_index(const Site *site, int reverse, int kind, int offset) {
  return OFF_BPWM_CWT + (kind * BOUNDARY_CWT_COUNT + offset) * KNUC +
         site->boundary_cwt_nuc[reverse][offset];
}

static void cwt_at_pos(const Contig *c, int g, float *power,
                       double derivative[DERIV_COUNT][WAVE_COUNT]) {
  if (g < 0 || g >= c->genome_length || base_index(c->seq[g]) < 0) {
    memset(power, 0, WAVE_COUNT * sizeof(*power));
    memset(derivative, 0, DERIV_COUNT * WAVE_COUNT * sizeof(double));
    return;
  }
  memcpy(power, c->cwt_power + (size_t)g * WAVE_COUNT, WAVE_COUNT * sizeof(*power));
  if (c->genome_length >= 128 && g >= 16 && g + 16 < c->genome_length &&
      base_index(c->seq[g - 16]) >= 0 && base_index(c->seq[g + 16]) >= 0) {
    const float *p_left = c->cwt_power + (size_t)(g - 16) * WAVE_COUNT;
    const float *p_right = c->cwt_power + (size_t)(g + 16) * WAVE_COUNT;
    for (int scale = 0; scale < WAVE_COUNT; scale++)
      power[scale] -= 0.5f * (p_left[scale] + p_right[scale]);
  }
  for (int slot = 0; slot < DERIV_COUNT; slot++) {
    const double *kernel = derivative_kernels[slot];
    int valid = 1;
    for (int j = -DERIV_RADIUS; j <= DERIV_RADIUS; j++)
      if (kernel[j + DERIV_RADIUS] &&
          (j < -g || j >= c->genome_length - g || base_index(c->seq[g + j]) < 0)) {
        valid = 0;
        break;
      }
    for (int scale = 0; scale < WAVE_COUNT; scale++) {
      double value = 0;
      if (valid)
        for (int j = -DERIV_RADIUS; j <= DERIV_RADIUS; j++)
          if (kernel[j + DERIV_RADIUS])
            value += kernel[j + DERIV_RADIUS] * c->cwt_power[(size_t)(g + j) * WAVE_COUNT + scale];
      derivative[slot][scale] = value;
    }
  }
}

static void boundary_cwt_at(const Contig *c, int position, int reverse, int offset,
                            float *power, double derivative[DERIV_COUNT][WAVE_COUNT]) {
  int displacement = boundary_cwt_offsets[offset];
  cwt_at_pos(c, reverse ? position - displacement - 1 : position + displacement, power, derivative);
  if (reverse)
    for (int slot = 0; slot < DERIV_COUNT; slot += 2)
      for (int scale = 0; scale < WAVE_COUNT; scale++)
        derivative[slot][scale] = -derivative[slot][scale];
}

static void boundary_kmer_cwt_at(const Contig *c, int position, int reverse, int offset,
                                 double *d1_mean, double *d2_mean) {
  int displacement = boundary_kmer_offsets[offset];
  float power[WAVE_COUNT];
  double derivative[DERIV_COUNT][WAVE_COUNT];
  cwt_at_pos(c, reverse ? position - displacement - 1 : position + displacement, power, derivative);
  double s1 = 0, s2 = 0;
  for (int scale = 0; scale < WAVE_COUNT; scale++) {
    s1 += reverse ? -derivative[0][scale] : derivative[0][scale];
    s2 += derivative[1][scale];
  }
  *d1_mean = s1 * (1.0 / WAVE_COUNT);
  *d2_mean = s2 * (1.0 / WAVE_COUNT);
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
  if (pair_to_edge[a][b] < 0) pair_to_edge[a][b] = e;
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
  memset(pair_to_edge, -1, sizeof(pair_to_edge));
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

  for (int k = 0; k < 16; k++)
    rev_k2[k] = (uint8_t)(((3 - (k & 3)) << 2) | (3 - ((k >> 2) & 3)));
  rev_k2[16] = 16;
  for (int k = 0; k < 64; k++) {
    int b0 = (k >> 4) & 3, b1 = (k >> 2) & 3, b2 = k & 3;
    k3_gc_count[k] = (uint8_t)((b0 == 1 || b0 == 2) + (b1 == 1 || b1 == 2) + (b2 == 1 || b2 == 2));
  }
  k3_gc_count[64] = 0;
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
  int es = emit_state(state), rev = state >= 7;
  double score = 0;
  for (int layer = 0; layer < EMIT_KMER_LAYERS; layer++)
    score += w[emission_index(site, layer, state)];
  if (site->gc_dev != 0.0f || site->gc_contrast != 0.0f) {
    score += w[OFF_GC_DIFF + es] * site->gc_contrast;
    score += w[OFF_GC_LVL + es] * site->gc_dev;
    score += w[OFF_E3_GC + es * K3 + (int)site->kmer_emit[0][rev]] * site->gc_dev;
  }
  return score;
}

static void emission_accumulate(double *values, const Site *site, int state, double value) {
  int es = emit_state(state), rev = state >= 7;
  for (int layer = 0; layer < EMIT_KMER_LAYERS; layer++)
    values[emission_index(site, layer, state)] += value;
  if (site->gc_dev != 0.0f || site->gc_contrast != 0.0f) {
    values[OFF_GC_DIFF + es] += value * site->gc_contrast;
    values[OFF_GC_LVL + es] += value * site->gc_dev;
    values[OFF_E3_GC + es * K3 + (int)site->kmer_emit[0][rev]] += value * site->gc_dev;
  }
}

typedef struct {
  float power[WAVE_COUNT];
  double deriv[DERIV_COUNT][WAVE_COUNT];
  double d1_mean, d2_mean;
} CwtCached;

typedef struct {
  const Contig *contig;
  uint8_t nuc_ring[64];
  uint8_t k2_ring[64];
  CwtCached cwt_ring[64];
  int gc_count, atgc_count;
} SiteStream;

static inline void site_stream_fill_cwt(SiteStream *s, int g) {
  CwtCached *ce = &s->cwt_ring[g & 63];
  cwt_at_pos(s->contig, g, ce->power, ce->deriv);
  double s1 = 0, s2 = 0;
  for (int sc = 0; sc < WAVE_COUNT; sc++) {
    s1 += ce->deriv[0][sc];
    s2 += ce->deriv[1][sc];
  }
  ce->d1_mean = s1 * (1.0 / WAVE_COUNT);
  ce->d2_mean = s2 * (1.0 / WAVE_COUNT);
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
  s->gc_count = s->atgc_count = 0;
  if (c->genome_length >= 128) {
    int g0 = g_start - 1;
    int lo = g0 - GC_RADIUS < 0 ? 0 : g0 - GC_RADIUS;
    int hi = g0 + GC_RADIUS >= c->genome_length ? c->genome_length - 1 : g0 + GC_RADIUS;
    for (int g = lo; g <= hi; g++) {
      int b = base_index(c->seq[g]);
      if (b >= 0) {
        s->atgc_count++;
        if (b == 1 || b == 2) s->gc_count++;
      }
    }
  }
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

  if (c->genome_length >= 128) {
    int g_out = g0 - GC_RADIUS - 1;
    if (g_out >= 0 && g_out < c->genome_length) {
      int bo = base_index(c->seq[g_out]);
      if (bo >= 0) {
        s->atgc_count--;
        if (bo == 1 || bo == 2) s->gc_count--;
      }
    }
    int g_in = g0 + GC_RADIUS;
    if (g_in >= 0 && g_in < c->genome_length) {
      int bi = base_index(c->seq[g_in]);
      if (bi >= 0) {
        s->atgc_count++;
        if (bi == 1 || bi == 2) s->gc_count++;
      }
    }
  }

  Site site;
  site.nuc = s->nuc_ring[g0 & 63];
  uint8_t bm1 = s->nuc_ring[(g0 - 1) & 63], b0_c = site.nuc, bp1 = s->nuc_ring[(g0 + 1) & 63];
  if (bm1 < 4 && b0_c < 4 && bp1 < 4) {
    site.kmer_emit[0][0] = (uint32_t)((bm1 << 4) | (b0_c << 2) | bp1);
    site.kmer_emit[0][1] = (uint32_t)(((3 - bp1) << 4) | ((3 - b0_c) << 2) | (3 - bm1));
  } else {
    site.kmer_emit[0][0] = site.kmer_emit[0][1] = 64;
  }
  for (int offset = 0; offset < BOUNDARY_PWM_COUNT; offset++) {
    int displacement = BOUNDARY_PWM_FROM + offset;
    uint8_t bf = s->nuc_ring[(g0 + displacement) & 63];
    uint8_t br = s->nuc_ring[(g0 - displacement - 1) & 63];
    site.boundary_nuc[0][offset] = bf;
    site.boundary_nuc[1][offset] = br < 4 ? (uint8_t)(3 - br) : 4;
  }
  for (int offset = 0; offset < BOUNDARY_CWT_COUNT; offset++) {
    int displacement = boundary_cwt_offsets[offset];
    uint8_t bf = s->nuc_ring[(g0 + displacement) & 63];
    uint8_t br = s->nuc_ring[(g0 - displacement - 1) & 63];
    site.boundary_cwt_nuc[0][offset] = bf;
    site.boundary_cwt_nuc[1][offset] = br < 4 ? (uint8_t)(3 - br) : 4;
  }
  if (c->genome_length >= 128) {
    double lgc = s->atgc_count > 0 ? (double)s->gc_count / (double)s->atgc_count : 0.40;
    uint32_t k3 = site.kmer_emit[0][0];
    site.gc_dev = (float)(lgc - 0.40);
    site.gc_contrast = k3 < 64 ? (float)((double)k3_gc_count[k3] / 3.0 - lgc) : 0.0f;
  } else {
    site.gc_dev = 0.0f;
    site.gc_contrast = 0.0f;
  }
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

  double b_score[2][4] = {{0}};
  for (int rev = 0; rev < 2; rev++) {
    const uint8_t *k2 = site.boundary_k2[rev];
    for (int offset = 0; offset < BOUNDARY_KMER_COUNT; offset++) {
      int disp = boundary_kmer_offsets[offset];
      int g = rev ? g0 - disp - 1 : g0 + disp;
      const CwtCached *ce = &s->cwt_ring[g & 63];
      double d1_m = rev ? -ce->d1_mean : ce->d1_mean;
      double d2_m = ce->d2_mean;
      int k2_code = k2[offset];
      for (int kind = 0; kind < 4; kind++) {
        int idx = (kind * BOUNDARY_KMER_COUNT + offset) * K2 + k2_code;
        b_score[rev][kind] += w[OFF_SPLICE_START_STOP + idx] +
                              w[OFF_BCROSS_D1 + idx] * d1_m +
                              w[OFF_BCROSS_D2 + idx] * d2_m;
      }
    }
    for (int kind = 0; kind < 4; kind++) {
      double score = 0;
      for (int offset = 0; offset < BOUNDARY_PWM_COUNT; offset++)
        score += w[boundary_pwm_index(&site, rev, kind, offset)];
      b_score[rev][kind] += score;
    }
    for (int offset = 0; offset < BOUNDARY_CWT_COUNT; offset++) {
      int displacement = boundary_cwt_offsets[offset];
      int g = rev ? g0 - displacement - 1 : g0 + displacement;
      const CwtCached *ce = &s->cwt_ring[g & 63];
      double d2_m = ce->d2_mean;
      uint8_t cnuc = site.boundary_cwt_nuc[rev][offset];
      for (int kind = 0; kind < 4; kind++) {
        const double *wbp = w + bp_index(kind, offset, 0);
        const double *wd0 = w + derivative_index(kind, 0, offset, 0);
        const double *wd1 = w + derivative_index(kind, 1, offset, 0);
        const double *wd2 = w + derivative_index(kind, 2, offset, 0);
        double cw = w[OFF_BPWM_CWT + (kind * BOUNDARY_CWT_COUNT + offset) * KNUC + cnuc] * d2_m;
        for (int scale = 0; scale < WAVE_COUNT; scale++) {
          cw += wbp[scale] * ce->power[scale];
          double d0 = ce->deriv[0][scale];
          double d1 = ce->deriv[1][scale];
          double d2 = ce->deriv[2][scale];
          if (rev) { d0 = -d0; d2 = -d2; }
          cw += wd0[scale] * d0 + wd1[scale] * d1 + wd2[scale] * d2;
        }
        b_score[rev][kind] += cw;
      }
    }
  }

  for (int edge = 0; edge < ntrans; edge++) {
    int reverse = edge_rev[edge], kind = edge_kind[edge];
    int t_id = edge_trans[edge] - OFF_TRANSITION;
    double score = w[edge_trans[edge]];
    if (site.gc_dev != 0.0f) score += w[OFF_TRANS_GC + t_id] * site.gc_dev;
    if (kind >= 0) score += b_score[reverse][kind];
    double p = score + emissions[trans_b[edge]];
    if (!isfinite(p))
      anno_fail("Non-finite CRF transition at position %d", position + 1);
    psi[edge] = p;
  }
  return site;
}

static int path_edge(const uint8_t *y, int length, int position) {
  int a = y[position - 1], b = y[position];
  int trans = motif_transition_try(motif_at(y, length, position - 1),
                                   motif_at(y, length, position), a >= 7 || b >= 7);
  if (trans >= 0)
    for (int e = 0; e < ntrans; e++)
      if (trans_a[e] == a && trans_b[e] == b && edge_trans[e] == trans)
        return e;
  return pair_to_edge[a][b];
}

static int gold_initial_hidden(int base, int state, int length) {
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
    if (hidden_state[hidden] == state && initial_allowed(base, hidden) &&
        (length > 1 || terminal_allowed(hidden)))
      return hidden;
  return -1;
}

static int gold_step_edge(int previous_hidden, int base, const uint8_t *y,
                          int length, int position, int *next_hidden) {
  int state = y[position], last = position == length - 1;
  int motif = motif_at(y, length, position);
  const HiddenEdge *edges = hidden_edges[base];
  int edge_count = hidden_edge_count[base];
  for (int e = 0; e < edge_count; e++) {
    const HiddenEdge *step = edges + e;
    if (step->previous != previous_hidden || trans_b[step->transition] != state)
      continue;
    int h = step->current;
    if (last && !terminal_allowed(h)) continue;
    if ((state == 1 || state == 9) && hidden_motif[h] != motif) continue;
    if (is_intron(state)) {
      if (hidden_motif[h] != motif) continue;
      if (hidden_sub[h] >= 3 && (hidden_sub[h] == 4) != (last || y[position + 1] != state))
        continue;
    }
    if (next_hidden) *next_hidden = h;
    return step->transition;
  }
  return -1;
}

static int gold_path_valid(const Contig *c) {
  const uint8_t *y = c->label + c->offset;
  int b0 = base_index(c->seq[c->offset]);
  int hidden = gold_initial_hidden(b0 < 0 ? 4 : b0, y[0], c->length);
  if (hidden < 0) return 0;
  for (int position = 1; position < c->length; position++) {
    int base = base_index(c->seq[c->offset + position]);
    if (gold_step_edge(hidden, base < 0 ? 4 : base, y, c->length, position, &hidden) < 0)
      return 0;
  }
  return terminal_allowed(hidden);
}

typedef struct {
  Site *sites;
  CwtCached *cwt;
  double *psi, *log_norm, *log_alpha, *aux;
  int capacity;
} Work;

static void work_ensure(Work *k, int length) {
  if (length <= k->capacity) return;
  free(k->sites), free(k->cwt), free(k->psi), free(k->log_norm), free(k->log_alpha), free(k->aux);
  k->sites = anno_alloc((size_t)length, sizeof(Site));
  k->cwt = anno_alloc((size_t)length + 64, sizeof(CwtCached));
  k->psi = anno_alloc((size_t)length, (size_t)ntrans * sizeof(double));
  k->log_norm = anno_alloc((size_t)length, sizeof(double));
  k->log_alpha = anno_alloc((size_t)length, HIDDEN_STATES * sizeof(double));
  k->aux = anno_alloc(2 * (size_t)length, HIDDEN_STATES * sizeof(double));
  k->capacity = length;
}

static void work_free(Work *k) {
  free(k->sites), free(k->cwt), free(k->psi), free(k->log_norm), free(k->log_alpha), free(k->aux);
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

static inline void feature_observe(double *grad, double *scale, double *freq,
                                   int f, double coeff, double value) {
  if (grad) grad[f] += coeff * value;
  if (scale) scale[f] += value * value;
  if (freq && value != 0) freq[f] += 1;
}

static void boundary_observe(const Contig *c, const Site *site, int genome_pos,
                             int rev, int kind, double coeff,
                             double *grad, double *scale, double *freq) {
  for (int offset = 0; offset < BOUNDARY_KMER_COUNT; offset++) {
    feature_observe(grad, scale, freq, boundary_kmer_index(site, rev, kind, offset), coeff, 1);
    double d1_m = 0, d2_m = 0;
    boundary_kmer_cwt_at(c, genome_pos, rev, offset, &d1_m, &d2_m);
    feature_observe(grad, scale, freq, bcross_d1_index(site, rev, kind, offset), coeff, d1_m);
    feature_observe(grad, scale, freq, bcross_d2_index(site, rev, kind, offset), coeff, d2_m);
  }
  for (int offset = 0; offset < BOUNDARY_PWM_COUNT; offset++)
    feature_observe(grad, scale, freq, boundary_pwm_index(site, rev, kind, offset), coeff, 1);
  for (int offset = 0; offset < BOUNDARY_CWT_COUNT; offset++) {
    float power[WAVE_COUNT];
    double derivative[DERIV_COUNT][WAVE_COUNT];
    boundary_cwt_at(c, genome_pos, rev, offset, power, derivative);
    double d2_m = 0;
    for (int sc = 0; sc < WAVE_COUNT; sc++) {
      feature_observe(grad, scale, freq, bp_index(kind, offset, sc), coeff, power[sc]);
      for (int slot = 0; slot < DERIV_COUNT; slot++)
        feature_observe(grad, scale, freq, derivative_index(kind, slot, offset, sc),
                        coeff, derivative[slot][sc]);
      d2_m += derivative[1][sc];
    }
    d2_m *= (1.0 / WAVE_COUNT);
    feature_observe(grad, scale, freq, bpwm_cwt_index(site, rev, kind, offset), coeff, d2_m);
  }
}

static void boundary_observe_cached(const CwtCached *cwt_center, const Site *site,
                                    int rev, int kind, double coeff, double *grad) {
  for (int offset = 0; offset < BOUNDARY_KMER_COUNT; offset++) {
    int disp = boundary_kmer_offsets[offset];
    const CwtCached *ce = cwt_center + (rev ? -disp - 1 : disp);
    double d1_m = rev ? -ce->d1_mean : ce->d1_mean;
    double d2_m = ce->d2_mean;
    grad[boundary_kmer_index(site, rev, kind, offset)] += coeff;
    grad[bcross_d1_index(site, rev, kind, offset)] += coeff * d1_m;
    grad[bcross_d2_index(site, rev, kind, offset)] += coeff * d2_m;
  }
  for (int offset = 0; offset < BOUNDARY_PWM_COUNT; offset++)
    grad[boundary_pwm_index(site, rev, kind, offset)] += coeff;
  for (int offset = 0; offset < BOUNDARY_CWT_COUNT; offset++) {
    int disp = boundary_cwt_offsets[offset];
    const CwtCached *ce = cwt_center + (rev ? -disp - 1 : disp);
    grad[bpwm_cwt_index(site, rev, kind, offset)] += coeff * ce->d2_mean;
    double *gbp = grad + bp_index(kind, offset, 0);
    double *gd0 = grad + derivative_index(kind, 0, offset, 0);
    double *gd1 = grad + derivative_index(kind, 1, offset, 0);
    double *gd2 = grad + derivative_index(kind, 2, offset, 0);
    for (int sc = 0; sc < WAVE_COUNT; sc++) {
      gbp[sc] += coeff * ce->power[sc];
      double d0 = ce->deriv[0][sc], d1 = ce->deriv[1][sc], d2 = ce->deriv[2][sc];
      if (rev) { d0 = -d0; d2 = -d2; }
      gd0[sc] += coeff * d0;
      gd1[sc] += coeff * d1;
      gd2[sc] += coeff * d2;
    }
  }
}

static inline double f1_pred_weight(int st, int gold_state, int length) {
  if (is_cds(st)) return 1.0;
  if (length >= 64 && gold_state == 0 && is_intron(st)) return 0.5;
  return 0.0;
}

static void backward_step(const Work *k, int length, int position, double *B, double *D,
                          const double *prev_Ap, const double *prev_At,
                          double c_pred, double c_tp, int gold_state, double r_bar,
                          double *edge_marginals) {
  const Site *site = k->sites + position;
  const double *L = k->log_alpha + (size_t)position * HIDDEN_STATES;
  const double *prev_L = L - HIDDEN_STATES;
  const double *psi = k->psi + (size_t)position * ntrans;
  double log_norm = k->log_norm[position];
  int last = position == length - 1;
  double B_prev[HIDDEN_STATES] = {0}, D_prev[HIDDEN_STATES] = {0};
  const HiddenEdge *edges = hidden_edges[site->nuc];
  int edge_count = hidden_edge_count[site->nuc];

  double psi_max = -INFINITY;
  for (int t = 0; t < ntrans; t++)
    if (psi[t] > psi_max) psi_max = psi[t];
  int fast = isfinite(psi_max);
  for (int t = 0; fast && t < ntrans; t++)
    if (psi[t] - psi_max < -600.0) fast = 0;
  if (fast) {
    double exp_prev[HIDDEN_STATES] = {0}, inv_L[HIDDEN_STATES] = {0}, exp_psi[STATES * STATES];
    for (int t = 0; t < ntrans; t++) exp_psi[t] = exp(psi[t] - psi_max);
    for (int h = 0; h < HIDDEN_STATES; h++) {
      if (prev_L[h] >= -600.0) exp_prev[h] = exp(prev_L[h]);
      if (B[h] > 0.0 && L[h] > -INFINITY && (!last || terminal_allowed(h))) {
        double shift = psi_max - log_norm - L[h];
        inv_L[h] = fabs(shift) <= 600.0 ? exp(shift) : -1.0;
      }
    }
    for (int edge = 0; edge < edge_count; edge++) {
      const HiddenEdge *step = edges + edge;
      int u = step->previous, v = step->current;
      double inv_v = inv_L[v];
      if (inv_v == 0.0) continue;
      double ep_u = exp_prev[u];
      double q;
      if (inv_v > 0.0 && ep_u > 0.0) {
        q = ep_u * exp_psi[step->transition] * inv_v;
      } else {
        if (!(prev_L[u] > -INFINITY)) continue;
        q = exp(prev_L[u] + psi[step->transition] - log_norm - L[v]);
      }
      double edge_prob = B[v] * q;
      B_prev[u] += edge_prob;
      if (edge_marginals) edge_marginals[step->transition] += edge_prob;
      if (D) {
        int st = hidden_state[v];
        double pw = f1_pred_weight(st, gold_state, length);
        double rv = pw * c_pred + (is_cds(st) && st == gold_state ? c_tp : 0.0);
        double au = c_pred * prev_Ap[u] + c_tp * prev_At[u];
        D_prev[u] += q * (B[v] * rv + D[v]);
        edge_marginals[step->transition] += q * (B[v] * (au + rv - r_bar) + D[v]);
      }
    }
  } else {
    for (int edge = 0; edge < edge_count; edge++) {
      const HiddenEdge *step = edges + edge;
      int u = step->previous, v = step->current;
      if (!(B[v] > 0.0) || !forward_edge_allowed(prev_L, step, last)) continue;
      double q = exp(prev_L[u] + psi[step->transition] - log_norm - L[v]);
      double edge_prob = B[v] * q;
      B_prev[u] += edge_prob;
      if (edge_marginals) edge_marginals[step->transition] += edge_prob;
      if (D) {
        int st = hidden_state[v];
        double pw = f1_pred_weight(st, gold_state, length);
        double rv = pw * c_pred + (is_cds(st) && st == gold_state ? c_tp : 0.0);
        double au = c_pred * prev_Ap[u] + c_tp * prev_At[u];
        D_prev[u] += q * (B[v] * rv + D[v]);
        edge_marginals[step->transition] += q * (B[v] * (au + rv - r_bar) + D[v]);
      }
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
 * Accumulates the exact observed-minus-expected gradient (and optional differentiable Soft F1
 * loss gradient) for one training window in a single forward-backward pass.
 */
static double contig_gradient_f1(const double *w, const Contig *c, Work *k,
                                 double *gradient, double f1_weight, double *f1_loss) {
  work_ensure(k, c->length);
  SiteStream stream;
  site_stream_init(&stream, c, 0);
  for (int d = -27; d <= 25; d++)
    k->cwt[32 + d] = stream.cwt_ring[(c->offset + d) & 63];
  const uint8_t *gold_y = c->label + c->offset;
  int use_f1 = f1_loss != NULL && f1_weight > 0.0;
  double logz = 0, gold_score = 0, gold_pos = 0;
  for (int position = 0; position < c->length; position++) {
    double emissions[STATES];
    double *psi = k->psi + (size_t)position * ntrans;
    Site site = site_scores_stream(w, &stream, position, emissions, psi);
    k->sites[position] = site;
    k->cwt[32 + position + 26] = stream.cwt_ring[(c->offset + position + 26) & 63];
    int gold_state = gold_y[position];
    gold_score += emissions[gold_state];
    if (use_f1 && is_cds(gold_state)) gold_pos += 1.0;

    double *L = k->log_alpha + (size_t)position * HIDDEN_STATES;
    double *Ap = k->aux + 2 * (size_t)position * HIDDEN_STATES;
    double *At = Ap + HIDDEN_STATES;
    if (!position) {
      for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) {
        int st = hidden_state[hidden];
        L[hidden] = initial_allowed(site.nuc, hidden) &&
                    (c->length > 1 || terminal_allowed(hidden))
          ? emissions[st] : -INFINITY;
        if (use_f1) {
          Ap[hidden] = L[hidden] > -INFINITY ? f1_pred_weight(st, gold_state, c->length) : 0.0;
          At[hidden] = L[hidden] > -INFINITY && is_cds(st) && st == gold_state;
        }
      }
    } else {
      int gold_edge = path_edge(gold_y, c->length, position);
      if (gold_edge < 0) anno_fail("Invalid transition in validated CRF training path");
      gold_score += psi[gold_edge] - emissions[gold_state];

      const double *prev_L = L - HIDDEN_STATES;
      const double *prev_Ap = Ap - 2 * HIDDEN_STATES, *prev_At = prev_Ap + HIDDEN_STATES;
      const HiddenEdge *edges = hidden_edges[site.nuc];
      int edge_count = hidden_edge_count[site.nuc];
      int last = position == c->length - 1;
      double incoming_sum[HIDDEN_STATES] = {0};
      for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) {
        L[hidden] = -INFINITY;
        if (use_f1) Ap[hidden] = At[hidden] = 0;
      }
      double psi_max = -INFINITY;
      for (int t = 0; t < ntrans; t++)
        if (psi[t] > psi_max) psi_max = psi[t];
      int fast = isfinite(psi_max);
      for (int t = 0; fast && t < ntrans; t++)
        if (psi[t] - psi_max < -600.0) fast = 0;

      if (fast) {
        double exp_prev[HIDDEN_STATES] = {0}, exp_psi[STATES * STATES];
        uint8_t reachable[HIDDEN_STATES] = {0}, slow[HIDDEN_STATES] = {0};
        int any_slow = 0;
        for (int t = 0; t < ntrans; t++) exp_psi[t] = exp(psi[t] - psi_max);
        for (int h = 0; h < HIDDEN_STATES; h++)
          if (prev_L[h] >= -600.0) exp_prev[h] = exp(prev_L[h]);
        for (int edge = 0; edge < edge_count; edge++) {
          const HiddenEdge *step = edges + edge;
          if (!forward_edge_allowed(prev_L, step, last)) continue;
          reachable[step->current] = 1;
          double w_edge = exp_prev[step->previous] * exp_psi[step->transition];
          incoming_sum[step->current] += w_edge;
          if (use_f1) {
            Ap[step->current] += w_edge * prev_Ap[step->previous];
            At[step->current] += w_edge * prev_At[step->previous];
          }
        }
        for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
          if (reachable[hidden]) {
            if (incoming_sum[hidden] >= 1e-250) {
              if (use_f1) {
                int st = hidden_state[hidden];
                Ap[hidden] = Ap[hidden] / incoming_sum[hidden] +
                             f1_pred_weight(st, gold_state, c->length);
                At[hidden] = At[hidden] / incoming_sum[hidden] + (is_cds(st) && st == gold_state);
              }
              L[hidden] = psi_max + log(incoming_sum[hidden]);
            } else {
              slow[hidden] = 1;
              any_slow = 1;
              incoming_sum[hidden] = 0.0;
              if (use_f1) Ap[hidden] = At[hidden] = 0.0;
            }
          }
        if (any_slow) {
          for (int edge = 0; edge < edge_count; edge++) {
            const HiddenEdge *step = edges + edge;
            if (!slow[step->current] || !forward_edge_allowed(prev_L, step, last)) continue;
            double score = prev_L[step->previous] + psi[step->transition];
            if (score > L[step->current]) L[step->current] = score;
          }
          for (int edge = 0; edge < edge_count; edge++) {
            const HiddenEdge *step = edges + edge;
            if (!slow[step->current] || !forward_edge_allowed(prev_L, step, last)) continue;
            double w_edge = exp(prev_L[step->previous] + psi[step->transition] - L[step->current]);
            incoming_sum[step->current] += w_edge;
            if (use_f1) {
              Ap[step->current] += w_edge * prev_Ap[step->previous];
              At[step->current] += w_edge * prev_At[step->previous];
            }
          }
          for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
            if (slow[hidden] && L[hidden] > -INFINITY) {
              if (use_f1) {
                int st = hidden_state[hidden];
                Ap[hidden] = Ap[hidden] / incoming_sum[hidden] +
                             f1_pred_weight(st, gold_state, c->length);
                At[hidden] = At[hidden] / incoming_sum[hidden] + (is_cds(st) && st == gold_state);
              }
              L[hidden] += log(incoming_sum[hidden]);
            }
        }
      } else {
        for (int edge = 0; edge < edge_count; edge++) {
          const HiddenEdge *step = edges + edge;
          if (!forward_edge_allowed(prev_L, step, last)) continue;
          double score = prev_L[step->previous] + psi[step->transition];
          if (score > L[step->current]) L[step->current] = score;
        }
        for (int edge = 0; edge < edge_count; edge++) {
          const HiddenEdge *step = edges + edge;
          if (!forward_edge_allowed(prev_L, step, last)) continue;
          double w_edge = exp(prev_L[step->previous] + psi[step->transition] - L[step->current]);
          incoming_sum[step->current] += w_edge;
          if (use_f1) {
            Ap[step->current] += w_edge * prev_Ap[step->previous];
            At[step->current] += w_edge * prev_At[step->previous];
          }
        }
        for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
          if (L[hidden] > -INFINITY) {
            if (use_f1) {
              int st = hidden_state[hidden];
              Ap[hidden] = Ap[hidden] / incoming_sum[hidden] +
                           f1_pred_weight(st, gold_state, c->length);
              At[hidden] = At[hidden] / incoming_sum[hidden] + (is_cds(st) && st == gold_state);
            }
            L[hidden] += log(incoming_sum[hidden]);
          }
      }
    }
    k->log_norm[position] = normalize_log_alpha(L, position);
    logz += k->log_norm[position];
  }

  double likelihood = gold_score - logz;
  const double *last_L = k->log_alpha + (size_t)(c->length - 1) * HIDDEN_STATES;
  double B[HIDDEN_STATES], D[HIDDEN_STATES] = {0};
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) B[hidden] = exp(last_L[hidden]);

  double c_pred = 0, c_tp = 0, r_bar = 0;
  if (use_f1) {
    const double *last_Ap = k->aux + 2 * (size_t)(c->length - 1) * HIDDEN_STATES;
    const double *last_At = last_Ap + HIDDEN_STATES;
    double pred_pos = 0, tp = 0;
    for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) {
      pred_pos += B[hidden] * last_Ap[hidden];
      tp += B[hidden] * last_At[hidden];
    }
    double s = w[OFF_F1_WEIGHT], exp_s = exp(s), exp_neg_s = exp(-s);
    double denom = pred_pos + gold_pos;
    double base_scale = gold_pos > 0.0 ? f1_weight * gold_pos : f1_weight;
    double raw_loss = gold_pos > 0.0 ? 1.0 - 2.0 * tp / denom : pred_pos;
    double scale = base_scale * exp_neg_s;
    *f1_loss = scale * raw_loss + base_scale * CRF_F1_TARGET_LOSS * exp_s;
    gradient[OFF_F1_WEIGHT] -= base_scale * (CRF_F1_TARGET_LOSS * exp_s - exp_neg_s * raw_loss);
    if (gold_pos > 0.0) {
      c_pred = 2.0 * scale * tp / (denom * denom);
      c_tp = -2.0 * scale / denom;
    } else {
      c_pred = scale;
    }
    r_bar = c_pred * pred_pos + c_tp * tp;
  } else if (f1_loss) {
    *f1_loss = 0;
  }

  for (int position = c->length - 1; position >= 0; position--) {
    const Site *site = k->sites + position;
    const double *Ap = k->aux + 2 * (size_t)position * HIDDEN_STATES;
    const double *At = Ap + HIDDEN_STATES;
    int state = gold_y[position];

    double state_marginals[STATES] = {0};
    for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) {
      int st = hidden_state[hidden];
      state_marginals[st] += B[hidden];
      if (use_f1) {
        double ah = c_pred * Ap[hidden] + c_tp * At[hidden];
        state_marginals[st] += B[hidden] * (ah - r_bar) + D[hidden];
      }
    }
    state_marginals[state] -= 1.0;

    for (int st = 0; st < STATES; st++)
      emission_accumulate(gradient, site, st, -state_marginals[st]);

    if (position) {
      int gold_edge = path_edge(gold_y, c->length, position);
      double edge_marginals[STATES * STATES] = {0};
      const double *prev_Ap = Ap - 2 * HIDDEN_STATES;
      backward_step(k, c->length, position, B, use_f1 ? D : NULL,
                    use_f1 ? prev_Ap : NULL, use_f1 ? prev_Ap + HIDDEN_STATES : NULL,
                    c_pred, c_tp, state, r_bar, edge_marginals);

      double b_coeff[4][2] = {{0}};
      int has_boundary = 0;
      edge_marginals[gold_edge] -= 1.0;

      for (int edge = 0; edge < ntrans; edge++) {
        double m = edge_marginals[edge];
        if (m != 0) {
          int t_id = edge_trans[edge] - OFF_TRANSITION;
          gradient[edge_trans[edge]] -= m;
          if (site->gc_dev != 0.0f) gradient[OFF_TRANS_GC + t_id] -= m * site->gc_dev;
          int reverse = edge_rev[edge], kind = edge_kind[edge];
          if (kind >= 0) {
            b_coeff[kind][reverse] -= m;
            has_boundary = 1;
          }
        }
      }
      if (has_boundary)
        for (int rev = 0; rev < 2; rev++)
          for (int kind = 0; kind < 4; kind++)
            if (b_coeff[kind][rev] != 0)
              boundary_observe_cached(k->cwt + 32 + position, site, rev, kind,
                                      b_coeff[kind][rev], gradient);
    }
  }
  return likelihood;
}

static inline double contig_gradient(const double *w, const Contig *c,
                                     Work *k, double *gradient) {
  return contig_gradient_f1(w, c, k, gradient, 0, NULL);
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
  atomic_size_t next_window;
  int threads;
} CrfTrain;

typedef struct {
  CrfTrain *train;
  int thread_id;
  const double *w;
  double likelihood, f1_loss;
  size_t positions;
  double *gradient;
  Work work;
} CrfWorker;

static void *crf_worker(void *arg) {
  CrfWorker *worker = arg;
  CrfTrain *train = worker->train;
  worker->likelihood = 0;
  worker->f1_loss = 0;
  worker->positions = 0;
  memset(worker->gradient, 0, WEIGHTS * sizeof(*worker->gradient));
  for (;;) {
    size_t slot = atomic_fetch_add_explicit(&train->next_window, 1, memory_order_relaxed);
    if (slot >= train->count) break;
    const CrfWindow *window = train->windows + slot;
    Contig instance = *window->contig;
    instance.offset = window->start;
    instance.length = window->length;
    double f1 = 0;
    double score = contig_gradient_f1(worker->w, &instance, &worker->work,
                                      worker->gradient, CRF_F1_BASE_WEIGHT, &f1);
    if (!isfinite(score) || !isfinite(f1)) anno_fail("Non-finite training likelihood");
    worker->likelihood += score;
    worker->f1_loss += f1;
    worker->positions += (size_t)instance.length;
  }
  return NULL;
}

typedef struct {
  CrfTrain train;
  CrfWorker *workers;
  pthread_t *threads;
  double likelihood, f1_loss;
  size_t positions;
  double *scale, *weights;
} CrfObjective;

static double crf_objective(const double *w, double *gradient, void *context) {
  CrfObjective *objective = context;
  CrfTrain *train = &objective->train;
  double likelihood = 0, f1_loss = 0;
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
    likelihood += objective->workers[t].likelihood;
    f1_loss += objective->workers[t].f1_loss;
    positions += objective->workers[t].positions;
    for (int f = 0; f < WEIGHTS; f++) gradient[f] += objective->workers[t].gradient[f];
  }
  if (!isfinite(likelihood) || !isfinite(f1_loss)) anno_fail("Non-finite training likelihood");

  objective->likelihood = likelihood;
  objective->f1_loss = f1_loss;
  objective->positions = positions;
  double penalty = 0;
  for (int f = 0; f < WEIGHTS; f++) {
    penalty += w[f] * w[f];
    gradient[f] = -gradient[f] + 2 * CRF_L2_C2 * w[f];
    if (!isfinite(gradient[f])) anno_fail("Non-finite training gradient");
  }
  double loss = -likelihood + f1_loss + CRF_L2_C2 * penalty;
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

static void feature_counts(const CrfWindow *windows, size_t count, double *scale, double *freq) {
  memset(scale, 0, WEIGHTS * sizeof(*scale));
  if (freq) memset(freq, 0, WEIGHTS * sizeof(*freq));
  double dummy_w[WEIGHTS] = {0}, dummy_e[STATES], dummy_psi[STATES * STATES];
  int has_genome_window = 0;
  for (size_t w = 0; w < count; w++) {
    Contig c = *windows[w].contig;
    c.offset = windows[w].start;
    c.length = windows[w].length;
    if (c.length >= 64) has_genome_window = 1;
    const uint8_t *y = c.label + c.offset;
    SiteStream stream;
    site_stream_init(&stream, &c, 0);
    for (int i = 0; i < c.length; i++) {
      Site site = site_scores_stream(dummy_w, &stream, i, dummy_e, dummy_psi);
      int state = y[i], es = emit_state(state), rev = state >= 7;
      for (int layer = 0; layer < EMIT_KMER_LAYERS; layer++)
        feature_observe(NULL, scale, freq, emission_index(&site, layer, state), 1, 1);
      if (site.gc_dev != 0.0f || site.gc_contrast != 0.0f) {
        feature_observe(NULL, scale, freq, OFF_GC_DIFF + es, 1, site.gc_contrast);
        feature_observe(NULL, scale, freq, OFF_GC_LVL + es, 1, site.gc_dev);
        feature_observe(NULL, scale, freq, OFF_E3_GC + es * K3 + (int)site.kmer_emit[0][rev],
                        1, site.gc_dev);
      }
      if (is_cds(state))
        feature_observe(NULL, scale, freq, OFF_F1_WEIGHT, 1, CRF_F1_BASE_WEIGHT);
      if (!i) continue;
      int edge = path_edge(y, c.length, i);
      if (edge < 0) anno_fail("Invalid transition in validated CRF training path");
      int t_id = edge_trans[edge] - OFF_TRANSITION;
      feature_observe(NULL, scale, freq, edge_trans[edge], 1, 1);
      if (site.gc_dev != 0.0f)
        feature_observe(NULL, scale, freq, OFF_TRANS_GC + t_id, 1, site.gc_dev);
      int reverse = edge_rev[edge], kind = edge_kind[edge];
      if (kind >= 0)
        boundary_observe(&c, &site, c.offset + i, reverse, kind, 1, NULL, scale, freq);
    }
  }
  if (has_genome_window && freq) {
    /* Activate non-N candidate features on genome windows and match interaction scales to base counts. */
    for (int s = 0; s < EMIT_STATES; s++) {
      for (int k = 0; k < 64; k++) {
        int idx_e3 = OFF_E3 + s * K3 + k, idx_gc = OFF_E3_GC + s * K3 + k;
        feature_observe(NULL, scale, freq, idx_e3, 1, 1);
        feature_observe(NULL, scale, freq, idx_gc, 1, 1);
        if (scale[idx_gc] < scale[idx_e3]) scale[idx_gc] = scale[idx_e3];
      }
      feature_observe(NULL, scale, freq, OFF_GC_DIFF + s, 1, 1);
      feature_observe(NULL, scale, freq, OFF_GC_LVL + s, 1, 1);
    }
    static const int pwm_min_disp[4] = {-3, -24, -6, -3};
    static const int pwm_max_disp[4] = { 6,   3,  4,  5};
    for (int kind = 0; kind < 4; kind++) {
      for (int offset = 0; offset < BOUNDARY_KMER_COUNT; offset++)
        for (int k = 0; k < 16; k++) {
          int idx = (kind * BOUNDARY_KMER_COUNT + offset) * K2 + k;
          feature_observe(NULL, scale, freq, OFF_SPLICE_START_STOP + idx, 1, 1);
          feature_observe(NULL, scale, freq, OFF_BCROSS_D1 + idx, 1, 1);
          feature_observe(NULL, scale, freq, OFF_BCROSS_D2 + idx, 1, 1);
          double base_s = fmax(scale[OFF_SPLICE_START_STOP + idx], 16.0);
          if (scale[OFF_BCROSS_D1 + idx] < base_s) scale[OFF_BCROSS_D1 + idx] = base_s;
          if (scale[OFF_BCROSS_D2 + idx] < base_s) scale[OFF_BCROSS_D2 + idx] = base_s;
        }
      for (int offset = 0; offset < BOUNDARY_PWM_COUNT; offset++) {
        int disp = BOUNDARY_PWM_FROM + offset;
        int active_disp = disp >= pwm_min_disp[kind] && disp <= pwm_max_disp[kind];
        for (int b = 0; b < 4; b++) {
          int idx = OFF_PWM + (kind * BOUNDARY_PWM_COUNT + offset) * KNUC + b;
          if (active_disp) feature_observe(NULL, scale, freq, idx, 1, 1);
          else freq[idx] = 0;
        }
      }
      for (int offset = 0; offset < BOUNDARY_CWT_COUNT; offset++)
        for (int b = 0; b < 4; b++) {
          int idx = OFF_BPWM_CWT + (kind * BOUNDARY_CWT_COUNT + offset) * KNUC + b;
          feature_observe(NULL, scale, freq, idx, 1, 4);
        }
    }
    for (int t = 0; t < TRANSITION_FEATURES; t++) {
      feature_observe(NULL, scale, freq, OFF_TRANSITION + t, 1, 1);
      feature_observe(NULL, scale, freq, OFF_TRANS_GC + t, 1, 1);
      if (scale[OFF_TRANS_GC + t] < scale[OFF_TRANSITION + t])
        scale[OFF_TRANS_GC + t] = scale[OFF_TRANSITION + t];
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
  {"boundary_pwm", OFF_PWM, OFF_BCROSS_D1},
  {"boundary_k2_cwt_d1", OFF_BCROSS_D1, OFF_BCROSS_D2},
  {"boundary_k2_cwt_d2", OFF_BCROSS_D2, OFF_BPWM_CWT},
  {"boundary_pwm_cwt", OFF_BPWM_CWT, OFF_GC_DIFF},
  {"gc_context", OFF_GC_DIFF, OFF_TRANSITION},
  {"transition", OFF_TRANSITION, OFF_F1_WEIGHT},
  {"f1_weight", OFF_F1_WEIGHT, WEIGHTS}
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

static void viterbi_span(const double *w, const Contig *c, int from, int to,
                         double *state_vec, uint8_t *back) {
  SiteStream stream;
  site_stream_init(&stream, c, from);
  double emissions[STATES], transitions[STATES * STATES], current[HIDDEN_STATES];
  uint8_t dummy_back[HIDDEN_STATES];
  for (int pos = from; pos < to; pos++) {
    Site site = site_scores_stream(w, &stream, pos, emissions, transitions);
    viterbi_relax(site.nuc, pos, emissions, transitions, state_vec, current,
                  back ? back + (size_t)pos * HIDDEN_STATES : dummy_back);
    memcpy(state_vec, current, sizeof(current));
  }
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
  if (chunk->chunk_id == 0) {
    memset(chunk->end_vec, 0, sizeof(chunk->end_vec));
    viterbi_span(w, c, 0, chunk->c_end, chunk->end_vec, chunk->back);
    return NULL;
  }

  int w_start = chunk->c_start > CRF_WINDOW_SIZE ? chunk->c_start - CRF_WINDOW_SIZE : 0;
  double previous[HIDDEN_STATES] = {0};
  if (w_start > 0)
    for (int h = 1; h < HIDDEN_STATES; h++) previous[h] = -INFINITY;
  viterbi_span(w, c, w_start, chunk->c_start, previous, NULL);
  memcpy(chunk->start_vec, previous, sizeof(previous));
  viterbi_span(w, c, chunk->c_start, chunk->c_end, previous, chunk->back);
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

  double previous[HIDDEN_STATES] = {0};
  if (n_chunks == 1) {
    viterbi_span(w, c, 0, length, previous, back);
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
  static const char *const transition_names[TRANSITION_FEATURES] = {
    "intergenic->intergenic", "intergenic->start_codon",
    "start_codon->start_codon", "start_codon->CDS", "start_codon->stop_codon", "start_codon->5ss",
    "CDS->CDS", "CDS->stop_codon", "CDS->5ss",
    "stop_codon->stop_codon", "stop_codon->intergenic", "stop_codon->5ss",
    "5ss->5ss", "5ss->intron", "5ss->3ss",
    "intron->intron", "intron->3ss",
    "3ss->3ss", "3ss->start_codon", "3ss->CDS", "3ss->stop_codon"
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
  for (int kind = 0; kind < 4; kind++)
    for (int offset = 0; offset < BOUNDARY_CWT_COUNT; offset++)
      for (int sc = 0; sc < WAVE_COUNT; sc++) {
        int idx = bp_index(kind, offset, sc);
        if (!(float)w[idx]) continue;
        fprintf(file, "%d\tboundary_cwt_power\t%s\toffset=%+d,scale=%d\t%.9g\n",
                idx, kind_names[kind], boundary_cwt_offsets[offset], wave_sizes[sc], (float)w[idx]);
      }
  for (int kind = 0; kind < 4; kind++)
    for (int slot = 0; slot < DERIV_COUNT; slot++)
      for (int offset = 0; offset < BOUNDARY_CWT_COUNT; offset++)
        for (int sc = 0; sc < WAVE_COUNT; sc++) {
          int idx = derivative_index(kind, slot, offset, sc);
          if (!(float)w[idx]) continue;
          fprintf(file, "%d\tboundary_cwt_deriv%d\t%s\toffset=%+d,scale=%d\t%.9g\n",
                  idx, slot + 1, kind_names[kind], boundary_cwt_offsets[offset], wave_sizes[sc],
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
  for (int kind = 0; kind < 4; kind++)
    for (int offset = 0; offset < BOUNDARY_KMER_COUNT; offset++)
      for (int k = 0; k < K2; k++) {
        int idx1 = OFF_BCROSS_D1 + (kind * BOUNDARY_KMER_COUNT + offset) * K2 + k;
        int idx2 = OFF_BCROSS_D2 + (kind * BOUNDARY_KMER_COUNT + offset) * K2 + k;
        format_kmer(k, 2, kbuf);
        if ((float)w[idx1])
          fprintf(file, "%d\tboundary_k2_cwt_d1\t%s\toffset=%+d,k2=%s\t%.9g\n",
                  idx1, kind_names[kind], boundary_kmer_offsets[offset], kbuf, (float)w[idx1]);
        if ((float)w[idx2])
          fprintf(file, "%d\tboundary_k2_cwt_d2\t%s\toffset=%+d,k2=%s\t%.9g\n",
                  idx2, kind_names[kind], boundary_kmer_offsets[offset], kbuf, (float)w[idx2]);
      }
  for (int kind = 0; kind < 4; kind++)
    for (int offset = 0; offset < BOUNDARY_CWT_COUNT; offset++)
      for (int b = 0; b < KNUC; b++) {
        int idx = OFF_BPWM_CWT + (kind * BOUNDARY_CWT_COUNT + offset) * KNUC + b;
        if (!(float)w[idx]) continue;
        fprintf(file, "%d\tboundary_pwm_cwt\t%s\toffset=%+d,base=%c\t%.9g\n",
                idx, kind_names[kind], boundary_cwt_offsets[offset], "ACGTN"[b], (float)w[idx]);
      }
  for (int s = 0; s < EMIT_STATES; s++) {
    int idx_d = OFF_GC_DIFF + s, idx_l = OFF_GC_LVL + s;
    if ((float)w[idx_d])
      fprintf(file, "%d\tgc_contrast\t%s\tk3_minus_local_gc\t%.9g\n",
              idx_d, state_names[s], (float)w[idx_d]);
    if ((float)w[idx_l])
      fprintf(file, "%d\tgc_level\t%s\tlocal_gc_dev\t%.9g\n",
              idx_l, state_names[s], (float)w[idx_l]);
    for (int k = 0; k < K3; k++) {
      int idx_k = OFF_E3_GC + s * K3 + k;
      if (!(float)w[idx_k]) continue;
      format_kmer(k, 3, kbuf);
      fprintf(file, "%d\temit_k3_gc\t%s\tk3=%s\t%.9g\n",
              idx_k, state_names[s], kbuf, (float)w[idx_k]);
    }
  }
  for (int transition = 0; transition < TRANSITION_FEATURES; transition++) {
    int idx = OFF_TRANS_GC + transition;
    if (!(float)w[idx]) continue;
    fprintf(file, "%d\ttransition_gc\t%s\tlocal_gc_dev\t%.9g\n",
            idx, transition_names[transition], (float)w[idx]);
  }

  for (int transition = 0; transition < TRANSITION_FEATURES; transition++) {
    int idx = OFF_TRANSITION + transition;
    if (!(float)w[idx]) continue;
    fprintf(file, "%d\ttransition\t%s\tbias\t%.9g\n",
            idx, transition_names[transition], (float)w[idx]);
  }
  if ((float)w[OFF_F1_WEIGHT])
    fprintf(file, "%d\tf1_weight\tCDS\tlog_scale\t%.9g\n",
            OFF_F1_WEIGHT, (float)w[OFF_F1_WEIGHT]);

  if (fclose(file)) anno_fail("Cannot finish writing weights text: %s", path);
}

#endif
