#ifdef __APPLE__
#define _DARWIN_C_SOURCE
#endif
#define _POSIX_C_SOURCE 200809L

#include <complex.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "anno.h"

/* Intercept anno_fail when testing expected failure conditions (e.g., no valid CRF path). */
static int g_expect_fail = 0;
static int g_fail_caught = 0;
static char g_fail_msg[512];
static jmp_buf g_fail_jmp;
static void *g_fail_allocs[16];
static int g_fail_alloc_count = 0;

void anno_fail(const char *format, ...) {
  va_list args;
  va_start(args, format);
  vsnprintf(g_fail_msg, sizeof(g_fail_msg), format, args);
  va_end(args);
  if (g_expect_fail) {
    for (int i = 0; i < g_fail_alloc_count; i++) free(g_fail_allocs[i]);
    g_fail_alloc_count = 0;
    g_fail_caught = 1;
    longjmp(g_fail_jmp, 1);
  }
  fprintf(stderr, "Error: %s\n", g_fail_msg);
  exit(EXIT_FAILURE);
}

void *anno_alloc(size_t count, size_t size) {
  if (size && (count ? count : 1) > SIZE_MAX / size) anno_fail("Allocation size overflow");
  void *memory = calloc(count ? count : 1, size);
  if (!memory) anno_fail("Out of memory");
  if (g_expect_fail && g_fail_alloc_count < 16)
    g_fail_allocs[g_fail_alloc_count++] = memory;
  return memory;
}

#define RFFT_IMPLEMENTATION
#include "rfft.h"
#define LBFGS_IMPLEMENTATION
#include "lbfgs.h"
#define CRF_IMPLEMENTATION
#include "crf.h"

const int wave_sizes[WAVE_COUNT] = {CWT_SCALES};
#define CWT_IMPLEMENTATION
#include "cwt.h"

static void test_cwt_transform(const Wavelets *wavelets, const char *sequence, int length,
                               int threads, double *features, float *complex_values,
                               float *powers) {
  size_t value_count = (size_t)length * CWT_CHANNELS;
  float *values = complex_values;
  if (features && !values) values = anno_alloc(value_count, sizeof(*values));
  cwt_transform(wavelets, sequence, length, threads, values, powers);
  if (features)
    for (size_t i = 0; i < value_count; i++) features[i] = values[i];
  if (features && !complex_values) free(values);
}

#define cwt_transform test_cwt_transform

/* =========================================================================
 * Lightweight Unity-style Test Harness
 * ========================================================================= */

static int g_tests_run = 0;
static int g_tests_failed = 0;
static int g_current_failed = 0;
static const char *g_current_test = NULL;

#define TEST_FAIL_AT(line, ...)                                                \
  do {                                                                         \
    fprintf(stderr, "test.c:%d:%s:FAIL: ", (line), g_current_test);            \
    fprintf(stderr, __VA_ARGS__);                                              \
    fputc('\n', stderr);                                                       \
    g_current_failed = 1;                                                      \
    return;                                                                    \
  } while (0)

#define TEST_ASSERT(cond)                                                      \
  do {                                                                         \
    if (!(cond)) TEST_FAIL_AT(__LINE__, "Assertion failed: %s", #cond);        \
  } while (0)

#define TEST_ASSERT_MSG(cond, ...)                                             \
  do {                                                                         \
    if (!(cond)) TEST_FAIL_AT(__LINE__, __VA_ARGS__);                          \
  } while (0)

#define TEST_ASSERT_INT_EQ(expected, actual)                                   \
  do {                                                                         \
    long long _e = (long long)(expected), _a = (long long)(actual);            \
    if (_e != _a)                                                              \
      TEST_FAIL_AT(__LINE__, "Expected %lld, got %lld (%s != %s)", _e, _a,     \
                   #expected, #actual);                                        \
  } while (0)

#define TEST_ASSERT_DOUBLE_WITHIN(tol, expected, actual)                       \
  do {                                                                         \
    double _tol = (double)(tol);                                               \
    double _exp = (double)(expected);                                          \
    double _act = (double)(actual);                                            \
    double _diff = fabs(_act - _exp);                                          \
    if (!isfinite(_act) || !isfinite(_exp) || _diff > _tol)                    \
      TEST_FAIL_AT(__LINE__,                                                   \
                   "Expected %.12g +/- %.3e, got %.12g (diff %.3e) [%s vs %s]",\
                   _exp, _tol, _act, _diff, #expected, #actual);               \
  } while (0)

#define TEST_ASSERT_FAILS(stmt)                                                \
  do {                                                                         \
    g_expect_fail = 1;                                                         \
    g_fail_caught = 0;                                                         \
    g_fail_alloc_count = 0;                                                    \
    g_fail_msg[0] = '\0';                                                      \
    if (setjmp(g_fail_jmp) == 0) {                                             \
      stmt;                                                                    \
    }                                                                          \
    g_expect_fail = 0;                                                         \
    g_fail_alloc_count = 0;                                                    \
    if (!g_fail_caught)                                                        \
      TEST_FAIL_AT(__LINE__, "Expected anno_fail to be triggered by: %s",      \
                   #stmt);                                                     \
  } while (0)

#define RUN_TEST(fn)                                                           \
  do {                                                                         \
    g_current_test = #fn;                                                      \
    g_current_failed = 0;                                                      \
    g_expect_fail = 0;                                                         \
    g_tests_run++;                                                             \
    fn();                                                                      \
    if (g_current_failed) {                                                    \
      g_tests_failed++;                                                        \
    } else {                                                                   \
      printf("test.c:%s:PASS\n", #fn);                                         \
    }                                                                          \
  } while (0)

static uint64_t rng_state = 0x853c49e6748fea9bULL;

static void rng_seed(uint64_t seed) {
  rng_state = seed ? seed : 0x853c49e6748fea9bULL;
}

static double rng_uniform(double low, double high) {
  uint64_t x = rng_state;
  x ^= x >> 12;
  x ^= x << 25;
  x ^= x >> 27;
  rng_state = x;
  double u = (double)((x * 2685821657736338717ULL) >> 11) * (1.0 / 9007199254740992.0);
  return low + (high - low) * u;
}

/* =========================================================================
 * 1. CRF Brute-Force Reference Oracle (No Dynamic Programming)
 * ========================================================================= */

enum { BF_MAX_T = 16, BF_MAX_PATHS = 256 };

typedef struct {
  int count;
  int length;
  uint8_t paths[BF_MAX_PATHS][BF_MAX_T];
  double scores[BF_MAX_PATHS];
  double log_probs[BF_MAX_PATHS];
  double probs[BF_MAX_PATHS];
  double logz;
  double max_score;
  int best_index;
  int tie_count;
  double node_marginals[BF_MAX_T][STATES];
  double edge_marginals[BF_MAX_T][STATES][STATES];
  double expected_features[WEIGHTS];
} BruteForceCrf;

/* Check if a candidate observed state sequence y[0..T-1] is valid in the active grammar. */
static int bf_path_valid_active(const Contig *c, const uint8_t *y) {
  int T = c->length;
  int b0 = base_index(c->seq[c->offset]);
  if (b0 < 0) b0 = 4;
  int current_hidden = -1;
  for (int h = 0; h < HIDDEN_STATES; h++) {
    if (hidden_state[h] == y[0] && initial_allowed(b0, h) &&
        (T > 1 || terminal_allowed(h))) {
      current_hidden = h;
      break;
    }
  }
  if (current_hidden < 0) return 0;

  for (int t = 1; t < T; t++) {
    int base = base_index(c->seq[c->offset + t]);
    if (base < 0) base = 4;
    const HiddenEdge *edges = hidden_edges[base];
    int edge_count = hidden_edge_count[base];
    int next_hidden = -1;
    for (int e = 0; e < edge_count; e++) {
      if (edges[e].previous == current_hidden &&
          trans_b[edges[e].transition] == y[t]) {
        if (t < T - 1 || terminal_allowed(edges[e].current))
          next_hidden = edges[e].current;
        break;
      }
    }
    if (next_hidden < 0) return 0;
    current_hidden = next_hidden;
  }
  return terminal_allowed(current_hidden);
}

/* Compute S(x, y) and accumulate a * F(x, y) directly from definitions (without DP). */
static double bf_score_and_accumulate_features(const double *w, const Contig *c,
                                               const uint8_t *y, double scale,
                                               double *feat_accum) {
  double total_score = 0.0;
  for (int t = 0; t < c->length; t++) {
    Site site = site_at(c, t);
    int st = y[t];
    int e_idx = emission_index(&site, st);
    total_score += w[e_idx];
    if (feat_accum) feat_accum[e_idx] += scale;
    int backoff_idx = OFF_E3 + emit_state(st) * K3 + site.k3_emit[st >= 7];
    total_score += w[backoff_idx];
    if (feat_accum) feat_accum[backoff_idx] += scale;

    if (t > 0) {
      int prev = y[t - 1];
      int edge = pair_to_edge[prev][st];
      int rev = edge_rev[edge];
      int kind = edge_kind[edge];
      if (kind >= 0) {
        for (int off = 0; off < BOUNDARY_KMER_COUNT; off++) {
          int k_idx = boundary_kmer_index(&site, rev, kind, off);
          total_score += w[k_idx];
          if (feat_accum) feat_accum[k_idx] += scale;
        }
        for (int off = 0; off < BOUNDARY_PWM_COUNT; off++) {
          int idx = OFF_PWM + (kind * BOUNDARY_PWM_COUNT + off) * KNUC +
                    site.boundary_nuc[rev][off];
          total_score += w[idx];
          if (feat_accum) feat_accum[idx] += scale;
        }
        if (kind == 1) {
          int gpos = c->offset + t;
          for (int off = 0; off < BOUNDARY_CWT_COUNT; off++) {
            float power[WAVE_COUNT];
            double deriv[DERIV_COUNT][WAVE_COUNT];
            boundary_cwt_at(c, gpos, rev, off, power, deriv);
            for (int sc = 0; sc < WAVE_COUNT; sc++) {
              int bp_idx = OFF_BP + off * WAVE_COUNT + sc;
              total_score += w[bp_idx] * (double)power[sc];
              if (feat_accum) feat_accum[bp_idx] += scale * (double)power[sc];
              for (int slot = 0; slot < DERIV_COUNT; slot++) {
                int bd_idx = derivative_index(slot, off, sc);
                total_score += w[bd_idx] * deriv[slot][sc];
                if (feat_accum) feat_accum[bd_idx] += scale * deriv[slot][sc];
              }
            }
          }
        }
      }
    }
  }
  return total_score;
}

static void bf_enum_rec(const double *w, const Contig *c, int t, int hidden,
                        uint8_t *cur_path, BruteForceCrf *out) {
  if (t == c->length) {
    if (!terminal_allowed(hidden)) return;
    if (out->count >= BF_MAX_PATHS) anno_fail("BruteForceCrf path capacity exceeded");
    int idx = out->count++;
    memcpy(out->paths[idx], cur_path, (size_t)c->length);
    out->scores[idx] = bf_score_and_accumulate_features(w, c, cur_path, 0.0, NULL);
    return;
  }

  int base = base_index(c->seq[c->offset + t]);
  if (base < 0) base = 4;

  if (t == 0) {
    for (int h = 0; h < HIDDEN_STATES; h++) {
      if (initial_allowed(base, h) && (c->length > 1 || terminal_allowed(h))) {
        cur_path[0] = hidden_state[h];
        bf_enum_rec(w, c, 1, h, cur_path, out);
      }
    }
  } else {
    const HiddenEdge *edges = hidden_edges[base];
    int edge_count = hidden_edge_count[base];
    for (int e = 0; e < edge_count; e++) {
      if (edges[e].previous == hidden) {
        cur_path[t] = trans_b[edges[e].transition];
        bf_enum_rec(w, c, t + 1, edges[e].current, cur_path, out);
      }
    }
  }
}

/* Brute-force enumeration without DP; for tiny T and K=2, also verifies against all States^T. */
static void bf_crf_solve(const double *w, const Contig *c, BruteForceCrf *out,
                         int verify_exhaustive_states) {
  memset(out, 0, sizeof(*out));
  out->length = c->length;
  uint8_t cur_path[BF_MAX_T] = {0};
  bf_enum_rec(w, c, 0, 0, cur_path, out);

  if (verify_exhaustive_states && c->length <= 4) {
    /* Cross-check by iterating over all 13^T sequences and checking bf_path_valid_active */
    int total_tuples = 1;
    for (int t = 0; t < c->length; t++) total_tuples *= STATES;
    int valid_tuples = 0;
    for (int mask = 0; mask < total_tuples; mask++) {
      int tmp = mask;
      for (int t = c->length - 1; t >= 0; t--) {
        cur_path[t] = (uint8_t)(tmp % STATES);
        tmp /= STATES;
      }
      if (bf_path_valid_active(c, cur_path)) valid_tuples++;
    }
    if (valid_tuples != out->count)
      anno_fail("Exhaustive tuple count (%d) != grammar enumeration (%d)",
                valid_tuples, out->count);
  }

  double max_s = -INFINITY;
  out->best_index = -1;
  for (int m = 0; m < out->count; m++) {
    if (out->scores[m] > max_s) {
      max_s = out->scores[m];
      out->best_index = m;
    }
  }
  out->max_score = max_s;
  out->tie_count = 0;
  for (int m = 0; m < out->count; m++) {
    if (fabs(out->scores[m] - max_s) <= 1e-12) out->tie_count++;
  }

  double sum_exp = 0.0;
  for (int m = 0; m < out->count; m++) sum_exp += exp(out->scores[m] - max_s);
  out->logz = max_s + log(sum_exp);

  for (int m = 0; m < out->count; m++) {
    out->log_probs[m] = out->scores[m] - out->logz;
    out->probs[m] = exp(out->log_probs[m]);
    double p = out->probs[m];
    for (int t = 0; t < c->length; t++) {
      out->node_marginals[t][out->paths[m][t]] += p;
      if (t > 0)
        out->edge_marginals[t][out->paths[m][t - 1]][out->paths[m][t]] += p;
    }
    bf_score_and_accumulate_features(w, c, out->paths[m], p, out->expected_features);
  }
}

/* Extract forward-backward logZ, node marginals, and edge marginals from crf.h's DP tables. */
static void dp_extract_marginals(const double *w, const Contig *c,
                                 double *out_logz,
                                 double node_marg[BF_MAX_T][STATES],
                                 double edge_marg[BF_MAX_T][STATES][STATES],
                                 double *out_grad) {
  Work work = {0};
  memset(out_grad, 0, WEIGHTS * sizeof(double));
  double ll = contig_gradient(w, c, &work, out_grad);
  double gold_s = bf_score_and_accumulate_features(w, c, c->label + c->offset, 0.0, NULL);
  *out_logz = gold_s - ll;

  memset(node_marg, 0, sizeof(double) * BF_MAX_T * STATES);
  memset(edge_marg, 0, sizeof(double) * BF_MAX_T * STATES * STATES);

  double B[HIDDEN_STATES];
  const double *last_L = work.log_alpha + (size_t)(c->length - 1) * HIDDEN_STATES;
  for (int h = 0; h < HIDDEN_STATES; h++) B[h] = exp(last_L[h]);

  for (int pos = c->length - 1; pos >= 0; pos--) {
    const Site *site = work.sites + pos;
    const double *L = work.log_alpha + (size_t)pos * HIDDEN_STATES;
    for (int h = 0; h < HIDDEN_STATES; h++)
      node_marg[pos][hidden_state[h]] += B[h];

    if (pos > 0) {
      const double *prev_L = L - HIDDEN_STATES;
      const double *psi = work.psi + (size_t)pos * ntrans;
      double log_norm = work.log_norm[pos];
      int last = pos == c->length - 1;
      double B_prev[HIDDEN_STATES] = {0};
      const HiddenEdge *edges = hidden_edges[site->nuc];
      int edge_count = hidden_edge_count[site->nuc];
      for (int e = 0; e < edge_count; e++) {
        const HiddenEdge *step = edges + e;
        if (!(B[step->current] > 0.0) || !forward_edge_allowed(prev_L, step, last)) continue;
        double ep_prob = B[step->current] *
          exp(prev_L[step->previous] + psi[step->transition] - log_norm - L[step->current]);
        B_prev[step->previous] += ep_prob;
        int a = trans_a[step->transition];
        int b = trans_b[step->transition];
        edge_marg[pos][a][b] += ep_prob;
      }
      double sum_b = 0.0;
      for (int h = 0; h < HIDDEN_STATES; h++) sum_b += B_prev[h];
      double inv_b = sum_b > 0.0 ? 1.0 / sum_b : 0.0;
      for (int h = 0; h < HIDDEN_STATES; h++) B[h] = B_prev[h] * inv_b;
    }
  }
  work_free(&work);
}

/* Configure a K=2 unrestricted binary chain over T=1..5 so all 2^T paths are valid.
   At t=0 (with base 'A'), hidden=0 (state 0) and hidden=3 (state 1) are initial-allowed.
   At t>=1, hidden=0 (state 0) and hidden=9 (state 3, context 1) are terminal-allowed.
   Transitions (0->0), (0->coding), (coding->0), (coding->coding) get distinct boundary kinds
   so all 4 transition entries A[i][j] are independently parameterized. */
static void setup_k2_binary_crf(void) {
  crf_init();
  ntrans = 0;
  memset(pair_to_edge, -1, sizeof(pair_to_edge));

  const int pairs[6][4] = {
    /* from, to, rev, kind */
    {0, 0, 0, 0},
    {0, 3, 0, 2},
    {1, 0, 0, 3},
    {1, 3, 0, 1},
    {3, 0, 0, 3},
    {3, 3, 0, 1},
  };
  for (int i = 0; i < 6; i++) {
    int a = pairs[i][0], b = pairs[i][1];
    pair_to_edge[a][b] = ntrans;
    trans_a[ntrans] = a;
    trans_b[ntrans] = b;
    edge_rev[ntrans] = pairs[i][2];
    edge_kind[ntrans] = pairs[i][3];
    ntrans++;
  }

  memset(hidden_edge_count, 0, sizeof(hidden_edge_count));
  for (int base = 0; base < KNUC; base++) {
    int c = 0;
    /* From t=0 states {h=0 (st 0), h=3 (st 1)} to t>=1 states {h=0 (st 0), h=9 (st 3)} */
    hidden_edges[base][c++] = (HiddenEdge){0, 0, (uint16_t)pair_to_edge[0][0]};
    hidden_edges[base][c++] = (HiddenEdge){0, 9, (uint16_t)pair_to_edge[0][3]};
    hidden_edges[base][c++] = (HiddenEdge){3, 0, (uint16_t)pair_to_edge[1][0]};
    hidden_edges[base][c++] = (HiddenEdge){3, 9, (uint16_t)pair_to_edge[1][3]};
    /* From t>=1 state h=9 (st 3) to {h=0 (st 0), h=9 (st 3)} */
    hidden_edges[base][c++] = (HiddenEdge){9, 0, (uint16_t)pair_to_edge[3][0]};
    hidden_edges[base][c++] = (HiddenEdge){9, 9, (uint16_t)pair_to_edge[3][3]};
    hidden_edge_count[base] = c;
  }
}

static Contig make_test_contig(const char *full_seq, int offset, int length,
                               const uint8_t *window_labels) {
  int glen = (int)strlen(full_seq);
  Contig c = {0};
  c.name = "test_contig";
  c.seq = (char *)full_seq;
  c.genome_length = glen;
  c.offset = offset;
  c.length = length;
  c.cwt_power = anno_alloc((size_t)glen * WAVE_COUNT, sizeof(float));
  c.label = anno_alloc((size_t)glen, sizeof(uint8_t));
  for (int i = 0; i < glen * WAVE_COUNT; i++)
    c.cwt_power[i] = (float)(0.15 + 0.07 * sin(0.37 * (i + 1)) + 0.04 * cos(0.19 * (i + 2)));
  if (window_labels)
    memcpy(c.label + offset, window_labels, (size_t)length);
  return c;
}

static void free_test_contig(Contig *c) {
  free(c->cwt_power);
  free(c->label);
}

/* =========================================================================
 * Part 1 Tests: CRF Module
 * ========================================================================= */

static void test_crf_context_coordinates_and_stream(void) {
  crf_init();
  const char *seq = "ACgtNATGCCATTAACCGGTANCGTTAGCACGTACGTACGTACGT";
  Contig c = make_test_contig(seq, 2, (int)strlen(seq) - 4, NULL);
  double *w = anno_alloc(WEIGHTS, sizeof(*w));
  for (int f = 0; f < WEIGHTS; f++) w[f] = sin(f * 0.17) * 0.1;
  for (int start = 0; start < c.length; start += 7) {
    SiteStream stream;
    site_stream_init(&stream, &c, start);
    for (int i = start; i < c.length; i++) {
      double emissions[STATES], psi[STATES * STATES];
      Site actual = site_scores_stream(w, &stream, i, emissions, psi);
      Site expected = site_at(&c, i);
      TEST_ASSERT_INT_EQ(expected.nuc, actual.nuc);
      for (int rev = 0; rev < 2; rev++) {
        TEST_ASSERT_INT_EQ(expected.k6_emit[rev], actual.k6_emit[rev]);
        TEST_ASSERT_INT_EQ(expected.k3_emit[rev], actual.k3_emit[rev]);
        for (int off = 0; off < BOUNDARY_KMER_COUNT; off++)
          TEST_ASSERT_INT_EQ(expected.boundary_k2[rev][off], actual.boundary_k2[rev][off]);
        for (int off = 0; off < BOUNDARY_PWM_COUNT; off++) {
          int d = BOUNDARY_PWM_FROM + off;
          int g = rev ? c.offset + i - d - 1 : c.offset + i + d;
          int b = g >= 0 && g < c.genome_length ? base_index(seq[g]) : -1;
          int code = b < 0 ? 4 : rev ? 3 - b : b;
          TEST_ASSERT_INT_EQ(code, actual.boundary_nuc[rev][off]);
        }
      }
      for (int st = 0; st < STATES; st++)
        TEST_ASSERT_DOUBLE_WITHIN(1e-12, emission_score(w, &expected, st), emissions[st]);
    }
  }
  free(w);
  free_test_contig(&c);
}

static void test_crf_reverse_complement_emissions(void) {
  crf_init();
  const char *seq = "ACGTAATGCGGCTANACGTTAGCCATAGC";
  int length = (int)strlen(seq);
  char reverse[64];
  for (int i = 0; i < length; i++) {
    int b = base_index(seq[length - i - 1]);
    reverse[i] = b < 0 ? 'N' : "TGCA"[b];
  }
  reverse[length] = '\0';
  Contig forward = make_test_contig(seq, 0, length, NULL);
  Contig backward = make_test_contig(reverse, 0, length, NULL);
  double *w = anno_alloc(WEIGHTS, sizeof(*w));
  for (int f = 0; f < WEIGHTS; f++) w[f] = sin(f * 0.13);
  for (int i = 0; i < length; i++) {
    Site a = site_at(&forward, i);
    Site b = site_at(&backward, length - i - 1);
    TEST_ASSERT_INT_EQ(a.k6_emit[0], b.k6_emit[1]);
    TEST_ASSERT_INT_EQ(a.k3_emit[0], b.k3_emit[1]);
    for (int st = 1; st <= 3; st++)
      TEST_ASSERT_DOUBLE_WITHIN(0, emission_score(w, &a, st), emission_score(w, &b, st + 6));
  }
  free(w);
  free_test_contig(&forward);
  free_test_contig(&backward);
}

static void test_crf_k2_exhaustive_t1_to_t5(void) {
  setup_k2_binary_crf();
  rng_seed(1001);

  double *w = anno_alloc(WEIGHTS, sizeof(double));
  for (int i = 0; i < WEIGHTS; i++) w[i] = rng_uniform(-0.6, 0.6);

  const char *full_seq =
    "ACGTACGTACGTACGTACGTACGTACGTACGT"
    "AAACGTACGTACGTACGTACGTACGTACGTAC";

  for (int T = 1; T <= 5; T++) {
    uint8_t gold[BF_MAX_T] = {0};
    if (T >= 2) gold[0] = 1;
    if (T >= 3) gold[T - 1] = 3;
    /* Ensure full_seq[offset] == 'A' so both h=0 and h=3 are initial-allowed */
    Contig c = make_test_contig(full_seq, 32, T, gold);

    BruteForceCrf bf;
    bf_crf_solve(w, &c, &bf, 1);
    int expected_paths = (T == 1) ? 1 : (1 << T);
    TEST_ASSERT_INT_EQ(expected_paths, bf.count);

    double dp_logz = 0.0;
    double dp_node[BF_MAX_T][STATES];
    double dp_edge[BF_MAX_T][STATES][STATES];
    double dp_grad[WEIGHTS];
    dp_extract_marginals(w, &c, &dp_logz, dp_node, dp_edge, dp_grad);

    /* 1. logZ matches brute-force log-sum-exp */
    TEST_ASSERT_DOUBLE_WITHIN(1e-11, bf.logz, dp_logz);

    /* 2. Label sequence log-probabilities match S(x, y) - logZ and sum to 1 */
    double prob_sum = 0.0;
    for (int m = 0; m < bf.count; m++) {
      Contig cy = c;
      uint8_t label_buf[128] = {0};
      memcpy(label_buf + c.offset, bf.paths[m], (size_t)T);
      cy.label = label_buf;
      Work wk = {0};
      double grad_tmp[WEIGHTS] = {0};
      double ll = contig_gradient(w, &cy, &wk, grad_tmp);
      work_free(&wk);
      TEST_ASSERT_DOUBLE_WITHIN(1e-11, bf.scores[m] - bf.logz, ll);
      prob_sum += exp(ll);
    }
    TEST_ASSERT_DOUBLE_WITHIN(1e-11, 1.0, prob_sum);

    /* 3. Viterbi path score matches brute-force maximum */
    uint8_t vpath[BF_MAX_T] = {0};
    crf_decode(w, &c, vpath);
    TEST_ASSERT(bf_path_valid_active(&c, vpath));
    double vscore = bf_score_and_accumulate_features(w, &c, vpath, 0.0, NULL);
    TEST_ASSERT_DOUBLE_WITHIN(1e-11, bf.max_score, vscore);
    if (bf.tie_count == 1) {
      for (int t = 0; t < T; t++)
        TEST_ASSERT_INT_EQ(bf.paths[bf.best_index][t], vpath[t]);
    }

    /* 4. Node marginals match brute-force */
    for (int t = 0; t < T; t++) {
      double marg_sum = 0.0;
      for (int s = 0; s < STATES; s++) {
        TEST_ASSERT_DOUBLE_WITHIN(1e-11, bf.node_marginals[t][s], dp_node[t][s]);
        marg_sum += dp_node[t][s];
      }
      TEST_ASSERT_DOUBLE_WITHIN(1e-11, 1.0, marg_sum);
    }

    /* 5. Edge marginals match brute-force and marginalize to node marginals */
    for (int t = 1; t < T; t++) {
      double edge_sum = 0.0;
      for (int a = 0; a < STATES; a++) {
        double row_sum = 0.0;
        for (int b = 0; b < STATES; b++) {
          TEST_ASSERT_DOUBLE_WITHIN(1e-11, bf.edge_marginals[t][a][b], dp_edge[t][a][b]);
          row_sum += dp_edge[t][a][b];
          edge_sum += dp_edge[t][a][b];
        }
        TEST_ASSERT_DOUBLE_WITHIN(1e-11, dp_node[t - 1][a], row_sum);
      }
      TEST_ASSERT_DOUBLE_WITHIN(1e-11, 1.0, edge_sum);
    }

    /* 6. Gradient matches brute-force feature expectations (F(x, y*) - E[F(x, y)]) */
    double gold_features[WEIGHTS] = {0};
    bf_score_and_accumulate_features(w, &c, gold, 1.0, gold_features);
    for (int f = 0; f < WEIGHTS; f++) {
      double expected_ll_grad = gold_features[f] - bf.expected_features[f];
      TEST_ASSERT_DOUBLE_WITHIN(1e-10, expected_ll_grad, dp_grad[f]);
    }

    free_test_contig(&c);
  }

  free(w);
  crf_init();
}

static void test_crf_grammar_exhaustive_both_strands(void) {
  crf_init();
  rng_seed(2026);

  double *w = anno_alloc(WEIGHTS, sizeof(double));
  for (int i = 0; i < WEIGHTS; i++) w[i] = rng_uniform(-0.35, 0.35);

  /* Case A: Plus-strand sequence with start codon, introns (5ss & 3ss + CWT), stop codon.
     Window is "ATGTTGTAA" (T=9) embedded at offset 30 inside a 72-bp contig so all CWT
     offsets (-27..+2) and derivative kernels are in-bounds and active, and introns in all
     3 reading frames (states 4, 5, 6) are valid. */
  const char *plus_genome =
    "ACGTACGTACGTACGTACGTACGTACGTAC" /* 30 bp upstream */
    "ATGTTGTAA"                      /* 9 bp window: 000000000, 123123123, 155523123, 126663123, 123444123, ... */
    "ACGTACGTACGTACGTACGTACGTACGTACGT";
  const uint8_t plus_gold[9] = {1, 2, 3, 4, 4, 4, 1, 2, 3}; /* includes 5ss (3->4) and 3ss (4->1) */
  Contig c_plus = make_test_contig(plus_genome, 30, 9, plus_gold);
  TEST_ASSERT(gold_path_valid(&c_plus));

  BruteForceCrf bf_plus;
  bf_crf_solve(w, &c_plus, &bf_plus, 0);
  TEST_ASSERT(bf_plus.count >= 6);

  double dp_logz = 0.0;
  double dp_node[BF_MAX_T][STATES];
  double dp_edge[BF_MAX_T][STATES][STATES];
  double dp_grad[WEIGHTS];
  dp_extract_marginals(w, &c_plus, &dp_logz, dp_node, dp_edge, dp_grad);

  TEST_ASSERT_DOUBLE_WITHIN(1e-11, bf_plus.logz, dp_logz);

  uint8_t vpath[BF_MAX_T] = {0};
  crf_decode(w, &c_plus, vpath);
  TEST_ASSERT(bf_path_valid_active(&c_plus, vpath));
  double vscore = bf_score_and_accumulate_features(w, &c_plus, vpath, 0.0, NULL);
  TEST_ASSERT_DOUBLE_WITHIN(1e-11, bf_plus.max_score, vscore);

  for (int t = 0; t < c_plus.length; t++) {
    for (int s = 0; s < STATES; s++)
      TEST_ASSERT_DOUBLE_WITHIN(1e-11, bf_plus.node_marginals[t][s], dp_node[t][s]);
    if (t > 0) {
      for (int a = 0; a < STATES; a++)
        for (int b = 0; b < STATES; b++)
          TEST_ASSERT_DOUBLE_WITHIN(1e-11, bf_plus.edge_marginals[t][a][b], dp_edge[t][a][b]);
    }
  }

  double gold_feat[WEIGHTS] = {0};
  bf_score_and_accumulate_features(w, &c_plus, plus_gold, 1.0, gold_feat);
  for (int f = 0; f < WEIGHTS; f++) {
    double expected_grad = gold_feat[f] - bf_plus.expected_features[f];
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, expected_grad, dp_grad[f]);
  }

  /* Verify finite differences against dp_grad across active emission, boundary_k2, CWT power, and CWT deriv weights */
  const double eps = 1e-6;
  int checked_groups[6] = {0};
  for (int f = 0; f < WEIGHTS; f++) {
    if (fabs(dp_grad[f]) < 1e-7) continue;
    int group = (f < OFF_SPLICE_START_STOP) ? 0 :
                (f < OFF_BP) ? 1 :
                (f < OFF_BD) ? 2 :
                (f < OFF_E3) ? 3 :
                (f < OFF_PWM) ? 4 : 5;
    if (checked_groups[group] >= 4) continue;
    checked_groups[group]++;

    double orig = w[f];
    Work wk = {0};
    double dummy[WEIGHTS];
    w[f] = orig + eps;
    double ll_p = contig_gradient(w, &c_plus, &wk, dummy);
    w[f] = orig - eps;
    double ll_m = contig_gradient(w, &c_plus, &wk, dummy);
    w[f] = orig;
    work_free(&wk);

    double fd = (ll_p - ll_m) / (2.0 * eps);
    TEST_ASSERT_DOUBLE_WITHIN(1e-7, fd, dp_grad[f]);
  }
  for (int g = 0; g < 6; g++) TEST_ASSERT(checked_groups[g] > 0);

  /* Case B: Minus-strand sequence with reverse-complement stop ("TTA"), intron (10..12), and start ("CAT") */
  const char *minus_genome =
    "ACGTACGTACGTACGTACGTACGTACGTAC"
    "TTATCACAT"
    "ACGTACGTACGTACGTACGTACGTACGTACGT";
  const uint8_t minus_gold[9] = {9, 8, 7, 12, 12, 12, 9, 8, 7};
  Contig c_minus = make_test_contig(minus_genome, 30, 9, minus_gold);
  TEST_ASSERT(gold_path_valid(&c_minus));

  BruteForceCrf bf_minus;
  bf_crf_solve(w, &c_minus, &bf_minus, 0);
  TEST_ASSERT(bf_minus.count >= 6);

  dp_extract_marginals(w, &c_minus, &dp_logz, dp_node, dp_edge, dp_grad);
  TEST_ASSERT_DOUBLE_WITHIN(1e-11, bf_minus.logz, dp_logz);

  memset(gold_feat, 0, sizeof(gold_feat));
  bf_score_and_accumulate_features(w, &c_minus, minus_gold, 1.0, gold_feat);
  for (int f = 0; f < WEIGHTS; f++)
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, gold_feat[f] - bf_minus.expected_features[f], dp_grad[f]);

  free_test_contig(&c_plus);
  free_test_contig(&c_minus);
  free(w);
}

static void test_crf_zero_scores(void) {
  double *w = anno_alloc(WEIGHTS, sizeof(double)); /* all zeros */
  const char *seq = "ACGTACGTACGTACGTACGTACGTACGTACGTAAAAAAAA";

  /* 1. K=2 unrestricted: logZ = T * log(2), path prob = 2^-T, node marg = 0.5, edge marg = 0.25 */
  setup_k2_binary_crf();
  for (int T = 2; T <= 5; T++) {
    uint8_t gold[BF_MAX_T] = {0};
    Contig c = make_test_contig(seq, 32, T, gold);
    double dp_logz = 0.0;
    double dp_node[BF_MAX_T][STATES];
    double dp_edge[BF_MAX_T][STATES][STATES];
    double dp_grad[WEIGHTS];
    dp_extract_marginals(w, &c, &dp_logz, dp_node, dp_edge, dp_grad);

    TEST_ASSERT_DOUBLE_WITHIN(1e-12, (double)T * log(2.0), dp_logz);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.5, dp_node[0][0]);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.5, dp_node[0][1]);
    for (int t = 1; t < T; t++) {
      TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.5, dp_node[t][0]);
      TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.5, dp_node[t][3]);
      int prev_coding = (t == 1) ? 1 : 3;
      TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.25, dp_edge[t][0][0]);
      TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.25, dp_edge[t][0][3]);
      TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.25, dp_edge[t][prev_coding][0]);
      TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.25, dp_edge[t][prev_coding][3]);
    }
    free_test_contig(&c);
  }

  /* 2. Full 13-state grammar with all-zero weights: logZ = log(#valid_paths), uniform 1/N */
  crf_init();
  const uint8_t gold6[6] = {1, 2, 3, 1, 2, 3};
  Contig c6 = make_test_contig("ATGTAA", 0, 6, gold6);
  BruteForceCrf bf6;
  bf_crf_solve(w, &c6, &bf6, 0);
  double dp_logz = 0.0;
  double dp_node[BF_MAX_T][STATES];
  double dp_edge[BF_MAX_T][STATES][STATES];
  double dp_grad[WEIGHTS];
  dp_extract_marginals(w, &c6, &dp_logz, dp_node, dp_edge, dp_grad);
  TEST_ASSERT_DOUBLE_WITHIN(1e-12, log((double)bf6.count), dp_logz);
  for (int m = 0; m < bf6.count; m++)
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 1.0 / (double)bf6.count, bf6.probs[m]);

  free_test_contig(&c6);
  free(w);
}

static void test_crf_length_one(void) {
  crf_init();
  double *w1 = anno_alloc(WEIGHTS, sizeof(double));
  double *w2 = anno_alloc(WEIGHTS, sizeof(double));
  rng_seed(303);
  for (int i = 0; i < WEIGHTS; i++) {
    w1[i] = rng_uniform(-1.0, 1.0);
    /* w2 has identical emission weights, but wildly different transition/boundary weights */
    w2[i] = (i < OFF_SPLICE_START_STOP || (i >= OFF_E3 && i < OFF_PWM))
      ? w1[i] : rng_uniform(-50.0, 50.0);
  }

  const uint8_t gold1[1] = {0};
  Contig c = make_test_contig("ACGTACGT", 3, 1, gold1);

  double logz1 = 0.0, logz2 = 0.0;
  double node1[BF_MAX_T][STATES], node2[BF_MAX_T][STATES];
  double edge1[BF_MAX_T][STATES][STATES], edge2[BF_MAX_T][STATES][STATES];
  double grad1[WEIGHTS], grad2[WEIGHTS];
  dp_extract_marginals(w1, &c, &logz1, node1, edge1, grad1);
  dp_extract_marginals(w2, &c, &logz2, node2, edge2, grad2);

  Site s0 = site_at(&c, 0);
  double expected_e0 = emission_score(w1, &s0, 0);
  TEST_ASSERT_DOUBLE_WITHIN(1e-12, expected_e0, logz1);
  TEST_ASSERT_DOUBLE_WITHIN(1e-12, logz1, logz2);
  TEST_ASSERT_DOUBLE_WITHIN(1e-12, 1.0, node1[0][0]);

  for (int f = 0; f < WEIGHTS; f++) {
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.0, grad1[f]);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.0, grad2[f]);
  }

  uint8_t path1[1] = {99}, path2[1] = {99};
  crf_decode(w1, &c, path1);
  crf_decode(w2, &c, path2);
  TEST_ASSERT_INT_EQ(0, path1[0]);
  TEST_ASSERT_INT_EQ(0, path2[0]);

  free_test_contig(&c);
  free(w1);
  free(w2);
}

static void test_crf_asymmetric_transitions(void) {
  setup_k2_binary_crf();
  double *w = anno_alloc(WEIGHTS, sizeof(double));

  /* Set 0->3 (kind 2) = +1.5 per offset (total +12.0), 3->3 (kind 1) = +0.5 per offset (+4.0),
     while 1/3->0 (kind 3) = -1.5 per offset (total -12.0). */
  for (int k = 0; k < K2; k++) {
    for (int off = 0; off < BOUNDARY_KMER_COUNT; off++) {
      w[OFF_SPLICE_START_STOP + (1 * BOUNDARY_KMER_COUNT + off) * K2 + k] = 0.5;
      w[OFF_SPLICE_START_STOP + (2 * BOUNDARY_KMER_COUNT + off) * K2 + k] = 1.5;
      w[OFF_SPLICE_START_STOP + (3 * BOUNDARY_KMER_COUNT + off) * K2 + k] = -1.5;
    }
  }

  const uint8_t gold[3] = {0, 3, 3};
  Contig c = make_test_contig("AAAAAAAA", 0, 3, gold);
  BruteForceCrf bf;
  bf_crf_solve(w, &c, &bf, 0);

  double dp_logz = 0.0;
  double dp_node[BF_MAX_T][STATES];
  double dp_edge[BF_MAX_T][STATES][STATES];
  double dp_grad[WEIGHTS];
  dp_extract_marginals(w, &c, &dp_logz, dp_node, dp_edge, dp_grad);

  TEST_ASSERT_DOUBLE_WITHIN(1e-11, bf.logz, dp_logz);
  /* 0->3 edge marginal at t=1 must be much larger than 1->0 edge marginal at t=1 */
  TEST_ASSERT(dp_edge[1][0][3] > 0.95);
  TEST_ASSERT(dp_edge[1][1][0] < 1e-5);
  TEST_ASSERT_DOUBLE_WITHIN(1e-11, bf.edge_marginals[1][0][3], dp_edge[1][0][3]);
  TEST_ASSERT_DOUBLE_WITHIN(1e-11, bf.edge_marginals[1][1][0], dp_edge[1][1][0]);

  uint8_t vpath[3] = {0};
  crf_decode(w, &c, vpath);
  TEST_ASSERT_INT_EQ(0, vpath[0]);
  TEST_ASSERT_INT_EQ(3, vpath[1]);
  TEST_ASSERT_INT_EQ(3, vpath[2]);

  free_test_contig(&c);
  free(w);
  crf_init();
}

static void test_crf_tied_optimal_paths(void) {
  crf_init();
  double *w = anno_alloc(WEIGHTS, sizeof(double));
  /* With all weights 0 on "ATGTAA", paths 000000 and 123123 tie at score 0.0. */
  const uint8_t gold[6] = {0, 0, 0, 0, 0, 0};
  Contig c = make_test_contig("ATGTAA", 0, 6, gold);

  BruteForceCrf bf;
  bf_crf_solve(w, &c, &bf, 0);
  TEST_ASSERT_INT_EQ(2, bf.count);
  TEST_ASSERT_INT_EQ(2, bf.tie_count);

  uint8_t vpath[6] = {0};
  crf_decode(w, &c, vpath);
  TEST_ASSERT(bf_path_valid_active(&c, vpath));
  double vscore = bf_score_and_accumulate_features(w, &c, vpath, 0.0, NULL);
  TEST_ASSERT_DOUBLE_WITHIN(1e-12, bf.max_score, vscore);
  /* Tie-break rule in crf_decode uses strict `>` starting from best=0, so state 0 wins ties */
  for (int t = 0; t < 6; t++) TEST_ASSERT_INT_EQ(0, vpath[t]);

  /* Now tie two non-zero intron paths in "ATGAAATAA" while penalizing state 0 */
  for (int k = 0; k < K6; k++) w[OFF_E6 + 0 * K6 + k] = -10.0;
  /* Make intron states 4, 5, 6 and CDS states 1, 2, 3 have score 0 so multiple coding/intron paths tie */
  Contig c9 = make_test_contig("ATGAAATAA", 0, 9, NULL);
  BruteForceCrf bf9;
  bf_crf_solve(w, &c9, &bf9, 0);
  TEST_ASSERT(bf9.tie_count >= 2);

  uint8_t vpath9[9] = {0};
  crf_decode(w, &c9, vpath9);
  TEST_ASSERT(bf_path_valid_active(&c9, vpath9));
  double vscore9 = bf_score_and_accumulate_features(w, &c9, vpath9, 0.0, NULL);
  TEST_ASSERT_DOUBLE_WITHIN(1e-12, bf9.max_score, vscore9);

  free_test_contig(&c);
  free_test_contig(&c9);
  free(w);
}

static void test_crf_forbidden_transitions(void) {
  crf_init();
  double *w = anno_alloc(WEIGHTS, sizeof(double));
  /* Strongly reward CDS states 1, 2, 3 and penalize intergenic 0 */
  for (int k = 0; k < K6; k++) {
    w[OFF_E6 + 0 * K6 + k] = -5.0;
    w[OFF_E6 + 1 * K6 + k] = 20.0;
    w[OFF_E6 + 2 * K6 + k] = 20.0;
    w[OFF_E6 + 3 * K6 + k] = 20.0;
  }

  /* 1. Sequence with in-frame stop codon "ATGTAGTAA":
     1,2,3,1,2,3,1,2,3 is forbidden because codon 2 is "TAG" (stop codon cannot continue in CDS). */
  const uint8_t forbidden_inframe_stop[9] = {1, 2, 3, 1, 2, 3, 1, 2, 3};
  Contig c_stop = make_test_contig("ATGTAGTAA", 0, 9, forbidden_inframe_stop);
  TEST_ASSERT(!gold_path_valid(&c_stop));

  /* make_windows must skip invalid gold paths */
  CrfWindow *windows = NULL;
  size_t skipped = 0;
  size_t win_count = make_windows(&c_stop, 1, &windows, &skipped);
  TEST_ASSERT_INT_EQ(0, win_count);
  TEST_ASSERT_INT_EQ(1, skipped);
  free(windows);

  /* crf_decode on "ATGTAGTAA" must still return a grammar-valid path (never the forbidden one) */
  uint8_t decoded[9] = {0};
  crf_decode(w, &c_stop, decoded);
  memcpy(c_stop.label, decoded, 9);
  TEST_ASSERT(gold_path_valid(&c_stop));
  int same_as_forbidden = 1;
  for (int i = 0; i < 9; i++)
    if (decoded[i] != forbidden_inframe_stop[i]) same_as_forbidden = 0;
  TEST_ASSERT(!same_as_forbidden);

  /* 2. Sequence with 'N' inside CDS ("ATGNAGTAA"): 'N' cannot be decoded as CDS 1, 2, 3 */
  Contig c_n = make_test_contig("ATGNAGTAA", 0, 9, NULL);
  crf_decode(w, &c_n, decoded);
  memcpy(c_n.label, decoded, 9);
  TEST_ASSERT(gold_path_valid(&c_n));
  TEST_ASSERT(decoded[3] == 0 || (decoded[3] >= 4 && decoded[3] <= 6));

  free_test_contig(&c_stop);
  free_test_contig(&c_n);
  free(w);
}

static void test_crf_no_valid_path_and_nonfinite(void) {
  crf_init();
  double *w = anno_alloc(WEIGHTS, sizeof(double));
  const uint8_t gold[6] = {1, 2, 3, 1, 2, 3};
  Contig c = make_test_contig("ATGTAA", 0, 6, gold);

  /* 1. Make all hidden edges empty (no valid path of length 6 exists) */
  memset(hidden_edge_count, 0, sizeof(hidden_edge_count));
  Work wk = {0};
  work_ensure(&wk, c.length);
  double grad[WEIGHTS] = {0};
  TEST_ASSERT_FAILS(contig_gradient(w, &c, &wk, grad));
  for (int f = 0; f < WEIGHTS; f++) TEST_ASSERT(!isnan(grad[f]));

  uint8_t path[6] = {0};
  TEST_ASSERT_FAILS(crf_decode(w, &c, path));

  /* Restore grammar and test non-finite weights (NaN and -INFINITY) */
  crf_init();
  Site s0 = site_at(&c, 0);
  w[emission_index(&s0, 0)] = NAN;
  TEST_ASSERT_FAILS(contig_gradient(w, &c, &wk, grad));
  TEST_ASSERT_FAILS(crf_decode(w, &c, path));
  work_free(&wk);

  free_test_contig(&c);
  free(w);
}

static void test_crf_extreme_scores_numerical_stability(void) {
  crf_init();
  double *w = anno_alloc(WEIGHTS, sizeof(double));
  const uint8_t gold[6] = {1, 2, 3, 1, 2, 3};
  Contig c = make_test_contig("ATGTAA", 0, 6, gold);

  /* Test both large positive emission weights (+160 per base -> path scores ~ +960, where
     naive exp(+960) overflows IEEE-754 double) and large negative emission weights (-160 per
     base -> path scores ~ -960, where naive exp(-960) underflows to 0.0). */
  const double scales[2] = {160.0, -160.0};
  for (int s = 0; s < 2; s++) {
    rng_seed(555 + (uint64_t)s);
    memset(w, 0, WEIGHTS * sizeof(double));
    for (int i = 0; i < OFF_SPLICE_START_STOP; i++)
      w[i] = scales[s] + rng_uniform(-2.0, 2.0);
    for (int i = OFF_SPLICE_START_STOP; i < WEIGHTS; i++)
      w[i] = rng_uniform(-0.5, 0.5);

    BruteForceCrf bf;
    bf_crf_solve(w, &c, &bf, 0);

    double dp_logz = 0.0;
    double dp_node[BF_MAX_T][STATES];
    double dp_edge[BF_MAX_T][STATES][STATES];
    double dp_grad[WEIGHTS];
    dp_extract_marginals(w, &c, &dp_logz, dp_node, dp_edge, dp_grad);

    TEST_ASSERT(isfinite(dp_logz));
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, bf.logz, dp_logz);

    for (int t = 0; t < c.length; t++) {
      for (int st = 0; st < STATES; st++) {
        TEST_ASSERT(isfinite(dp_node[t][st]));
        TEST_ASSERT_DOUBLE_WITHIN(1e-10, bf.node_marginals[t][st], dp_node[t][st]);
      }
    }

    uint8_t vpath[6] = {0};
    crf_decode(w, &c, vpath);
    double vscore = bf_score_and_accumulate_features(w, &c, vpath, 0.0, NULL);
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, bf.max_score, vscore);
  }

  free_test_contig(&c);
  free(w);
}

static void test_crf_shift_invariance(void) {
  crf_init();
  rng_seed(777);
  double *w1 = anno_alloc(WEIGHTS, sizeof(double));
  double *w2 = anno_alloc(WEIGHTS, sizeof(double));
  for (int i = 0; i < WEIGHTS; i++) w1[i] = w2[i] = rng_uniform(-0.5, 0.5);

  /* Keep all forward and reverse-complement 6-mer contexts unique so shared emission
     parameters do not couple the position-specific score shifts. */
  const char *seq = "CAGGAACTTGCGTAC";
  const uint8_t gold[9] = {1, 2, 3, 1, 2, 3, 1, 2, 3};
  Contig c = make_test_contig(seq, 3, 9, gold);

  /* Verify that all 9 positions have distinct forward & reverse 6-mer codes */
  double shift_sum = 0.0;
  for (int t = 0; t < c.length; t++) {
    Site st = site_at(&c, t);
    TEST_ASSERT(st.k6_emit[0] < 4096);
    TEST_ASSERT(st.k6_emit[1] < 4096);
    for (int previous = 0; previous < t; previous++) {
      Site prior = site_at(&c, previous);
      TEST_ASSERT(st.k6_emit[0] != prior.k6_emit[0]);
      TEST_ASSERT(st.k6_emit[1] != prior.k6_emit[1]);
    }
    double ct = (t - 4) * 3.25 + 1.1;
    shift_sum += ct;
    for (int state = 0; state < STATES; state++)
      w2[emission_index(&st, state)] += ct;
  }

  double logz1 = 0.0, logz2 = 0.0;
  double node1[BF_MAX_T][STATES], node2[BF_MAX_T][STATES];
  double edge1[BF_MAX_T][STATES][STATES], edge2[BF_MAX_T][STATES][STATES];
  double grad1[WEIGHTS], grad2[WEIGHTS];
  dp_extract_marginals(w1, &c, &logz1, node1, edge1, grad1);
  dp_extract_marginals(w2, &c, &logz2, node2, edge2, grad2);

  /* logZ must shift by exactly sum_t c_t */
  TEST_ASSERT_DOUBLE_WITHIN(1e-11, logz1 + shift_sum, logz2);

  /* Node marginals, edge marginals, gradients, and Viterbi path must be invariant */
  for (int t = 0; t < c.length; t++) {
    for (int s = 0; s < STATES; s++)
      TEST_ASSERT_DOUBLE_WITHIN(1e-11, node1[t][s], node2[t][s]);
    if (t > 0) {
      for (int a = 0; a < STATES; a++)
        for (int b = 0; b < STATES; b++)
          TEST_ASSERT_DOUBLE_WITHIN(1e-11, edge1[t][a][b], edge2[t][a][b]);
    }
  }
  for (int f = 0; f < WEIGHTS; f++)
    TEST_ASSERT_DOUBLE_WITHIN(1e-11, grad1[f], grad2[f]);

  uint8_t vpath1[9] = {0}, vpath2[9] = {0};
  crf_decode(w1, &c, vpath1);
  crf_decode(w2, &c, vpath2);
  for (int t = 0; t < c.length; t++)
    TEST_ASSERT_INT_EQ(vpath1[t], vpath2[t]);

  free_test_contig(&c);
  free(w1);
  free(w2);
}

static void test_crf_objective_sign_regularization_and_scaling(void) {
  crf_init();
  rng_seed(909);

  const uint8_t gold1[9] = {1, 2, 3, 4, 4, 4, 1, 2, 3};
  const uint8_t gold2[6] = {9, 8, 7, 9, 8, 7};
  Contig contigs[2];
  contigs[0] = make_test_contig("ATGAAATAA", 0, 9, gold1);
  contigs[1] = make_test_contig("TTACAT", 0, 6, gold2);

  CrfWindow *windows = NULL;
  size_t skipped = 0;
  size_t window_count = make_windows(contigs, 2, &windows, &skipped);
  TEST_ASSERT_INT_EQ(2, window_count);
  TEST_ASSERT_INT_EQ(0, skipped);

  CrfWorker worker = {
    .train = NULL,
    .thread_id = 0,
    .gradient = anno_alloc(WEIGHTS, sizeof(double))
  };
  CrfObjective obj = {
    .train = {.windows = windows, .count = window_count, .skipped = skipped, .threads = 1},
    .workers = &worker,
    .threads = NULL,
    .scale = anno_alloc(WEIGHTS, sizeof(double)),
    .weights = anno_alloc(WEIGHTS, sizeof(double))
  };
  worker.train = &obj.train;
  feature_scales(windows, window_count, obj.scale);

  double *w = anno_alloc(WEIGHTS, sizeof(double));
  double *grad_obj = anno_alloc(WEIGHTS, sizeof(double));
  double *grad_unreg = anno_alloc(WEIGHTS, sizeof(double));
  for (int f = 0; f < WEIGHTS; f++) w[f] = rng_uniform(-0.25, 0.25);

  /* Unregularized sum over windows via contig_gradient */
  Work wk = {0};
  double ll_sum = 0.0;
  ll_sum += contig_gradient(w, &contigs[0], &wk, grad_unreg);
  ll_sum += contig_gradient(w, &contigs[1], &wk, grad_unreg);
  work_free(&wk);

  /* Regularized objective via crf_objective */
  double loss = crf_objective(w, grad_obj, &obj);
  TEST_ASSERT_INT_EQ(15, obj.positions);
  TEST_ASSERT_DOUBLE_WITHIN(1e-12, ll_sum, obj.likelihood);

  double norm_sq = 0.0;
  for (int f = 0; f < WEIGHTS; f++) norm_sq += w[f] * w[f];
  /* Check minimization sign and L2 penalty: loss = -ll_sum + CRF_L2_C2 * ||w||^2 */
  TEST_ASSERT(ll_sum <= 0.0);
  TEST_ASSERT_DOUBLE_WITHIN(1e-11, -ll_sum + CRF_L2_C2 * norm_sq, loss);

  for (int f = 0; f < WEIGHTS; f++) {
    double expected_g = -grad_unreg[f] + 2.0 * CRF_L2_C2 * w[f];
    TEST_ASSERT_DOUBLE_WITHIN(1e-11, expected_g, grad_obj[f]);
    TEST_ASSERT(obj.scale[f] > 0.0 && obj.scale[f] <= 1.0);
  }

  /* Check preconditioned crf_scaled_objective: w = scale * v, grad_v = scale * grad_w */
  double *v = anno_alloc(WEIGHTS, sizeof(double));
  double *grad_v = anno_alloc(WEIGHTS, sizeof(double));
  for (int f = 0; f < WEIGHTS; f++) v[f] = w[f] / obj.scale[f];
  double scaled_loss = crf_scaled_objective(v, grad_v, &obj);
  TEST_ASSERT_DOUBLE_WITHIN(1e-11, loss, scaled_loss);
  for (int f = 0; f < WEIGHTS; f++)
    TEST_ASSERT_DOUBLE_WITHIN(1e-11, grad_obj[f] * obj.scale[f], grad_v[f]);

  work_free(&worker.work);
  free(worker.gradient);
  free(obj.scale);
  free(obj.weights);
  free(windows);
  free(w);
  free(grad_obj);
  free(grad_unreg);
  free(v);
  free(grad_v);
  free_test_contig(&contigs[0]);
  free_test_contig(&contigs[1]);
}

/* =========================================================================
 * 2. CWT Reference Implementation & Module Tests
 * ========================================================================= */

/* Slow direct time-domain discrete convolution reference for CWT (no FFT). */
static void cwt_direct_reference(const Wavelets *wavelets, const char *seq, int length,
                                 double *features, float *complex_values, float *powers) {
  for (int pos = 0; pos < length; pos++) {
    for (int sc = 0; sc < WAVE_COUNT; sc++) {
      int width = wavelets->widths[sc];
      int radius = width / 2;
      double complex sum = 0.0;
      for (int tap = 0; tap < width; tap++) {
        int disp = tap - radius;
        int q = pos + disp;
        if (q >= 0 && q < length)
          sum += cwt_base_signal(seq[q]) * wavelets->kernel[sc][tap];
      }
      double re = creal(sum), im = cimag(sum);
      size_t feat_idx = (size_t)pos * CWT_CHANNELS + 2 * (size_t)sc;
      if (features) {
        features[feat_idx] = (double)(float)re;
        features[feat_idx + 1] = (double)(float)im;
      }
      if (complex_values) {
        complex_values[feat_idx] = (float)re;
        complex_values[feat_idx + 1] = (float)im;
      }
      if (powers) {
        powers[(size_t)pos * WAVE_COUNT + (size_t)sc] = (float)(re * re + im * im);
      }
    }
  }
}

static void test_cwt_zero_signal(void) {
  Wavelets wavelets;
  wavelets_init(&wavelets);

  const int len = 130;
  char *seq = anno_alloc((size_t)len + 1, 1);
  memset(seq, 'N', (size_t)len);

  double *features = anno_alloc((size_t)len * CWT_CHANNELS, sizeof(double));
  float *cvals = anno_alloc((size_t)len * CWT_CHANNELS, sizeof(float));
  float *powers = anno_alloc((size_t)len * WAVE_COUNT, sizeof(float));

  /* Fill with non-zero sentinel values to detect stale output buffer bugs */
  for (int i = 0; i < len * CWT_CHANNELS; i++) {
    features[i] = 999.0;
    cvals[i] = 999.0f;
  }
  for (int i = 0; i < len * WAVE_COUNT; i++) powers[i] = 999.0f;

  cwt_transform(&wavelets, seq, len, 1, features, cvals, powers);

  for (int i = 0; i < len * CWT_CHANNELS; i++) {
    TEST_ASSERT_DOUBLE_WITHIN(1e-15, 0.0, features[i]);
    TEST_ASSERT_DOUBLE_WITHIN(1e-15, 0.0, cvals[i]);
  }
  for (int i = 0; i < len * WAVE_COUNT; i++)
    TEST_ASSERT_DOUBLE_WITHIN(1e-15, 0.0, powers[i]);

  free(seq);
  free(features);
  free(cvals);
  free(powers);
  wavelets_free(&wavelets);
}

static void test_cwt_linearity(void) {
  Wavelets wavelets;
  wavelets_init(&wavelets);

  const int len = 96;
  char seq_x[97], seq_y[97], seq_xy[97], seq_neg[97], seq_rot[97];
  const char bases[4] = {'A', 'C', 'G', 'T'};
  for (int i = 0; i < len; i++) {
    char b = bases[(i * 3 + 1) & 3];
    /* Disjoint supports: x has bases in first half, y has bases in second half */
    seq_x[i] = (i < len / 2) ? b : 'N';
    seq_y[i] = (i >= len / 2) ? b : 'N';
    seq_xy[i] = b;
    /* In cwt_base_signal: A=(1+i)/sqrt(2), C=(-1+i)/sqrt(2), G=(-1-i)/sqrt(2), T=(1-i)/sqrt(2).
       Swapping A<->G and C<->T multiplies signal by -1.
       Rotating A->C->G->T->A multiplies signal by +i. */
    seq_neg[i] = (b == 'A') ? 'G' : (b == 'G') ? 'A' : (b == 'C') ? 'T' : 'C';
    seq_rot[i] = (b == 'A') ? 'C' : (b == 'C') ? 'G' : (b == 'G') ? 'T' : 'A';
  }
  seq_x[len] = seq_y[len] = seq_xy[len] = seq_neg[len] = seq_rot[len] = '\0';

  double *fx = anno_alloc((size_t)len * CWT_CHANNELS, sizeof(double));
  double *fy = anno_alloc((size_t)len * CWT_CHANNELS, sizeof(double));
  double *fxy = anno_alloc((size_t)len * CWT_CHANNELS, sizeof(double));
  double *fneg = anno_alloc((size_t)len * CWT_CHANNELS, sizeof(double));
  double *frot = anno_alloc((size_t)len * CWT_CHANNELS, sizeof(double));
  float *pxy = anno_alloc((size_t)len * WAVE_COUNT, sizeof(float));
  float *pneg = anno_alloc((size_t)len * WAVE_COUNT, sizeof(float));
  float *prot = anno_alloc((size_t)len * WAVE_COUNT, sizeof(float));

  cwt_transform(&wavelets, seq_x, len, 1, fx, NULL, NULL);
  cwt_transform(&wavelets, seq_y, len, 1, fy, NULL, NULL);
  cwt_transform(&wavelets, seq_xy, len, 1, fxy, NULL, pxy);
  cwt_transform(&wavelets, seq_neg, len, 1, fneg, NULL, pneg);
  cwt_transform(&wavelets, seq_rot, len, 1, frot, NULL, prot);

  for (int i = 0; i < len * CWT_CHANNELS; i++) {
    /* Additivity: W(x + y) = W(x) + W(y) */
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, fx[i] + fy[i], fxy[i]);
    /* Homogeneity by -1: W(-z) = -W(z) */
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, -fxy[i], fneg[i]);
  }
  /* Homogeneity by +i: W(i * z) = i * W(z) => Re' = -Im, Im' = Re */
  for (int pos = 0; pos < len; pos++) {
    for (int sc = 0; sc < WAVE_COUNT; sc++) {
      size_t idx = (size_t)pos * CWT_CHANNELS + 2 * (size_t)sc;
      TEST_ASSERT_DOUBLE_WITHIN(1e-6, -fxy[idx + 1], frot[idx]);
      TEST_ASSERT_DOUBLE_WITHIN(1e-6, fxy[idx], frot[idx + 1]);
      TEST_ASSERT_DOUBLE_WITHIN(1e-6, pxy[pos * WAVE_COUNT + sc], pneg[pos * WAVE_COUNT + sc]);
      TEST_ASSERT_DOUBLE_WITHIN(1e-6, pxy[pos * WAVE_COUNT + sc], prot[pos * WAVE_COUNT + sc]);
    }
  }

  free(fx);
  free(fy);
  free(fxy);
  free(fneg);
  free(frot);
  free(pxy);
  free(pneg);
  free(prot);
  wavelets_free(&wavelets);
}

static void test_cwt_center_impulse(void) {
  Wavelets wavelets;
  wavelets_init(&wavelets);

  const int len = 201;
  const int center = 100;
  char *seq = anno_alloc((size_t)len + 1, 1);
  memset(seq, 'N', (size_t)len);
  seq[center] = 'A';

  double *features = anno_alloc((size_t)len * CWT_CHANNELS, sizeof(double));
  float *powers = anno_alloc((size_t)len * WAVE_COUNT, sizeof(float));
  cwt_transform(&wavelets, seq, len, 1, features, NULL, powers);

  double complex base_a = (1.0 + I) / sqrt(2.0);
  for (int sc = 0; sc < WAVE_COUNT; sc++) {
    double dilation = (double)wave_sizes[sc];
    int width = wavelets.widths[sc];
    int radius = width / 2;
    double complex kernel_sum = 0.0;
    double center_power = powers[center * WAVE_COUNT + sc];
    TEST_ASSERT(center_power > 0.0);

    for (int pos = 0; pos < len; pos++) {
      int d = center - pos; /* tap = radius + (center - pos) */
      size_t idx = (size_t)pos * CWT_CHANNELS + 2 * (size_t)sc;
      double actual_re = features[idx];
      double actual_im = features[idx + 1];

      if (abs(d) > radius) {
        TEST_ASSERT_DOUBLE_WITHIN(1e-13, 0.0, actual_re);
        TEST_ASSERT_DOUBLE_WITHIN(1e-13, 0.0, actual_im);
      } else {
        double complex k_val = sqrt(dilation) * cwt_morlet_integral(
          ((double)d - 0.5) / dilation, ((double)d + 0.5) / dilation);
        kernel_sum += k_val;
        double complex expected = base_a * k_val;
        TEST_ASSERT_DOUBLE_WITHIN(1e-6, creal(expected), actual_re);
        TEST_ASSERT_DOUBLE_WITHIN(1e-6, cimag(expected), actual_im);

        /* Check conjugate symmetry of Morlet kernel: k(-d) == conj(k(d)) */
        double complex k_neg = wavelets.kernel[sc][radius - d];
        TEST_ASSERT_DOUBLE_WITHIN(1e-12, creal(k_val), creal(k_neg));
        TEST_ASSERT_DOUBLE_WITHIN(1e-12, -cimag(k_val), cimag(k_neg));

        /* Impulse response power must peak at the impulse center (d == 0) */
        if (d != 0)
          TEST_ASSERT(powers[pos * WAVE_COUNT + sc] <= (float)(center_power + 1e-7));
      }
    }
    /* Admissibility: Morlet wavelet has zero mean */
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, 0.0, cabs(kernel_sum));
  }

  free(seq);
  free(features);
  free(powers);
  wavelets_free(&wavelets);
}

static void test_cwt_boundary_impulses(void) {
  Wavelets wavelets;
  wavelets_init(&wavelets);

  const int len = 160;
  char *seq_start = anno_alloc((size_t)len + 1, 1);
  char *seq_end = anno_alloc((size_t)len + 1, 1);
  memset(seq_start, 'N', (size_t)len);
  memset(seq_end, 'N', (size_t)len);
  seq_start[0] = 'C';
  seq_end[len - 1] = 'G';

  double *f_start = anno_alloc((size_t)len * CWT_CHANNELS, sizeof(double));
  double *f_end = anno_alloc((size_t)len * CWT_CHANNELS, sizeof(double));
  double *ref_start = anno_alloc((size_t)len * CWT_CHANNELS, sizeof(double));
  double *ref_end = anno_alloc((size_t)len * CWT_CHANNELS, sizeof(double));

  cwt_transform(&wavelets, seq_start, len, 1, f_start, NULL, NULL);
  cwt_transform(&wavelets, seq_end, len, 1, f_end, NULL, NULL);
  cwt_direct_reference(&wavelets, seq_start, len, ref_start, NULL, NULL);
  cwt_direct_reference(&wavelets, seq_end, len, ref_end, NULL, NULL);

  for (int i = 0; i < len * CWT_CHANNELS; i++) {
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, ref_start[i], f_start[i]);
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, ref_end[i], f_end[i]);
  }

  /* Verify no circular convolution wrap-around: impulse at 0 must be 0 at len-1 and vice versa */
  for (int sc = 0; sc < WAVE_COUNT; sc++) {
    size_t last_idx = (size_t)(len - 1) * CWT_CHANNELS + 2 * (size_t)sc;
    TEST_ASSERT_DOUBLE_WITHIN(1e-13, 0.0, f_start[last_idx]);
    TEST_ASSERT_DOUBLE_WITHIN(1e-13, 0.0, f_start[last_idx + 1]);
    size_t first_idx = 2 * (size_t)sc;
    TEST_ASSERT_DOUBLE_WITHIN(1e-13, 0.0, f_end[first_idx]);
    TEST_ASSERT_DOUBLE_WITHIN(1e-13, 0.0, f_end[first_idx + 1]);
  }

  free(seq_start);
  free(seq_end);
  free(f_start);
  free(f_end);
  free(ref_start);
  free(ref_end);
  wavelets_free(&wavelets);
}

static void test_cwt_direct_vs_fft_and_window_boundaries(void) {
  Wavelets wavelets;
  wavelets_init(&wavelets);
  rng_seed(404);

  /* Test lengths: 1, 2, 3, odd/even, and CWT_WINDOW_SIZE (B=1024) boundaries B-1, B, B+1, 2B+1 */
  const int lengths[] = {
    1, 2, 3, 7, 8, 15, 16, 63, 64, 65,
    CWT_WINDOW_SIZE - 1, CWT_WINDOW_SIZE, CWT_WINDOW_SIZE + 1,
    2 * CWT_WINDOW_SIZE + 1
  };
  const char alphabet[] = "ACGTN";
  int num_lengths = (int)(sizeof(lengths) / sizeof(lengths[0]));

  for (int idx = 0; idx < num_lengths; idx++) {
    int L = lengths[idx];
    char *seq = anno_alloc((size_t)L + 1, 1);
    for (int i = 0; i < L; i++) {
      int r = (int)(rng_uniform(0.0, 100.0));
      seq[i] = (r < 92) ? alphabet[r & 3] : 'N';
    }
    seq[L] = '\0';

    double *fft_feat = anno_alloc((size_t)L * CWT_CHANNELS, sizeof(double));
    float *fft_cval = anno_alloc((size_t)L * CWT_CHANNELS, sizeof(float));
    float *fft_pow = anno_alloc((size_t)L * WAVE_COUNT, sizeof(float));

    double *ref_feat = anno_alloc((size_t)L * CWT_CHANNELS, sizeof(double));
    float *ref_cval = anno_alloc((size_t)L * CWT_CHANNELS, sizeof(float));
    float *ref_pow = anno_alloc((size_t)L * WAVE_COUNT, sizeof(float));

    cwt_transform(&wavelets, seq, L, 1, fft_feat, fft_cval, fft_pow);
    cwt_direct_reference(&wavelets, seq, L, ref_feat, ref_cval, ref_pow);

    for (int i = 0; i < L * CWT_CHANNELS; i++) {
      TEST_ASSERT_DOUBLE_WITHIN(1e-5, ref_feat[i], fft_feat[i]);
      TEST_ASSERT_DOUBLE_WITHIN(1e-5, ref_cval[i], fft_cval[i]);
    }
    for (int i = 0; i < L * WAVE_COUNT; i++)
      TEST_ASSERT_DOUBLE_WITHIN(1e-5, ref_pow[i], fft_pow[i]);

    free(seq);
    free(fft_feat);
    free(fft_cval);
    free(fft_pow);
    free(ref_feat);
    free(ref_cval);
    free(ref_pow);
  }

  /* Also verify rfft.h radix-2 and Bluestein FFT against direct O(N^2) DFT */
  const size_t fft_sizes[] = {8, 16, 17, 23};
  for (size_t s = 0; s < sizeof(fft_sizes) / sizeof(fft_sizes[0]); s++) {
    size_t N = fft_sizes[s];
    double complex *orig = anno_alloc(N, sizeof(double complex));
    double complex *vec = anno_alloc(N, sizeof(double complex));
    for (size_t i = 0; i < N; i++)
      orig[i] = vec[i] = rng_uniform(-1.0, 1.0) + I * rng_uniform(-1.0, 1.0);

    fft_transform(vec, N, false);
    for (size_t k = 0; k < N; k++) {
      double complex dft_k = 0.0;
      for (size_t n = 0; n < N; n++) {
        double angle = -2.0 * acos(-1.0) * (double)(k * n) / (double)N;
        dft_k += orig[n] * (cos(angle) + I * sin(angle));
      }
      TEST_ASSERT_DOUBLE_WITHIN(1e-11, creal(dft_k), creal(vec[k]));
      TEST_ASSERT_DOUBLE_WITHIN(1e-11, cimag(dft_k), cimag(vec[k]));
    }

    fft_transform(vec, N, true);
    for (size_t i = 0; i < N; i++) {
      TEST_ASSERT_DOUBLE_WITHIN(1e-11, creal(orig[i]), creal(vec[i]) / (double)N);
      TEST_ASSERT_DOUBLE_WITHIN(1e-11, cimag(orig[i]), cimag(vec[i]) / (double)N);
    }
    free(orig);
    free(vec);
  }

  wavelets_free(&wavelets);
}

static void test_cwt_sinusoid_scale_response(void) {
  Wavelets wavelets;
  wavelets_init(&wavelets);

  /* In cwt_base_signal, A -> C -> G -> T has phase step +pi/2 per base (period 4, omega = pi/2).
     A -> C -> C -> G -> T -> T has period 6 (omega ~= pi/3).
     For Morlet wavelet with omega_0 = 6.0, effective angular frequency at scale s is ~ 6.0 / s:
     - scale 4 has peak frequency 6/4 = 1.50 rad/sample (closest to period 4: pi/2 ~= 1.571)
     - scale 6 has peak frequency 6/6 = 1.00 rad/sample (closer to longer periods: period 6, pi/3 ~= 1.047)
     Test period-4 and period-6 sequences in the matching positive-frequency direction. */
  const int len = 240;
  char seq_p4[241], seq_p6[241];
  const char pat4[4] = {'A', 'C', 'G', 'T'};
  const char pat6[6] = {'A', 'C', 'C', 'G', 'T', 'T'};
  for (int i = 0; i < len; i++) {
    seq_p4[i] = pat4[i % 4];
    seq_p6[i] = pat6[i % 6];
  }
  seq_p4[len] = seq_p6[len] = '\0';

  float *pow_p4 = anno_alloc((size_t)len * WAVE_COUNT, sizeof(float));
  float *pow_p6 = anno_alloc((size_t)len * WAVE_COUNT, sizeof(float));
  cwt_transform(&wavelets, seq_p4, len, 1, NULL, NULL, pow_p4);
  cwt_transform(&wavelets, seq_p6, len, 1, NULL, NULL, pow_p6);

  /* Measure average interior power (avoiding boundary cone of influence [60, 180)) */
  double mean_p4[WAVE_COUNT] = {0}, mean_p6[WAVE_COUNT] = {0};
  for (int pos = 60; pos < 180; pos++) {
    for (int sc = 0; sc < WAVE_COUNT; sc++) {
      mean_p4[sc] += pow_p4[pos * WAVE_COUNT + sc];
      mean_p6[sc] += pow_p6[pos * WAVE_COUNT + sc];
    }
  }

  /* For period-4 sinusoid (omega = pi/2 ~= 1.57), scale=4 (index 0) must dominate scale=5 and scale=6 */
  TEST_ASSERT(mean_p4[0] > mean_p4[1]);
  TEST_ASSERT(mean_p4[1] > mean_p4[2]);

  /* For period-6 signal (omega = pi/3 ~= 1.05), scale=6 (index 2) must dominate scale=5 and scale=4 */
  TEST_ASSERT(mean_p6[2] > mean_p6[1]);
  TEST_ASSERT(mean_p6[1] > mean_p6[0]);

  free(pow_p4);
  free(pow_p6);
  wavelets_free(&wavelets);
}

static void test_cwt_repeat_calls_and_threads(void) {
  Wavelets wavelets;
  wavelets_init(&wavelets);

  const int len_a = 2300;
  const int len_b = 1500;
  char *seq_a = anno_alloc((size_t)len_a + 1, 1);
  char *seq_b = anno_alloc((size_t)len_b + 1, 1);
  const char bases[4] = {'A', 'C', 'G', 'T'};
  for (int i = 0; i < len_a; i++) seq_a[i] = bases[(i * 7 + 3) & 3];
  for (int i = 0; i < len_b; i++) seq_b[i] = bases[(i * 5 + 1) & 3];

  double *feat_a1 = anno_alloc((size_t)len_a * CWT_CHANNELS, sizeof(double));
  float *cval_a1 = anno_alloc((size_t)len_a * CWT_CHANNELS, sizeof(float));
  float *pow_a1 = anno_alloc((size_t)len_a * WAVE_COUNT, sizeof(float));

  double *feat_a2 = anno_alloc((size_t)len_a * CWT_CHANNELS, sizeof(double));
  float *cval_a2 = anno_alloc((size_t)len_a * CWT_CHANNELS, sizeof(float));
  float *pow_a2 = anno_alloc((size_t)len_a * WAVE_COUNT, sizeof(float));

  /* Call A -> B (overwriting buffers) -> A again */
  cwt_transform(&wavelets, seq_a, len_a, 1, feat_a1, cval_a1, pow_a1);
  cwt_transform(&wavelets, seq_b, len_b, 1, feat_a2, cval_a2, pow_a2);
  cwt_transform(&wavelets, seq_a, len_a, 1, feat_a2, cval_a2, pow_a2);

  TEST_ASSERT_INT_EQ(0, memcmp(feat_a1, feat_a2, (size_t)len_a * CWT_CHANNELS * sizeof(double)));
  TEST_ASSERT_INT_EQ(0, memcmp(cval_a1, cval_a2, (size_t)len_a * CWT_CHANNELS * sizeof(float)));
  TEST_ASSERT_INT_EQ(0, memcmp(pow_a1, pow_a2, (size_t)len_a * WAVE_COUNT * sizeof(float)));

  /* Multi-threaded call (4 threads across 3 windows) must reproduce 1-thread output */
  memset(feat_a2, 0, (size_t)len_a * CWT_CHANNELS * sizeof(double));
  memset(cval_a2, 0, (size_t)len_a * CWT_CHANNELS * sizeof(float));
  memset(pow_a2, 0, (size_t)len_a * WAVE_COUNT * sizeof(float));
  cwt_transform(&wavelets, seq_a, len_a, 4, feat_a2, cval_a2, pow_a2);

  for (int i = 0; i < len_a * CWT_CHANNELS; i++) {
    TEST_ASSERT_DOUBLE_WITHIN(1e-14, feat_a1[i], feat_a2[i]);
    TEST_ASSERT_DOUBLE_WITHIN(1e-7, cval_a1[i], cval_a2[i]);
  }
  for (int i = 0; i < len_a * WAVE_COUNT; i++)
    TEST_ASSERT_DOUBLE_WITHIN(1e-7, pow_a1[i], pow_a2[i]);

  /* Also verify mmap output helper (cwt_map_output / cwt_unmap_output) */
  FILE *tmp_cval = tmpfile();
  FILE *tmp_pow = tmpfile();
  TEST_ASSERT(tmp_cval != NULL);
  TEST_ASSERT(tmp_pow != NULL);
  void *cval_map = NULL, *pow_map = NULL;
  size_t cval_map_sz = 0, pow_map_sz = 0;
  off_t cval_pos = 0, pow_pos = 0;
  size_t cval_cnt = (size_t)len_a * CWT_CHANNELS;
  size_t pow_cnt = (size_t)len_a * WAVE_COUNT;
  float *mapped_cval = cwt_map_output(tmp_cval, cval_cnt, &cval_map, &cval_map_sz, &cval_pos);
  float *mapped_pow = cwt_map_output(tmp_pow, pow_cnt, &pow_map, &pow_map_sz, &pow_pos);
  cwt_transform(&wavelets, seq_a, len_a, 2, NULL, mapped_cval, mapped_pow);
  TEST_ASSERT_INT_EQ(0, memcmp(cval_a1, mapped_cval, cval_cnt * sizeof(float)));
  TEST_ASSERT_INT_EQ(0, memcmp(pow_a1, mapped_pow, pow_cnt * sizeof(float)));
  cwt_unmap_output(tmp_cval, mapped_cval, cval_cnt, cval_map, cval_map_sz, cval_pos);
  cwt_unmap_output(tmp_pow, mapped_pow, pow_cnt, pow_map, pow_map_sz, pow_pos);
  fclose(tmp_cval);
  fclose(tmp_pow);

  free(seq_a);
  free(seq_b);
  free(feat_a1);
  free(cval_a1);
  free(pow_a1);
  free(feat_a2);
  free(cval_a2);
  free(pow_a2);
  wavelets_free(&wavelets);
}

/* =========================================================================
 * 3. L-BFGS Module Tests
 * ========================================================================= */

static void test_lbfgs_crf_directional_derivatives(void) {
  crf_init();
  rng_seed(20261020);

  const char *seq1 =
    "ACGTACGTACGTACGTACGTACGTACGTAC"
    "ATGAAATAA"
    "ACGTACGTACGTACGTACGTACGTACGTACGT";
  const uint8_t gold1[9] = {1, 2, 3, 4, 4, 4, 1, 2, 3};
  const char *seq2 =
    "ACGTACGTACGTACGTACGTACGTACGTAC"
    "TTATCACAT"
    "ACGTACGTACGTACGTACGTACGTACGTACGT";
  const uint8_t gold2[9] = {9, 8, 7, 12, 12, 12, 9, 8, 7};

  Contig contigs[2];
  contigs[0] = make_test_contig(seq1, 30, 9, gold1);
  contigs[1] = make_test_contig(seq2, 30, 9, gold2);

  CrfWindow *windows = NULL;
  size_t skipped = 0;
  size_t win_count = make_windows(contigs, 2, &windows, &skipped);
  TEST_ASSERT_INT_EQ(2, win_count);

  CrfWorker worker = {
    .train = NULL,
    .thread_id = 0,
    .gradient = anno_alloc(WEIGHTS, sizeof(double))
  };
  CrfObjective obj = {
    .train = {.windows = windows, .count = win_count, .skipped = 0, .threads = 1},
    .workers = &worker,
    .threads = NULL,
    .scale = anno_alloc(WEIGHTS, sizeof(double)),
    .weights = anno_alloc(WEIGHTS, sizeof(double))
  };
  worker.train = &obj.train;
  feature_scales(windows, win_count, obj.scale);

  double *theta = anno_alloc(WEIGHTS, sizeof(double));
  double *theta_p = anno_alloc(WEIGHTS, sizeof(double));
  double *theta_m = anno_alloc(WEIGHTS, sizeof(double));
  double *grad = anno_alloc(WEIGHTS, sizeof(double));
  double *grad_repeat = anno_alloc(WEIGHTS, sizeof(double));
  double *grad_dummy = anno_alloc(WEIGHTS, sizeof(double));
  double *dir = anno_alloc(WEIGHTS, sizeof(double));

  const double point_scales[3] = {0.0, 0.15, 0.8};
  const double steps[] = {1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7};
  const int step_count = (int)(sizeof(steps) / sizeof(steps[0]));

  for (int pt = 0; pt < 3; pt++) {
    for (int f = 0; f < WEIGHTS; f++)
      theta[f] = point_scales[pt] * rng_uniform(-1.0, 1.0);

    /* Repeatability check: evaluating the same point twice yields identical loss & gradient */
    double loss1 = crf_objective(theta, grad, &obj);
    double loss2 = crf_objective(theta, grad_repeat, &obj);
    TEST_ASSERT_DOUBLE_WITHIN(1e-15, loss1, loss2);
    for (int f = 0; f < WEIGHTS; f++)
      TEST_ASSERT_DOUBLE_WITHIN(1e-15, grad[f], grad_repeat[f]);

    /* Also compute unregularized NLL gradient for single window and multi-window */
    Work wk = {0};
    double grad_unreg[WEIGHTS] = {0};
    double ll_single = contig_gradient(theta, &contigs[0], &wk, grad_unreg);
    (void)ll_single;

    /* Test 4 random unit directions per parameter point */
    for (int d = 0; d < 4; d++) {
      double norm = 0.0;
      for (int f = 0; f < WEIGHTS; f++) {
        dir[f] = rng_uniform(-1.0, 1.0);
        norm = hypot(norm, dir[f]);
      }
      for (int f = 0; f < WEIGHTS; f++) dir[f] /= norm;

      double directional_unreg = 0.0;
      double directional_reg = 0.0;
      for (int f = 0; f < WEIGHTS; f++) {
        directional_unreg += (-grad_unreg[f]) * dir[f];
        directional_reg += grad[f] * dir[f];
      }

      /* 1. Check unregularized NLL directional derivative across step sizes h */
      double best_unreg_err = INFINITY;
      double best_reg_err = INFINITY;
      for (int s = 0; s < step_count; s++) {
        double h = steps[s];
        for (int f = 0; f < WEIGHTS; f++) {
          theta_p[f] = theta[f] + h * dir[f];
          theta_m[f] = theta[f] - h * dir[f];
        }
        double nll_p = -contig_gradient(theta_p, &contigs[0], &wk, grad_dummy);
        double nll_m = -contig_gradient(theta_m, &contigs[0], &wk, grad_dummy);
        double fd_unreg = (nll_p - nll_m) / (2.0 * h);
        double err_unreg = fabs(fd_unreg - directional_unreg);
        if (err_unreg < best_unreg_err) best_unreg_err = err_unreg;

        /* 2. Check regularized multi-window crf_objective directional derivative */
        double obj_p = crf_objective(theta_p, grad_dummy, &obj);
        double obj_m = crf_objective(theta_m, grad_dummy, &obj);
        double fd_reg = (obj_p - obj_m) / (2.0 * h);
        double err_reg = fabs(fd_reg - directional_reg);
        if (err_reg < best_reg_err) best_reg_err = err_reg;
      }
      TEST_ASSERT_DOUBLE_WITHIN(1e-8, 0.0, best_unreg_err);
      TEST_ASSERT_DOUBLE_WITHIN(1e-7, 0.0, best_reg_err);
    }

    /* 3. Check scaled objective directional derivative across step sizes */
    double loss_scaled = crf_scaled_objective(theta, grad_repeat, &obj);
    (void)loss_scaled;
    double dir_scaled = 0.0;
    for (int f = 0; f < WEIGHTS; f++)
      dir_scaled += grad_repeat[f] * dir[f];
    double best_scaled_err = INFINITY;
    for (int hi = 0; hi < step_count; hi++) {
      double h = steps[hi];
      for (int f = 0; f < WEIGHTS; f++) {
        theta_p[f] = theta[f] + h * dir[f];
        theta_m[f] = theta[f] - h * dir[f];
      }
      double sp = crf_scaled_objective(theta_p, grad_dummy, &obj);
      double sm = crf_scaled_objective(theta_m, grad_dummy, &obj);
      double err_scaled = fabs((sp - sm) / (2.0 * h) - dir_scaled);
      if (err_scaled < best_scaled_err) best_scaled_err = err_scaled;
    }
    TEST_ASSERT_DOUBLE_WITHIN(1e-8, 0.0, best_scaled_err);

    work_free(&wk);
  }

  work_free(&worker.work);
  free(worker.gradient);
  free(obj.scale);
  free(obj.weights);
  free(windows);
  free(theta);
  free(theta_p);
  free(theta_m);
  free(grad);
  free(grad_repeat);
  free(grad_dummy);
  free(dir);
  free_test_contig(&contigs[0]);
  free_test_contig(&contigs[1]);
}

typedef struct {
  size_t n;
  const double *center;
  const double *diag;
  const double *matrix;
  int eval_count;
  double last_progress_obj;
  double last_progress_gnorm;
  int last_progress_iter;
} QuadContext;

static double eval_isotropic_quad(const double *x, double *g, void *ctx) {
  QuadContext *q = ctx;
  q->eval_count++;
  double f = 0.0;
  for (size_t i = 0; i < q->n; i++) {
    double d = x[i] - q->center[i];
    f += 0.5 * d * d;
    g[i] = d;
  }
  return f;
}

static double eval_diag_quad(const double *x, double *g, void *ctx) {
  QuadContext *q = ctx;
  q->eval_count++;
  double f = 0.0;
  for (size_t i = 0; i < q->n; i++) {
    double d = x[i] - q->center[i];
    f += 0.5 * q->diag[i] * d * d;
    g[i] = q->diag[i] * d;
  }
  return f;
}

static double eval_coupled_quad(const double *x, double *g, void *ctx) {
  QuadContext *q = ctx;
  q->eval_count++;
  size_t n = q->n;
  double f = 0.0;
  for (size_t i = 0; i < n; i++) {
    double gi = 0.0;
    for (size_t j = 0; j < n; j++)
      gi += q->matrix[i * n + j] * (x[j] - q->center[j]);
    g[i] = gi;
    f += 0.5 * (x[i] - q->center[i]) * gi;
  }
  return f;
}

static void record_progress(void *ctx, int iteration, double objective, double step,
                            int trials, size_t n, const double *gradient) {
  (void)step;
  (void)trials;
  QuadContext *q = ctx;
  q->last_progress_iter = iteration;
  q->last_progress_obj = objective;
  q->last_progress_gnorm = lbfgs_norm(n, gradient);
}

static void test_lbfgs_quadratic_isotropic_diagonal_coupled(void) {
  /* 1. Basic isotropic quadratic: 1/2 ||x - c||^2 */
  const double c4[4] = {1.5, -2.25, 0.75, 4.0};
  double x[4] = {-3.0, 5.0, -1.0, 0.0};
  double g[4] = {0};
  QuadContext q1 = {.n = 4, .center = c4};
  int completed = 0;
  LbfgsStatus st = lbfgs_minimize(4, x, eval_isotropic_quad, record_progress, &q1, 50, &completed);
  TEST_ASSERT_INT_EQ(LBFGS_CONVERGED, st);
  double f_final = eval_isotropic_quad(x, g, &q1);
  TEST_ASSERT_DOUBLE_WITHIN(1e-12, q1.last_progress_obj, f_final);
  TEST_ASSERT_DOUBLE_WITHIN(1e-12, q1.last_progress_gnorm, lbfgs_norm(4, g));
  for (int i = 0; i < 4; i++) TEST_ASSERT_DOUBLE_WITHIN(1e-5, c4[i], x[i]);

  /* 2. Diagonal quadratic with widely varying curvatures d_i in [0.1, 100] */
  const double d6[6] = {0.1, 0.5, 2.0, 10.0, 40.0, 100.0};
  const double c6[6] = {1.0, -2.0, 3.0, -0.5, 0.25, -1.5};
  double x6[6] = {5.0, 4.0, -3.0, 2.0, -4.0, 3.0};
  double g6[6] = {0};
  QuadContext q2 = {.n = 6, .center = c6, .diag = d6};
  st = lbfgs_minimize(6, x6, eval_diag_quad, record_progress, &q2, 100, &completed);
  TEST_ASSERT_INT_EQ(LBFGS_CONVERGED, st);
  f_final = eval_diag_quad(x6, g6, &q2);
  TEST_ASSERT_DOUBLE_WITHIN(1e-12, q2.last_progress_obj, f_final);
  TEST_ASSERT_DOUBLE_WITHIN(1e-12, q2.last_progress_gnorm, lbfgs_norm(6, g6));
  for (int i = 0; i < 6; i++) TEST_ASSERT_DOUBLE_WITHIN(1e-4, c6[i], x6[i]);

  /* 3. Coupled non-diagonal SPD quadratic: 1/2 (x - c)^T A (x - c), A > 0 */
  const double A[16] = {
    6.0, 2.0, 1.0, 0.5,
    2.0, 5.0, 1.5, 0.2,
    1.0, 1.5, 4.0, 1.0,
    0.5, 0.2, 1.0, 3.0
  };
  double x_coupled[4] = {-2.0, 3.0, 4.0, -5.0};
  QuadContext q3 = {.n = 4, .center = c4, .matrix = A};
  st = lbfgs_minimize(4, x_coupled, eval_coupled_quad, record_progress, &q3, 100, &completed);
  TEST_ASSERT_INT_EQ(LBFGS_CONVERGED, st);
  f_final = eval_coupled_quad(x_coupled, g, &q3);
  TEST_ASSERT_DOUBLE_WITHIN(1e-12, q3.last_progress_obj, f_final);
  TEST_ASSERT_DOUBLE_WITHIN(1e-12, q3.last_progress_gnorm, lbfgs_norm(4, g));
  for (int i = 0; i < 4; i++) TEST_ASSERT_DOUBLE_WITHIN(1e-5, c4[i], x_coupled[i]);
}

static double eval_rosenbrock(const double *x, double *g, void *ctx) {
  QuadContext *q = ctx;
  q->eval_count++;
  double t1 = 1.0 - x[0];
  double t2 = x[1] - x[0] * x[0];
  g[0] = -2.0 * t1 - 400.0 * x[0] * t2;
  g[1] = 200.0 * t2;
  return t1 * t1 + 100.0 * t2 * t2;
}

static void test_lbfgs_rosenbrock(void) {
  double x[2] = {-1.2, 1.0};
  double g[2] = {0};
  QuadContext q = {.n = 2};
  int completed = 0;
  LbfgsStatus st = lbfgs_minimize(2, x, eval_rosenbrock, record_progress, &q, 200, &completed);
  TEST_ASSERT_INT_EQ(LBFGS_CONVERGED, st);
  TEST_ASSERT_INT_EQ(completed, q.last_progress_iter);

  double f_final = eval_rosenbrock(x, g, &q);
  TEST_ASSERT_DOUBLE_WITHIN(1e-12, q.last_progress_obj, f_final);
  TEST_ASSERT_DOUBLE_WITHIN(1e-12, q.last_progress_gnorm, lbfgs_norm(2, g));
  TEST_ASSERT(f_final < 1e-10);
  TEST_ASSERT(lbfgs_norm(2, g) < 1e-5);
  TEST_ASSERT_DOUBLE_WITHIN(1e-5, 1.0, x[0]);
  TEST_ASSERT_DOUBLE_WITHIN(1e-5, 1.0, x[1]);
}

/* Helper objectives for termination state tests */
static double eval_flat_offset_quad(const double *x, double *g, void *ctx) {
  (void)ctx;
  /* 12D diagonal quadratic with distinct curvatures + 1e8 constant offset so
     L-BFGS (history m=6 < 12) needs > 10 steps for gnorm, while relative objective
     change over LBFGS_PAST (10) steps is < lbfgs_delta (1e-5). */
  double f = 1e8;
  for (int i = 0; i < 12; i++) {
    double scale = 1.0 + (double)i * 2.5;
    f += 0.5 * scale * x[i] * x[i];
    g[i] = scale * x[i];
  }
  return f;
}

static double eval_bad_linesearch(const double *x, double *g, void *ctx) {
  (void)ctx;
  /* Reports negative slope at origin, but function increases in that direction */
  g[0] = 1.0;
  return (x[0] < 0.0) ? (10.0 - x[0]) : 1.0;
}

static double eval_nonfinite(const double *x, double *g, void *ctx) {
  (void)x;
  (void)ctx;
  g[0] = NAN;
  return NAN;
}

static void test_lbfgs_termination_states(void) {
  const double c2[2] = {2.0, -3.0};
  QuadContext q = {.n = 2, .center = c2};
  int completed = 99;

  /* 1. Initial point is already optimal (grad = 0 at start) */
  double x_opt[2] = {2.0, -3.0};
  LbfgsStatus st = lbfgs_minimize(2, x_opt, eval_isotropic_quad, NULL, &q, 50, &completed);
  TEST_ASSERT_INT_EQ(LBFGS_CONVERGED, st);
  TEST_ASSERT_INT_EQ(0, completed);
  TEST_ASSERT_INT_EQ(1, q.eval_count);

  /* 2. Iteration limit reached (LBFGS_MAX_ITERATIONS) */
  double x_rosen[2] = {-1.2, 1.0};
  st = lbfgs_minimize(2, x_rosen, eval_rosenbrock, NULL, &q, 3, &completed);
  TEST_ASSERT_INT_EQ(LBFGS_MAX_ITERATIONS, st);
  TEST_ASSERT_INT_EQ(3, completed);

  /* 3. Past objective stopping criterion (LBFGS_STOPPING_CRITERION) */
  double x_flat[12];
  for (int i = 0; i < 12; i++) x_flat[i] = (i % 2 == 0 ? 1.0 : -1.0) * (double)(i + 2);
  st = lbfgs_minimize(12, x_flat, eval_flat_offset_quad, NULL, NULL, 100, &completed);
  TEST_ASSERT_INT_EQ(LBFGS_STOPPING_CRITERION, st);
  TEST_ASSERT(completed >= LBFGS_PAST);

  /* 4. Line search failure (LBFGS_LINE_SEARCH_FAILURE) and restoration of previous x */
  double x_ls[1] = {0.0};
  st = lbfgs_minimize(1, x_ls, eval_bad_linesearch, NULL, NULL, 10, &completed);
  TEST_ASSERT_INT_EQ(LBFGS_LINE_SEARCH_FAILURE, st);
  TEST_ASSERT_INT_EQ(0, completed);
  TEST_ASSERT_DOUBLE_WITHIN(1e-15, 0.0, x_ls[0]);

  /* 5. Non-finite objective / gradient (LBFGS_NONFINITE) */
  double x_nan[1] = {0.0};
  st = lbfgs_minimize(1, x_nan, eval_nonfinite, NULL, NULL, 10, &completed);
  TEST_ASSERT_INT_EQ(LBFGS_NONFINITE, st);

  /* 6. Invalid arguments (LBFGS_INVALID_ARGUMENT) */
  TEST_ASSERT_INT_EQ(LBFGS_INVALID_ARGUMENT,
                     lbfgs_minimize(0, x_opt, eval_isotropic_quad, NULL, &q, 10, &completed));
  TEST_ASSERT_INT_EQ(LBFGS_INVALID_ARGUMENT,
                     lbfgs_minimize(2, NULL, eval_isotropic_quad, NULL, &q, 10, &completed));
  TEST_ASSERT_INT_EQ(LBFGS_INVALID_ARGUMENT,
                     lbfgs_minimize(2, x_opt, NULL, NULL, &q, 10, &completed));
  TEST_ASSERT_INT_EQ(LBFGS_INVALID_ARGUMENT,
                     lbfgs_minimize(2, x_opt, eval_isotropic_quad, NULL, &q, 0, &completed));
  TEST_ASSERT_INT_EQ(LBFGS_INVALID_ARGUMENT,
                     lbfgs_minimize(2, x_opt, eval_isotropic_quad, NULL, &q, 10, NULL));
}

/* Compute explicit dense BFGS inverse Hessian direction -H_k * g from kept (s_i, y_i) pairs. */
static void explicit_bfgs_direction(size_t n, const double *g, const double *s,
                                    const double *y, int kept, double *out_dir) {
  if (kept <= 0) {
    for (size_t i = 0; i < n; i++) out_dir[i] = -g[i];
    return;
  }
  const double *s_last = s + (size_t)(kept - 1) * n;
  const double *y_last = y + (size_t)(kept - 1) * n;
  double ys_last = lbfgs_dot(n, s_last, y_last);
  double yy_last = lbfgs_dot(n, y_last, y_last);
  if (!(ys_last > 0.0 && yy_last > 0.0)) {
    for (size_t i = 0; i < n; i++) out_dir[i] = -g[i];
    return;
  }
  double gamma = ys_last / yy_last;

  double H[16] = {0}, next_H[16] = {0}, V[16] = {0}, tmp[16] = {0};
  for (size_t i = 0; i < n; i++) H[i * n + i] = gamma;

  for (int k = 0; k < kept; k++) {
    const double *sk = s + (size_t)k * n;
    const double *yk = y + (size_t)k * n;
    double rho_k = 1.0 / lbfgs_dot(n, sk, yk);
    /* V = I - rho_k * y_k * s_k^T, so H_{k+1} = V^T * H_k * V + rho_k * s_k * s_k^T */
    for (size_t r = 0; r < n; r++)
      for (size_t c = 0; c < n; c++)
        V[r * n + c] = (r == c ? 1.0 : 0.0) - rho_k * yk[r] * sk[c];

    for (size_t r = 0; r < n; r++) {
      for (size_t c = 0; c < n; c++) {
        double sum = 0.0;
        for (size_t m = 0; m < n; m++) sum += V[m * n + r] * H[m * n + c];
        tmp[r * n + c] = sum;
      }
    }
    for (size_t r = 0; r < n; r++) {
      for (size_t c = 0; c < n; c++) {
        double sum = 0.0;
        for (size_t m = 0; m < n; m++) sum += tmp[r * n + m] * V[m * n + c];
        next_H[r * n + c] = sum + rho_k * sk[r] * sk[c];
      }
    }
    memcpy(H, next_H, n * n * sizeof(double));
  }

  for (size_t r = 0; r < n; r++) {
    double sum = 0.0;
    for (size_t c = 0; c < n; c++) sum += H[r * n + c] * g[c];
    out_dir[r] = -sum;
  }
}

static void test_lbfgs_two_loop_history_and_wolfe(void) {
  rng_seed(8080);
  const size_t n = 4;
  double s[LBFGS_HISTORY * 4] = {0};
  double y[LBFGS_HISTORY * 4] = {0};
  double rho[LBFGS_HISTORY] = {0};
  double alpha[LBFGS_HISTORY] = {0};
  double g[4] = {0.7, -1.3, 2.1, -0.4};
  double dir_two_loop[4] = {0};
  double dir_explicit[4] = {0};

  /* 1. Simulate 10 history updates (exceeding LBFGS_HISTORY = 6 so buffer shifts/wraps)
     and compare two-loop recursion against explicit dense BFGS matrix at every step. */
  int kept = 0;
  for (int step = 0; step < 10; step++) {
    double sk[4], yk[4];
    for (size_t j = 0; j < n; j++) {
      sk[j] = rng_uniform(-1.0, 1.0);
      yk[j] = 2.0 * sk[j] + rng_uniform(-0.2, 0.2); /* ensure s^T y > 0 */
    }
    double ys = lbfgs_dot(n, sk, yk);
    TEST_ASSERT(ys > 0.0);

    if (kept == LBFGS_HISTORY) {
      memmove(s, s + n, (LBFGS_HISTORY - 1) * n * sizeof(double));
      memmove(y, y + n, (LBFGS_HISTORY - 1) * n * sizeof(double));
      memmove(rho, rho + 1, (LBFGS_HISTORY - 1) * sizeof(double));
      kept--;
    }
    memcpy(s + (size_t)kept * n, sk, n * sizeof(double));
    memcpy(y + (size_t)kept * n, yk, n * sizeof(double));
    rho[kept++] = 1.0 / ys;

    int ret_kept = lbfgs_direction(n, g, s, y, rho, kept, alpha, dir_two_loop);
    TEST_ASSERT_INT_EQ(kept, ret_kept);
    explicit_bfgs_direction(n, g, s, y, kept, dir_explicit);

    for (size_t j = 0; j < n; j++)
      TEST_ASSERT_DOUBLE_WITHIN(1e-12, dir_explicit[j], dir_two_loop[j]);
  }

  /* 2. Negative curvature policy: if latest s^T y <= 0, lbfgs_direction resets to -g and returns kept = 0 */
  for (size_t j = 0; j < n; j++)
    y[(size_t)(kept - 1) * n + j] = -s[(size_t)(kept - 1) * n + j];
  int ret_kept = lbfgs_direction(n, g, s, y, rho, kept, alpha, dir_two_loop);
  TEST_ASSERT_INT_EQ(0, ret_kept);
  for (size_t j = 0; j < n; j++)
    TEST_ASSERT_DOUBLE_WITHIN(1e-15, -g[j], dir_two_loop[j]);

  /* 3. Verify More-Thuente line search satisfies Armijo and strong Wolfe conditions */
  double base_x[2] = {-1.2, 1.0};
  double x[2] = {-1.2, 1.0};
  double grad[2] = {0};
  QuadContext q = {.n = 2};
  double f0 = eval_rosenbrock(base_x, grad, &q);
  double d[2] = {-grad[0], -grad[1]};
  double dginit = lbfgs_dot(2, grad, d);
  TEST_ASSERT(dginit < 0.0);

  double obj = f0;
  double step = 1.0 / lbfgs_norm(2, grad);
  int trials = 0;
  LbfgsStatus ls = lbfgs_morethuente(2, x, grad, base_x, d, &obj, &step,
                                     eval_rosenbrock, &q, &trials);
  TEST_ASSERT_INT_EQ(LBFGS_CONVERGED, ls);
  TEST_ASSERT(step > 0.0 && trials >= 1);

  double dg_new = lbfgs_dot(2, grad, d);
  /* Armijo sufficient decrease: f(x + step*d) <= f(x) + step * ftol * dginit */
  TEST_ASSERT(obj <= f0 + step * lbfgs_ftol * dginit);
  /* Strong Wolfe curvature: |g(x + step*d)^T d| <= gtol * |dginit| */
  TEST_ASSERT(fabs(dg_new) <= lbfgs_gtol * (-dginit));
}

/* =========================================================================
 * 4. End-to-End Integration Tests (CWT + CRF + L-BFGS)
 * ========================================================================= */

/* Independent reference optimizer (gradient descent with Armijo backtracking line search)
   for cross-checking L-BFGS on small strongly convex CRF objectives. */
static double reference_gd_optimize(size_t n, double *x, lbfgs_evaluate evaluate,
                                    void *context, int max_iters) {
  double *g = anno_alloc(n, sizeof(double));
  double *x_try = anno_alloc(n, sizeof(double));
  double *g_try = anno_alloc(n, sizeof(double));
  double f = evaluate(x, g, context);

  for (int it = 0; it < max_iters; it++) {
    double gnorm_sq = lbfgs_dot(n, g, g);
    if (sqrt(gnorm_sq) < 1e-8) break;
    double step = 0.5;
    while (step > 1e-12) {
      for (size_t i = 0; i < n; i++) x_try[i] = x[i] - step * g[i];
      double f_try = evaluate(x_try, g_try, context);
      if (f_try <= f - 1e-4 * step * gnorm_sq) {
        memcpy(x, x_try, n * sizeof(double));
        memcpy(g, g_try, n * sizeof(double));
        f = f_try;
        break;
      }
      step *= 0.5;
    }
  }
  free(g);
  free(x_try);
  free(g_try);
  return f;
}

static void test_integration_fixed_features_crf_lbfgs_and_reference_solver(void) {
  crf_init();

  /* Step 1: Fixed small CWT features -> verify CRF objective & gradient by brute-force */
  const char *seq =
    "ACGTACGTACGTACGTACGTACGTACGTAC"
    "ATGAAATAA"
    "ACGTACGTACGTACGTACGTACGTACGTACGT";
  const uint8_t gold[9] = {1, 2, 3, 4, 4, 4, 1, 2, 3};
  Contig c = make_test_contig(seq, 30, 9, gold);

  CrfWindow *windows = NULL;
  size_t skipped = 0;
  size_t win_count = make_windows(&c, 1, &windows, &skipped);
  TEST_ASSERT_INT_EQ(1, win_count);

  CrfWorker worker = {
    .train = NULL,
    .thread_id = 0,
    .gradient = anno_alloc(WEIGHTS, sizeof(double))
  };
  CrfObjective obj = {
    .train = {.windows = windows, .count = win_count, .skipped = 0, .threads = 1},
    .workers = &worker,
    .threads = NULL
  };
  worker.train = &obj.train;

  double *w_init = anno_alloc(WEIGHTS, sizeof(double));
  double *grad_init = anno_alloc(WEIGHTS, sizeof(double));
  double loss0 = crf_objective(w_init, grad_init, &obj);

  BruteForceCrf bf0;
  bf_crf_solve(w_init, &c, &bf0, 0);
  double gold_s0 = bf_score_and_accumulate_features(w_init, &c, gold, 0.0, NULL);
  TEST_ASSERT_DOUBLE_WITHIN(1e-11, bf0.logz - gold_s0, loss0);

  /* Step 2: Optimize CRF objective with L-BFGS (crf_train) */
  double *w_lbfgs = anno_alloc(WEIGHTS, sizeof(double));
  crf_train(&c, 1, w_lbfgs, 60, 1);

  double *grad_lbfgs = anno_alloc(WEIGHTS, sizeof(double));
  double loss_lbfgs = crf_objective(w_lbfgs, grad_lbfgs, &obj);
  TEST_ASSERT(loss_lbfgs < loss0);
  TEST_ASSERT(lbfgs_norm(WEIGHTS, grad_lbfgs) < 1e-4);

  /* Step 3: Optimize the exact same CRF objective with an independent solver and compare
     final objective value, gradient norm, and predicted posterior path probabilities */
  double *w_ref = anno_alloc(WEIGHTS, sizeof(double));
  double *grad_ref = anno_alloc(WEIGHTS, sizeof(double));
  double loss_ref = reference_gd_optimize(WEIGHTS, w_ref, crf_objective, &obj, 400);
  loss_ref = crf_objective(w_ref, grad_ref, &obj);

  TEST_ASSERT_DOUBLE_WITHIN(1e-6, loss_ref, loss_lbfgs);
  TEST_ASSERT(lbfgs_norm(WEIGHTS, grad_ref) < 1e-4);

  BruteForceCrf bf_lbfgs, bf_ref;
  bf_crf_solve(w_lbfgs, &c, &bf_lbfgs, 0);
  bf_crf_solve(w_ref, &c, &bf_ref, 0);
  for (int m = 0; m < bf_lbfgs.count; m++)
    TEST_ASSERT_DOUBLE_WITHIN(1e-4, bf_ref.probs[m], bf_lbfgs.probs[m]);

  /* Decoded path must recover the gold path */
  uint8_t decoded[9] = {0};
  crf_decode(w_lbfgs, &c, decoded);
  for (int t = 0; t < 9; t++) TEST_ASSERT_INT_EQ(gold[t], decoded[t]);

  work_free(&worker.work);
  free(worker.gradient);
  free(windows);
  free(w_init);
  free(grad_init);
  free(w_lbfgs);
  free(grad_lbfgs);
  free(w_ref);
  free(grad_ref);
  free_test_contig(&c);
}

static void test_integration_cwt_to_crf_end_to_end_and_threads(void) {
  Wavelets wavelets;
  wavelets_init(&wavelets);
  crf_init();

  /* Step 4: Generate real CWT features from sequence, connect to CRF, verify brute-force + L-BFGS */
  const char *seq1 =
    "GCGTACGTACGTACGTACGTACGTACGTAC"
    "ATGAAATAA"
    "GCGTACGTACGTACGTACGTACGTACGTACGT";
  const uint8_t gold1[9] = {1, 2, 3, 4, 4, 4, 1, 2, 3};

  const char *seq2 =
    "TCGTACGTACGTACGTACGTACGTACGTAC"
    "TTATCACAT"
    "TCGTACGTACGTACGTACGTACGTACGTACGT";
  const uint8_t gold2[9] = {9, 8, 7, 12, 12, 12, 9, 8, 7};

  Contig contigs_1t[2], contigs_4t[2];
  contigs_1t[0] = make_test_contig(seq1, 30, 9, gold1);
  contigs_1t[1] = make_test_contig(seq2, 30, 9, gold2);
  contigs_4t[0] = make_test_contig(seq1, 30, 9, gold1);
  contigs_4t[1] = make_test_contig(seq2, 30, 9, gold2);

  /* Compute real CWT powers with 1 thread vs 4 threads */
  for (int i = 0; i < 2; i++) {
    cwt_transform(&wavelets, contigs_1t[i].seq, contigs_1t[i].genome_length,
                  1, NULL, NULL, contigs_1t[i].cwt_power);
    cwt_transform(&wavelets, contigs_4t[i].seq, contigs_4t[i].genome_length,
                  4, NULL, NULL, contigs_4t[i].cwt_power);
    for (int k = 0; k < contigs_1t[i].genome_length * WAVE_COUNT; k++)
      TEST_ASSERT_DOUBLE_WITHIN(1e-7, contigs_1t[i].cwt_power[k], contigs_4t[i].cwt_power[k]);
  }

  /* Verify brute-force oracle on real CWT-connected contig */
  double *w_probe = anno_alloc(WEIGHTS, sizeof(double));
  rng_seed(9999);
  for (int f = 0; f < WEIGHTS; f++) w_probe[f] = rng_uniform(-0.2, 0.2);

  BruteForceCrf bf;
  bf_crf_solve(w_probe, &contigs_1t[0], &bf, 0);
  double dp_logz = 0.0;
  double dp_node[BF_MAX_T][STATES];
  double dp_edge[BF_MAX_T][STATES][STATES];
  double dp_grad[WEIGHTS];
  dp_extract_marginals(w_probe, &contigs_1t[0], &dp_logz, dp_node, dp_edge, dp_grad);
  TEST_ASSERT_DOUBLE_WITHIN(1e-11, bf.logz, dp_logz);

  /* Step 5: Train and decode with 1 thread vs multi-thread and compare stage by stage */
  double *w_1t = anno_alloc(WEIGHTS, sizeof(double));
  double *w_2t = anno_alloc(WEIGHTS, sizeof(double));
  crf_train(contigs_1t, 2, w_1t, 50, 1);
  crf_train(contigs_4t, 2, w_2t, 50, 2);

  for (int f = 0; f < WEIGHTS; f++)
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, w_1t[f], w_2t[f]);

  for (int i = 0; i < 2; i++) {
    uint8_t path_1t[9] = {0}, path_2t[9] = {0};
    crf_decode(w_1t, &contigs_1t[i], path_1t);
    crf_decode(w_2t, &contigs_4t[i], path_2t);
    const uint8_t *expected_gold = (i == 0) ? gold1 : gold2;
    for (int t = 0; t < 9; t++) {
      TEST_ASSERT_INT_EQ( path_1t[t], path_2t[t] );
      TEST_ASSERT_INT_EQ( expected_gold[t], path_1t[t] );
    }
  }

  free(w_probe);
  free(w_1t);
  free(w_2t);
  for (int i = 0; i < 2; i++) {
    free_test_contig(&contigs_1t[i]);
    free_test_contig(&contigs_4t[i]);
  }
  wavelets_free(&wavelets);
}

int main(void) {
  printf("=== Running CRF, CWT, L-BFGS & Integration Test Suite ===\n");

  /* 1. CRF Module */
  RUN_TEST(test_crf_context_coordinates_and_stream);
  RUN_TEST(test_crf_reverse_complement_emissions);
  RUN_TEST(test_crf_k2_exhaustive_t1_to_t5);
  RUN_TEST(test_crf_grammar_exhaustive_both_strands);
  RUN_TEST(test_crf_zero_scores);
  RUN_TEST(test_crf_length_one);
  RUN_TEST(test_crf_asymmetric_transitions);
  RUN_TEST(test_crf_tied_optimal_paths);
  RUN_TEST(test_crf_forbidden_transitions);
  RUN_TEST(test_crf_no_valid_path_and_nonfinite);
  RUN_TEST(test_crf_extreme_scores_numerical_stability);
  RUN_TEST(test_crf_shift_invariance);
  RUN_TEST(test_crf_objective_sign_regularization_and_scaling);

  /* 2. CWT Module */
  RUN_TEST(test_cwt_zero_signal);
  RUN_TEST(test_cwt_linearity);
  RUN_TEST(test_cwt_center_impulse);
  RUN_TEST(test_cwt_boundary_impulses);
  RUN_TEST(test_cwt_direct_vs_fft_and_window_boundaries);
  RUN_TEST(test_cwt_sinusoid_scale_response);
  RUN_TEST(test_cwt_repeat_calls_and_threads);

  /* 3. L-BFGS Module */
  RUN_TEST(test_lbfgs_crf_directional_derivatives);
  RUN_TEST(test_lbfgs_quadratic_isotropic_diagonal_coupled);
  RUN_TEST(test_lbfgs_rosenbrock);
  RUN_TEST(test_lbfgs_termination_states);
  RUN_TEST(test_lbfgs_two_loop_history_and_wolfe);

  /* 4. Integration (CWT + CRF + L-BFGS) */
  RUN_TEST(test_integration_fixed_features_crf_lbfgs_and_reference_solver);
  RUN_TEST(test_integration_cwt_to_crf_end_to_end_and_threads);

  printf("---------------------------------------------------------\n");
  printf("%d Tests %d Failures\n", g_tests_run, g_tests_failed);
  printf("%s\n", g_tests_failed == 0 ? "OK" : "FAIL");
  return g_tests_failed == 0 ? 0 : 1;
}
