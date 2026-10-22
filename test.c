#ifdef __APPLE__
#define _DARWIN_C_SOURCE
#endif
#define _POSIX_C_SOURCE 200809L

#include <complex.h>
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

const int wave_sizes[WAVE_COUNT] = {CWT_SCALES};
#define CWT_IMPLEMENTATION
#include "cwt.h"
#define CRF_IMPLEMENTATION
#include "crf.h"

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

static void init_test_gff_cwt_model(GffCwtModel *m) {
  rng_seed(987654321ULL);
  for (int r = 0; r < GFF_REGIONS; r++)
    for (int ch = 0; ch < CWT_CHANNELS; ch++)
      m->region[r][ch] = (float)rng_uniform(-0.8, 0.8);
  for (int rev = 0; rev < 2; rev++)
    for (int k = 0; k < GFF_BOUND_KINDS; k++)
      for (int s = 0; s < GFF_CWT_WINDOW; s++)
        for (int ch = 0; ch < CWT_CHANNELS; ch++)
          m->bound[rev][k][s][ch] = (float)rng_uniform(-0.5, 0.5);
}

static Contig make_test_contig(const char *seq, const uint8_t *labels) {
  int len = (int)strlen(seq);
  Contig c = {
    .name = strdup("test"),
    .seq = strdup(seq),
    .length = len,
    .cwt = anno_alloc((size_t)len * CWT_CHANNELS, sizeof(float)),
    .label = labels ? anno_alloc((size_t)len, sizeof(uint8_t)) : NULL
  };
  for (int i = 0; i < len * CWT_CHANNELS; i++)
    c.cwt[i] = (float)(0.25 * sin(0.17 * (i + 1)) - 0.15 * cos(0.11 * (i + 3)));
  if (labels) memcpy(c.label, labels, (size_t)len);
  return c;
}

static void free_test_contig(Contig *c) {
  free(c->name);
  free(c->seq);
  free(c->cwt);
  free(c->label);
  memset(c, 0, sizeof(*c));
}

/* 1. FASTA -> CWT Tests */

static void cwt_direct_reference(const Wavelets *wavelets, const char *seq, int length,
                                 float *complex_values) {
  for (int pos = 0; pos < length; pos++)
    for (int sc = 0; sc < WAVE_COUNT; sc++) {
      int width = wavelets->widths[sc], radius = width / 2;
      double complex sum = 0.0;
      for (int tap = 0; tap < width; tap++) {
        int q = pos + tap - radius;
        if (q >= 0 && q < length)
          sum += cwt_base_signal(seq[q]) * wavelets->kernel[sc][tap];
      }
      size_t idx = (size_t)pos * CWT_CHANNELS + 2 * (size_t)sc;
      complex_values[idx] = (float)creal(sum);
      complex_values[idx + 1] = (float)cimag(sum);
    }
}

static void test_cwt_transform_and_chunking(void) {
  Wavelets wavelets;
  wavelets_init(&wavelets);
  rng_seed(404);

  char n_seq[65];
  memset(n_seq, 'N', 64);
  n_seq[64] = '\0';
  float n_cwt[64 * CWT_CHANNELS];
  cwt_transform(&wavelets, n_seq, 64, 1, n_cwt);
  for (int i = 0; i < 64 * CWT_CHANNELS; i++)
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.0, n_cwt[i]);

  const int lengths[] = {1, 7, 64, CWT_WINDOW_SIZE - 1, CWT_WINDOW_SIZE, CWT_WINDOW_SIZE + 1, 2055};
  const char alphabet[] = "ACGTN";
  for (size_t idx = 0; idx < sizeof(lengths) / sizeof(lengths[0]); idx++) {
    int L = lengths[idx];
    char *seq = anno_alloc((size_t)L + 1, 1);
    for (int i = 0; i < L; i++) {
      int r = (int)rng_uniform(0.0, 100.0);
      seq[i] = (r < 92) ? alphabet[r & 3] : 'N';
    }
    float *fft_1t = anno_alloc((size_t)L * CWT_CHANNELS, sizeof(float));
    float *fft_4t = anno_alloc((size_t)L * CWT_CHANNELS, sizeof(float));
    float *ref = anno_alloc((size_t)L * CWT_CHANNELS, sizeof(float));
    float *chunked = anno_alloc((size_t)L * CWT_CHANNELS, sizeof(float));

    cwt_transform(&wavelets, seq, L, 1, fft_1t);
    cwt_transform(&wavelets, seq, L, 4, fft_4t);
    cwt_direct_reference(&wavelets, seq, L, ref);

    int chunk_size = 701, halo = wavelets.max_width / 2;
    for (int start = 0; start < L;) {
      int chunk_len = L - start > chunk_size ? chunk_size : L - start;
      int src_start = start > halo ? start - halo : 0;
      int core_end = start + chunk_len;
      int src_end = L - core_end > halo ? core_end + halo : L;
      int src_len = src_end - src_start;
      float *buf = anno_alloc((size_t)src_len * CWT_CHANNELS, sizeof(float));
      cwt_transform(&wavelets, seq + src_start, src_len, 1, buf);
      memcpy(chunked + (size_t)start * CWT_CHANNELS,
             buf + (size_t)(start - src_start) * CWT_CHANNELS,
             (size_t)chunk_len * CWT_CHANNELS * sizeof(float));
      free(buf);
      start = core_end;
    }

    for (int i = 0; i < L * CWT_CHANNELS; i++) {
      TEST_ASSERT_DOUBLE_WITHIN(1e-5, ref[i], fft_1t[i]);
      TEST_ASSERT_DOUBLE_WITHIN(1e-6, fft_1t[i], fft_4t[i]);
      TEST_ASSERT_DOUBLE_WITHIN(1e-5, fft_1t[i], chunked[i]);
    }
    free(seq), free(fft_1t), free(fft_4t), free(ref), free(chunked);
  }
  wavelets_free(&wavelets);
}

/* 2. GFF -> CWT Model, Adjoint, and Accumulator Tests */

static double eval_hard_residual(const Contig *c, const GffCwtModel *model) {
  double *reg_m = anno_alloc((size_t)c->length * GFF_REGIONS, sizeof(double));
  double *bnd_m = anno_alloc((size_t)c->length * (2 * GFF_BOUND_KINDS), sizeof(double));
  double *pred = anno_alloc((size_t)c->length * CWT_CHANNELS, sizeof(double));
  for (int i = 0; i < c->length; i++) {
    reg_m[i * GFF_REGIONS + gff_region(c->label[i])] = 1.0;
    if (i > 0) {
      int kind = boundary_kind(c->label[i - 1], c->label[i]);
      if (kind >= 0) {
        int rev = boundary_rev(c->label[i - 1], c->label[i]);
        bnd_m[i * (2 * GFF_BOUND_KINDS) + rev * GFF_BOUND_KINDS + kind] = 1.0;
      }
    }
  }
  gff_cwt_predict_soft(reg_m, bnd_m, c->length, model, pred);
  double energy = 0.0;
  for (int i = 0; i < c->length; i++) {
    if (base_index(c->seq[i]) < 0) continue;
    for (int ch = 0; ch < CWT_CHANNELS; ch++) {
      double d = pred[i * CWT_CHANNELS + ch] - (double)c->cwt[i * CWT_CHANNELS + ch];
      energy += 0.5 * d * d;
    }
  }
  free(reg_m), free(bnd_m), free(pred);
  return energy;
}

static void test_gff_cwt_model_and_adjoint(void) {
  Wavelets wavelets;
  wavelets_init(&wavelets);
  rng_seed(4242);

  GffCwtModel model;
  init_test_gff_cwt_model(&model);
  const int L = 80;
  double *reg_m = anno_alloc((size_t)L * GFF_REGIONS, sizeof(double));
  double *bnd_m = anno_alloc((size_t)L * (2 * GFF_BOUND_KINDS), sizeof(double));
  double *pred_cwt = anno_alloc((size_t)L * CWT_CHANNELS, sizeof(double));
  double *res_cwt = anno_alloc((size_t)L * CWT_CHANNELS, sizeof(double));
  double *grad_reg = anno_alloc((size_t)L * GFF_REGIONS, sizeof(double));
  double *grad_bnd = anno_alloc((size_t)L * (2 * GFF_BOUND_KINDS), sizeof(double));

  for (int i = 0; i < L * GFF_REGIONS; i++) reg_m[i] = rng_uniform(0.0, 1.0);
  for (int i = 2 * GFF_BOUND_KINDS; i < L * (2 * GFF_BOUND_KINDS); i++)
    bnd_m[i] = rng_uniform(0.0, 0.5);
  for (int i = 0; i < L * CWT_CHANNELS; i++) res_cwt[i] = rng_uniform(-1.0, 1.0);

  gff_cwt_predict_soft(reg_m, bnd_m, L, &model, pred_cwt);
  gff_cwt_residual_adjoint(res_cwt, L, &model, grad_reg, grad_bnd);

  double lhs = 0.0, rhs = 0.0;
  for (int i = 0; i < L * CWT_CHANNELS; i++) lhs += res_cwt[i] * pred_cwt[i];
  for (int i = 0; i < L * GFF_REGIONS; i++) rhs += grad_reg[i] * reg_m[i];
  for (int i = 0; i < L * (2 * GFF_BOUND_KINDS); i++) rhs += grad_bnd[i] * bnd_m[i];
  TEST_ASSERT_DOUBLE_WITHIN(1e-11, lhs, rhs);

  const char *seq = "ACGTACGTACGTATGTTGTAATGCAACGTACGT";
  uint8_t labels[33] = {0};
  const uint8_t gene[9] = {1, 2, 3, 4, 4, 4, 1, 2, 3};
  memcpy(labels + 12, gene, sizeof(gene));
  Contig c = make_test_contig(seq, labels);
  cwt_transform(&wavelets, c.seq, c.length, 1, c.cwt);

  GffCwtModel zero_model = {0}, full_model, split_model;
  double zero_res = eval_hard_residual(&c, &zero_model);

  GffCwtAccumulator acc;
  gff_cwt_accumulator_init(&acc);
  gff_cwt_accumulate_range(&acc, &c, 0, c.length, c.cwt);
  gff_cwt_accumulator_finish(&acc, &full_model);
  double fit_res = eval_hard_residual(&c, &full_model);
  TEST_ASSERT(fit_res < zero_res);

  gff_cwt_accumulator_init(&acc);
  gff_cwt_accumulate_range(&acc, &c, 0, 15, c.cwt);
  gff_cwt_accumulate_range(&acc, &c, 15, c.length - 15, c.cwt + 15 * CWT_CHANNELS);
  gff_cwt_accumulator_finish(&acc, &split_model);
  TEST_ASSERT(memcmp(&full_model, &split_model, sizeof(full_model)) == 0);

  free(reg_m), free(bnd_m), free(pred_cwt), free(res_cwt), free(grad_reg), free(grad_bnd);
  free_test_contig(&c);
  wavelets_free(&wavelets);
}

/* 3. CRF Grammar, Decoding, Residual Gradient, and Self-Training Tests */

static void test_crf_reverse_complement_and_decode(void) {
  crf_init();
  const char *seq = "AAACCCGGGTTTACGT";
  char rev_seq[17];
  for (int i = 0; i < 16; i++)
    rev_seq[i] = "TGCA"[base_index(seq[15 - i])];
  rev_seq[16] = '\0';

  Contig cf = make_test_contig(seq, NULL);
  Contig cr = make_test_contig(rev_seq, NULL);

  double *w = anno_alloc(WEIGHTS, sizeof(double));
  rng_seed(202);
  for (int f = 0; f < WEIGHTS; f++) w[f] = rng_uniform(-0.5, 0.5);

  for (int i = 0; i < 16; i++) {
    Site sf = site_at(&cf, i);
    Site sr = site_at(&cr, 15 - i);
    TEST_ASSERT_INT_EQ(sf.k3[0], sr.k3[1]);
    TEST_ASSERT_INT_EQ(sf.k3[1], sr.k3[0]);
    for (int f = 0; f < 3; f++) {
      TEST_ASSERT_DOUBLE_WITHIN(1e-12, w[emission_index(&sf, 1 + f)], w[emission_index(&sr, 7 + f)]);
      TEST_ASSERT_DOUBLE_WITHIN(1e-12, w[emission_index(&sf, 4 + f)], w[emission_index(&sr, 10 + f)]);
    }
  }

  uint8_t path[16] = {0};
  crf_decode(w, &cf, path);
  w[0] = NAN;
  TEST_ASSERT_FAILS(crf_decode(w, &cf, path));

  free(w);
  free_test_contig(&cf);
  free_test_contig(&cr);
}

static void test_crf_residual_gradient_and_objective(void) {
  crf_init();
  GffCwtModel gmodel;
  init_test_gff_cwt_model(&gmodel);
  crf_set_gff_cwt(&gmodel);
  rng_seed(1337);

  double *w = anno_alloc(WEIGHTS, sizeof(double));
  for (int f = 0; f < WEIGHTS; f++) w[f] = rng_uniform(-0.3, 0.3);

  const char *genome =
    "ACGTACGTACGTACGTACGTACGTACGTAC"
    "ATGTTGTAA"
    "ACGTACGTACGTACGTACGTACGTACGTACGT";
  Contig c = make_test_contig(genome, NULL);

  Work wk = {0};
  double grad[WEIGHTS] = {0}, dummy[WEIGHTS] = {0};
  double res_weight = 0.25;
  double loss0 = window_residual_gradient(w, &c, &wk, grad, res_weight, 25, 19, 3, 12);
  TEST_ASSERT(loss0 > 0.0);

  double direction[WEIGHTS], dot = 0.0, orig[WEIGHTS];
  memcpy(orig, w, sizeof(orig));
  for (int f = 0; f < WEIGHTS; f++) {
    direction[f] = rng_uniform(-0.3, 0.3);
    dot += grad[f] * direction[f];
  }
  for (int f = 0; f < WEIGHTS; f++) w[f] = orig[f] + 1e-6 * direction[f];
  double loss_plus = window_residual_gradient(w, &c, &wk, dummy, res_weight, 25, 19, 3, 12);
  for (int f = 0; f < WEIGHTS; f++) w[f] = orig[f] - 1e-6 * direction[f];
  double loss_minus = window_residual_gradient(w, &c, &wk, dummy, res_weight, 25, 19, 3, 12);
  memcpy(w, orig, sizeof(orig));
  TEST_ASSERT_DOUBLE_WITHIN(1e-6, dot, (loss_plus - loss_minus) / 2e-6);

  CrfWindow win = {&c, 25, 19, 3, 12};
  CrfWorker worker = {.gradient = anno_alloc(WEIGHTS, sizeof(double))};
  CrfObjective obj = {
    .train = {.windows = &win, .count = 1, .threads = 1},
    .workers = &worker
  };
  worker.train = &obj.train;
  double obj_grad[WEIGHTS] = {0}, obj_dot = 0.0;
  crf_objective(w, obj_grad, &obj);
  for (int f = 0; f < WEIGHTS; f++) obj_dot += obj_grad[f] * direction[f];
  for (int f = 0; f < WEIGHTS; f++) w[f] = orig[f] + 1e-6 * direction[f];
  double obj_plus = crf_objective(w, dummy, &obj);
  for (int f = 0; f < WEIGHTS; f++) w[f] = orig[f] - 1e-6 * direction[f];
  double obj_minus = crf_objective(w, dummy, &obj);
  memcpy(w, orig, sizeof(orig));
  TEST_ASSERT_DOUBLE_WITHIN(1e-6, obj_dot, (obj_plus - obj_minus) / 2e-6);

  work_free(&worker.work);
  free(worker.gradient);
  work_free(&wk);
  free_test_contig(&c);
  crf_set_gff_cwt(NULL);
  free(w);
}

static void test_self_training_end_to_end_and_threads(void) {
  Wavelets wavelets;
  wavelets_init(&wavelets);
  crf_init();

  const char *seq1 =
    "ACGTACGTACGTACGTACGTACGTACGTAC"
    "ATGAAATAA"
    "ACGTACGTACGTACGTACGTACGTACGTACGT";
  const char *seq2 =
    "ACGTACGTACGTACGTACGTACGTACGTAC"
    "TTATTTCAT"
    "ACGTACGTACGTACGTACGTACGTACGTACGT";
  const uint8_t labels1[9] = {1, 2, 3, 1, 2, 3, 1, 2, 3};
  const uint8_t labels2[9] = {9, 8, 7, 9, 8, 7, 9, 8, 7};

  Contig labeled[2] = {make_test_contig(seq1, NULL), make_test_contig(seq2, NULL)};
  for (int i = 0; i < 2; i++) {
    cwt_transform(&wavelets, labeled[i].seq, labeled[i].length, 1, labeled[i].cwt);
    labeled[i].label = anno_alloc((size_t)labeled[i].length, 1);
  }
  memcpy(labeled[0].label + 30, labels1, 9);
  memcpy(labeled[1].label + 30, labels2, 9);

  GffCwtAccumulator acc;
  GffCwtModel gmodel;
  gff_cwt_accumulator_init(&acc);
  for (int i = 0; i < 2; i++) {
    gff_cwt_accumulate_range(&acc, &labeled[i], 0, labeled[i].length, labeled[i].cwt);
    free(labeled[i].label);
    labeled[i].label = NULL;
  }
  gff_cwt_accumulator_finish(&acc, &gmodel);
  crf_set_gff_cwt(&gmodel);

  Work wk = {0};
  double dummy_grad[WEIGHTS] = {0};
  double *w_1t = anno_alloc(WEIGHTS, sizeof(double));
  double *w_2t = anno_alloc(WEIGHTS, sizeof(double));
  double init_res = window_residual_gradient(w_1t, &labeled[0], &wk, dummy_grad,
                                             CRF_RESIDUAL_WEIGHT, 0, labeled[0].length,
                                             0, labeled[0].length);

  crf_train(labeled, 2, w_1t, 15, 1);
  crf_train(labeled, 2, w_2t, 15, 2);

  double trained_res = window_residual_gradient(w_1t, &labeled[0], &wk, dummy_grad,
                                                CRF_RESIDUAL_WEIGHT, 0, labeled[0].length,
                                                0, labeled[0].length);
  TEST_ASSERT(trained_res < init_res);
  for (int f = 0; f < WEIGHTS; f++)
    TEST_ASSERT_DOUBLE_WITHIN(1e-8, w_1t[f], w_2t[f]);

  work_free(&wk);
  free(w_1t), free(w_2t);
  free_test_contig(&labeled[0]);
  free_test_contig(&labeled[1]);
  crf_set_gff_cwt(NULL);
  wavelets_free(&wavelets);
}

/* 4. L-BFGS Optimizer Tests */

static double rosenbrock_eval(const double *x, double *g, void *ctx) {
  (void)ctx;
  double a = 1.0 - x[0], b = x[1] - x[0] * x[0];
  g[0] = -2.0 * a - 400.0 * x[0] * b;
  g[1] = 200.0 * b;
  return a * a + 100.0 * b * b;
}

static void test_lbfgs_rosenbrock_and_validation(void) {
  double x[2] = {-1.2, 1.0};
  int completed = 0;
  LbfgsStatus st = lbfgs_minimize(2, x, rosenbrock_eval, NULL, NULL, 200, &completed);
  TEST_ASSERT(st == LBFGS_CONVERGED || st == LBFGS_STOPPING_CRITERION);
  TEST_ASSERT_DOUBLE_WITHIN(1e-3, 1.0, x[0]);
  TEST_ASSERT_DOUBLE_WITHIN(1e-3, 1.0, x[1]);
  TEST_ASSERT_INT_EQ(LBFGS_INVALID_ARGUMENT,
                     lbfgs_minimize(0, x, rosenbrock_eval, NULL, NULL, 10, &completed));
}

int main(void) {
  printf("=== Running CWT, GFF->CWT, CRF & L-BFGS Test Suite ===\n");
  RUN_TEST(test_cwt_transform_and_chunking);
  RUN_TEST(test_gff_cwt_model_and_adjoint);
  RUN_TEST(test_crf_reverse_complement_and_decode);
  RUN_TEST(test_crf_residual_gradient_and_objective);
  RUN_TEST(test_self_training_end_to_end_and_threads);
  RUN_TEST(test_lbfgs_rosenbrock_and_validation);
  printf("------------------------------------------------------\n");
  printf("%d Tests %d Failures\n", g_tests_run, g_tests_failed);
  printf("%s\n", g_tests_failed == 0 ? "OK" : "FAIL");
  return g_tests_failed == 0 ? 0 : 1;
}
