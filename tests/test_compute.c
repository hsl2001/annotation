#define main anno_cli_main
#include "../anno.c"
#undef main
#include <time.h>

static double wall_seconds(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static int check_close(const char *name, const double *actual, const double *ref,
                       size_t n, double atol, double rtol) {
  double max_abs = 0.0, max_rel = 0.0;
  size_t failures = 0;
  for (size_t i = 0; i < n; i++) {
    double diff = fabs(actual[i] - ref[i]);
    double limit = atol + rtol * fabs(ref[i]);
    if (diff > max_abs) max_abs = diff;
    if (fabs(ref[i]) > 1e-12) {
      double rel = diff / fabs(ref[i]);
      if (rel > max_rel) max_rel = rel;
    }
    if (!isfinite(actual[i]) || diff > limit) failures++;
  }
  printf("  [%-16s] n=%-8zu max_abs=%-10.3e max_rel=%-10.3e -> %s\n",
         name, n, max_abs, max_rel, failures == 0 ? "PASS" : "FAIL");
  return failures == 0 ? 0 : 1;
}

static int check_close_f32(const char *name, const float *actual, const float *ref,
                           size_t n, double atol, double rtol) {
  double max_abs = 0.0, max_rel = 0.0;
  size_t failures = 0;
  for (size_t i = 0; i < n; i++) {
    double diff = fabs((double)actual[i] - (double)ref[i]);
    double limit = atol + rtol * fabs((double)ref[i]);
    if (diff > max_abs) max_abs = diff;
    if (fabs((double)ref[i]) > 1e-7) {
      double rel = diff / fabs((double)ref[i]);
      if (rel > max_rel) max_rel = rel;
    }
    if (!isfinite(actual[i]) || diff > limit) failures++;
  }
  printf("  [%-16s] n=%-8zu max_abs=%-10.3e max_rel=%-10.3e -> %s\n",
         name, n, max_abs, max_rel, failures == 0 ? "PASS" : "FAIL");
  return failures == 0 ? 0 : 1;
}

int main(int argc, char **argv) {
  const char *mode = argc >= 3 ? argv[1] : "--self";
  const char *ref_path = argc >= 3 ? argv[2] : NULL;

#ifdef _OPENACC
  const char *target = "GPU (OpenACC)";
#else
  const char *target = "CPU";
#endif

  Wavelets wavelets;
  wavelets_init(&wavelets);
  crf_init();

  int count = 0;
  Contig *contigs = read_fasta("truth.nip.chr5.3m.fasta", &count);
  /* Use a 262,144 bp (4-window) slice for fast deterministic verification */
  if (contigs[0].length > 262144) {
    contigs[0].length = 262144;
    contigs[0].genome_length = 262144;
  }
  int len = contigs[0].length;
  size_t power_count = (size_t)len * WAVE_COUNT;
  size_t complex_count = (size_t)len * CWT_CHANNELS;

  float *powers = anno_alloc(power_count, sizeof(float));
  float *complex_vals = anno_alloc(complex_count, sizeof(float));

  double t0 = wall_seconds();
  cwt_transform(&wavelets, contigs[0].seq, len, 8, NULL, complex_vals, powers);
  double t_cwt = wall_seconds() - t0;

  contigs[0].cwt_power = powers;
  contigs[0].label = anno_alloc((size_t)len, 1);
  int skipped = 0;
  label_cds(contigs, count, "truth.nip.chr5.3m.gff3", &skipped);

  /* Single objective + gradient evaluation */
  CrfWindow *windows = NULL;
  size_t window_count = make_windows(contigs, count, &windows);
  CrfObjective objective = {
    .train = {.windows = windows, .count = window_count, .threads = 8}
  };
  atomic_init(&objective.train.next_instance, 0);
#ifdef _OPENACC
  crf_gpu_init(&objective.gpu, contigs, count, windows, window_count);
#else
  objective.workers = anno_alloc((size_t)objective.train.threads, sizeof(*objective.workers));
  objective.threads = anno_alloc((size_t)objective.train.threads, sizeof(*objective.threads));
  for (int t = 0; t < objective.train.threads; t++)
    objective.workers[t] = (CrfWorker){
      .train = &objective.train,
      .gradient = anno_alloc(WEIGHTS, sizeof(double))
    };
#endif

  double *w_init = anno_alloc(WEIGHTS, sizeof(double));
  double *grad = anno_alloc(WEIGHTS, sizeof(double));
  for (int f = 0; f < WEIGHTS; f++) w_init[f] = 0.001 * ((f % 7) - 3);

  double t1 = wall_seconds();
  double loss = crf_objective(w_init, grad, &objective);
  double t_obj = wall_seconds() - t1;

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

  /* L-BFGS 2-iteration training + Viterbi decode */
  double *w_trained = anno_alloc(WEIGHTS, sizeof(double));
  uint8_t *decoded = anno_alloc((size_t)len, 1);
  double t2 = wall_seconds();
  crf_train(contigs, count, w_trained, 2, 8);
  double t_train = wall_seconds() - t2;

  double t3 = wall_seconds();
  crf_decode(w_trained, &contigs[0], decoded);
  double t_decode = wall_seconds() - t3;

  printf("=== %s Execution Summary (N=%d bp) ===\n", target, len);
  printf("  CWT (power+complex): %.4fs | CRF obj+grad: %.4fs | Train(2 it): %.4fs | Decode: %.4fs\n",
         t_cwt, t_obj, t_train, t_decode);

  int status = 0;
  if (!isfinite(loss)) {
    fprintf(stderr, "Error: non-finite loss %g\n", loss);
    status = 1;
  }

  if (!strcmp(mode, "--write") && ref_path) {
    FILE *fp = fopen(ref_path, "wb");
    if (!fp) anno_fail("Cannot write reference file: %s", ref_path);
    if (fwrite(powers, sizeof(float), power_count, fp) != power_count ||
        fwrite(complex_vals, sizeof(float), complex_count, fp) != complex_count ||
        fwrite(&loss, sizeof(double), 1, fp) != 1 ||
        fwrite(grad, sizeof(double), WEIGHTS, fp) != WEIGHTS ||
        fwrite(w_trained, sizeof(double), WEIGHTS, fp) != WEIGHTS ||
        fwrite(decoded, 1, (size_t)len, fp) != (size_t)len ||
        fclose(fp))
      anno_fail("Failed writing reference file: %s", ref_path);
    printf("  Saved reference outputs to %s\n", ref_path);
  } else if (!strcmp(mode, "--verify") && ref_path) {
    float *ref_powers = anno_alloc(power_count, sizeof(float));
    float *ref_complex = anno_alloc(complex_count, sizeof(float));
    double ref_loss = 0.0;
    double *ref_grad = anno_alloc(WEIGHTS, sizeof(double));
    double *ref_weights = anno_alloc(WEIGHTS, sizeof(double));
    uint8_t *ref_decoded = anno_alloc((size_t)len, 1);

    FILE *fp = fopen(ref_path, "rb");
    if (!fp) anno_fail("Cannot open reference file: %s", ref_path);
    if (fread(ref_powers, sizeof(float), power_count, fp) != power_count ||
        fread(ref_complex, sizeof(float), complex_count, fp) != complex_count ||
        fread(&ref_loss, sizeof(double), 1, fp) != 1 ||
        fread(ref_grad, sizeof(double), WEIGHTS, fp) != WEIGHTS ||
        fread(ref_weights, sizeof(double), WEIGHTS, fp) != WEIGHTS ||
        fread(ref_decoded, 1, (size_t)len, fp) != (size_t)len)
      anno_fail("Corrupted reference file: %s", ref_path);
    fclose(fp);

    printf("=== Numerical Verification vs Reference (|y - y_ref| <= atol + rtol*|y_ref|) ===\n");
    status |= check_close_f32("cwt_powers", powers, ref_powers, power_count, 1e-5, 1e-5);
    status |= check_close_f32("cwt_complex", complex_vals, ref_complex, complex_count, 1e-5, 1e-5);
    status |= check_close("crf_loss", &loss, &ref_loss, 1, 1e-6, 1e-9);
    status |= check_close("crf_gradient", grad, ref_grad, WEIGHTS, 1e-5, 1e-7);
    status |= check_close("trained_weights", w_trained, ref_weights, WEIGHTS, 1e-5, 1e-6);

    size_t path_diff = 0;
    for (int i = 0; i < len; i++)
      if (decoded[i] != ref_decoded[i]) path_diff++;
    printf("  [%-16s] n=%-8d mismatches=%-8zu -> %s\n",
           "viterbi_decode", len, path_diff, path_diff == 0 ? "PASS" : "FAIL");
    if (path_diff != 0) status |= 1;

    free(ref_powers), free(ref_complex), free(ref_grad), free(ref_weights), free(ref_decoded);
  }

  free(powers), free(complex_vals), free(w_init), free(grad), free(w_trained), free(decoded);
  wavelets_free(&wavelets);
  free_contigs(contigs, count);
  return status;
}
