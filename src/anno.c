/**
 * @file anno.c
 * @brief Command-line entry point for CRF training and genome prediction.
 */
#define _POSIX_C_SOURCE 200809L

#include "include/anno_internal.h"

#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/* ------------------------------------------------------------------------ */
/* Command line                                                             */
/* ------------------------------------------------------------------------ */

/**
 * @brief Print command-line usage to standard error.
 * @return None.
 */
static void usage(void) {
  fprintf(stderr,
          "Usage: anno [options] genome.fasta[.gz] [annotation.gff3[.gz]]\n"
          "  with an annotation: train a CRF with L-BFGS and write it to -m\n"
          "  without:            predict genes (GFF3 on stdout unless -o)\n"
          "Options:\n"
          "  -p N       threads [1]\n"
          "  -i N       training iterations [%d]\n"
          "  -m FILE    CRF model file [anno.crf]\n"
          "  -o PREFIX  write PREFIX.gff3 and PREFIX.{class}.bedgraph\n"
          "  -h         this help\n",
          ANNO_DEFAULT_ITERATIONS);
}

/**
 * @brief Parse a positive integer option value.
 * @param[in] text Numeric option text.
 * @param[in] name Option name used in diagnostics.
 * @return Parsed positive integer.
 */
static int positive(const char *text, const char *name) {
  char *end = NULL;
  errno = 0;
  long value = strtol(text, &end, 10);
  if (errno || end == text || *end || value < 1 || value > INT_MAX)
    anno_fail("invalid %s: %s", name, text);
  return (int)value;
}

/**
 * @brief Run the anno command-line application.
 * @param[in] argc Number of command-line arguments.
 * @param[in] argv Command-line argument vector.
 * @return EXIT_SUCCESS on success, otherwise EXIT_FAILURE.
 */
int main(int argc, char **argv) {
  int threads = 1, iterations = ANNO_DEFAULT_ITERATIONS;
  const char *model_path = "anno.crf", *prefix = NULL;
  int opt;
  opterr = 0;
  while ((opt = getopt(argc, argv, "hp:i:m:o:")) != -1) {
    if (opt == 'h') {
      usage();
      return EXIT_SUCCESS;
    }
    if (opt == '?') {
      usage();
      return EXIT_FAILURE;
    }
    if (opt == 'p')
      threads = positive(optarg, "thread count");
    if (opt == 'i')
      iterations = positive(optarg, "iteration count");
    if (opt == 'm')
      model_path = optarg;
    if (opt == 'o')
      prefix = optarg;
  }
#ifdef _OPENMP
  omp_set_num_threads(threads);
#endif
  int a = optind;
  if (argc - a < 1 || argc - a > 2) {
    usage();
    return EXIT_FAILURE;
  }
  int count = 0;
  AnnoSequence *sequences = anno_read_fasta(argv[a], &count);
  if (argc - a == 2) {
    double started = anno_now_seconds();
    anno_read_labels(argv[a + 1], sequences, count);
    CrfModel *model = train_crf(sequences, count, iterations, threads);
    int status = crf_save_model(model_path, model);
    if (status)
      anno_fail("cannot write CRF model %s: %s", model_path, strerror(status));
    anno_log("trained in %.1fs; CRF model written to %s", anno_now_seconds() - started,
        model_path);
    free(model);
  } else {
    CrfModel *model = anno_xcalloc(1, sizeof(*model));
    int status = crf_load_model(model_path, model);
    if (status == EINVAL)
      anno_fail("invalid or incompatible CRF model %s; retrain with this version of anno",
                model_path);
    if (status)
      anno_fail("cannot read CRF model %s: %s", model_path, strerror(status));
    FILE *gff = stdout;
    if (prefix) {
      char path[ANNO_PATH_BUFFER];
      if (snprintf(path, sizeof(path), "%s.gff3", prefix) >= (int)sizeof(path))
        anno_fail("output prefix too long");
      if (!(gff = fopen(path, "w")))
        anno_fail("cannot create %s", path);
    }
    fputs("##gff-version 3\n", gff);
    unsigned long genes =
        predict_genome(model, sequences, count, gff, prefix);
    if (fflush(gff) || ferror(gff) || (gff != stdout && fclose(gff)))
      anno_fail("cannot write GFF3");
    anno_log("predicted %lu genes", genes);
    free(model);
  }
  anno_sequences_free(sequences, count);
  return EXIT_SUCCESS;
}
