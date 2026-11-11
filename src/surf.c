/**
 * @file anno.c
 * @brief CLI entry-point for CRF training and genome prediction.
 */
#include "include/surf.h"
#include "include/crf.h"
#include "include/utils.h"

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
          "Usage: surf [options] genome.fasta[.gz] [annotation.gff3[.gz]]\n"
          "  with an annotation: train a CRF with L-BFGS and write it to -m\n"
          "  without:            predict genes (GFF3 on stdout unless -o)\n"
          "Options:\n"
          "  -p N       threads [1]\n"
          "  -i N       training iterations [%d]\n"
          "  -m FILE    CRF model file [surf.crf]\n"
          "  -o PREFIX  write PREFIX.gff3 and PREFIX.{class}.bedgraph\n"
          "  -h         this help\n",
          DEFAULT_ITERATIONS);
}

/**
 * @brief Parse a positive integer option value.
 * @param[in] text: Numeric option text.
 * @param[in] name: Option name used in diagnostics.
 * @return Parsed positive integer.
 */
static int parse_positive(const char *text, const char *name) {
  char *end = NULL;
  errno = 0;
  long value = strtol(text, &end, 10);
  if (errno || end == text || *end || value < 1 || value > INT_MAX)
    bye("invalid %s: %s", name, text);
  return (int)value;
}

/**
 * @brief Run the surf command-line application.
 * @param[in] argc Number of command-line arguments.
 * @param[in] argv Command-line argument vector.
 * @return EXIT_SUCCESS on success, otherwise EXIT_FAILURE.
 */
int main(int argc, char **argv) {
  int threads = 1, iterations = DEFAULT_ITERATIONS;
  const char *model_path = "surf.crf", *prefix = NULL;
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
      threads = parse_positive(optarg, "thread count");
    if (opt == 'i')
      iterations = parse_positive(optarg, "iteration count");
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
  Sequence *sequences = 0;
  read_fasta(argv[a], &count, &sequences);
  if (argc - a == 2) {
    double started = whattime();
    read_gff3(argv[a + 1], sequences, count);
    Model *model = train_pipeline(sequences, count, iterations, threads);
    int status = save_model(model_path, model);
    if (status)
      bye("cannot write CRF model %s: %s", model_path, strerror(status));
    cry("trained in %.1fs; CRF model written to %s", whattime() - started,
        model_path);
    free(model);
  } else {
    Model *model = scalloc(1, sizeof(*model));
    int status = load_model(model_path, model);
    if (status == EINVAL)
      bye("invalid or incompatible CRF model %s; retrain with this version of "
          "surf",
          model_path);
    if (status)
      bye("cannot read CRF model %s: %s", model_path, strerror(status));
    FILE *gff = stdout;
    if (prefix) {
      char path[PATH_BUFFER];
      if (snprintf(path, sizeof(path), "%s.gff3", prefix) >= (int)sizeof(path))
        bye("output prefix too long");
      if (!(gff = fopen(path, "w")))
        bye("cannot create %s", path);
    }
    fputs("##gff-version 3\n", gff);
    unsigned long genes = predict(model, sequences, count, gff, prefix);
    if (fflush(gff) || ferror(gff) || (gff != stdout && fclose(gff)))
      bye("cannot write GFF3");
    cry("predicted %lu genes", genes);
    free(model);
  }
  free_sequences(sequences, count);
  return EXIT_SUCCESS;
}
