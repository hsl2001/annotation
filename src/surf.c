#include "include/model.h"
#include "include/utils.h"

#include <errno.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#ifdef _OPENMP
#include <omp.h>
#endif

static void usage(void) {
  fprintf(
      stderr,
      "Usage: surf [options] genome.fasta[.gz] [annotation.gff3[.gz]]\n"
      "Options:\n"
      "  -p N       threads [1]\n"
      "  -i N       training iterations [%d]\n"
      "  -m FILE    point model file [surf.model]\n"
      "  -o PREFIX  write GFF3 and strand-specific point bedGraphs\n"
      "  -h         this help\n",
      DEFAULT_ITERATIONS);
}

static int positive(const char *text) {
  char *end;
  errno = 0;
  long value = strtol(text, &end, 10);
  if (errno || end == text || *end || value < 1 || value > INT_MAX)
    bye("invalid positive integer: %s", text);
  return (int)value;
}

int main(int argc, char **argv) {
  int threads = 1, iterations = DEFAULT_ITERATIONS, opt;
  const char *model_path = "surf.model", *prefix = NULL;
  opterr = 0;
  while ((opt = getopt(argc, argv, "hp:i:m:o:")) != -1) {
    if (opt == 'h' || opt == '?') {
      usage();
      return opt == 'h' ? EXIT_SUCCESS : EXIT_FAILURE;
    }
    if (opt == 'p')
      threads = positive(optarg);
    if (opt == 'i')
      iterations = positive(optarg);
    if (opt == 'm')
      model_path = optarg;
    if (opt == 'o')
      prefix = optarg;
  }
  if (argc - optind < 1 || argc - optind > 2) {
    usage();
    return EXIT_FAILURE;
  }
#ifdef _OPENMP
  omp_set_num_threads(threads);
#endif
  int count;
  Sequence *sequences;
  read_fasta(argv[optind], &count, &sequences);
  Model *model;
  if (argc - optind == 2) {
    double started = whattime();
    read_gff3(argv[optind + 1], sequences, count);
    model = train_model(sequences, count, iterations, threads);
    int status = save_model(model_path, model);
    if (status)
      bye("cannot write model %s: %s", model_path, strerror(status));
    cry("trained in %.2fs; model written to %s", whattime() - started,
        model_path);
  } else {
    model = scalloc(1, sizeof(*model));
    int status = load_model(model_path, model);
    if (status == EINVAL)
      bye("incompatible point model %s; retrain with this version", model_path);
    if (status)
      bye("cannot read model %s: %s", model_path, strerror(status));
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
  }
  free(model);
  free_sequences(sequences, count);
  return EXIT_SUCCESS;
}
