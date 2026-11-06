/**
 * @file test.c
 * @brief Small, hand-checkable unit tests for every public surf function.
 */

#include "../src/include/crf.h"
#include "../src/include/lbfgs.h"
#include "../src/include/surf.h"
#include "../src/include/utils.h"

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int is_stop(const char *);
int workspace_plan(size_t, size_t, size_t, WorkspacePlan *);

static int failures;
static int checks;

#define CHECK(condition)                                                       \
  do {                                                                         \
    checks++;                                                                  \
    if (!(condition)) {                                                        \
      failures++;                                                              \
      fprintf(stderr, "%s:%d: CHECK failed: %s\n", __FILE__, __LINE__,         \
              #condition);                                                     \
    }                                                                          \
  } while (0)

static FILE *make_temp(char path[], const char *contents) {
  strcpy(path, "/tmp/surf-test-XXXXXX");
  int descriptor = mkstemp(path);
  CHECK(descriptor >= 0);
  if (descriptor < 0)
    return NULL;
  FILE *file = fdopen(descriptor, "w+");
  CHECK(file != NULL);
  if (!file) {
    close(descriptor);
    unlink(path);
    return NULL;
  }
  fputs(contents, file);
  fflush(file);
  rewind(file);
  return file;
}

static int close_enough(double actual, double expected) {
  return fabs(actual - expected) < 1e-10;
}

static double quadratic(const double *x, double *gradient, void *context) {
  (void)context;
  gradient[0] = 2.0 * (x[0] - 3.0);
  return (x[0] - 3.0) * (x[0] - 3.0);
}

static void test_utils(void) {
  CHECK(base_index('A') == 0);
  CHECK(base_index('C') == 1);
  CHECK(base_index('G') == 2);
  CHECK(base_index('T') == 3);
  CHECK(base_index('N') == -1);
  CHECK(complement('A') == 'T');
  CHECK(complement('C') == 'G');
  CHECK(complement('G') == 'C');
  CHECK(complement('T') == 'A');
  CHECK(complement('N') == 'N');
  CHECK(smul(6, 7) == 42);
  int *zeroes = scalloc(3, sizeof(*zeroes));
  CHECK(zeroes[0] == 0 && zeroes[1] == 0 && zeroes[2] == 0);
  zeroes = srealloc(zeroes, 4, sizeof(*zeroes));
  zeroes[3] = 9;
  CHECK(zeroes[3] == 9);
  free(zeroes);
  char *copy = sstrdup("surf");
  CHECK(!strcmp(copy, "surf"));
  free(copy);
  double before = whattime();
  double after = whattime();
  CHECK(isfinite(before) && after >= before);
  cry("utility test message");
}

static void test_fasta_and_gff_input(void) {
  char fasta_path[32], gff_path[32];
  FILE *fasta = make_temp(fasta_path, ">chr1 description\natg aaa taa\n");
  FILE *gff = make_temp(gff_path,
                        "##gff-version 3\n"
                        "chr1\tanno\texon\t1\t9\t.\t+\t.\tID=e1;Parent=tx1\n");
  CHECK(fasta && gff);
  if (!fasta || !gff)
    return;
  fclose(fasta);
  fclose(gff);

  int count = 0;
  Sequence *sequences = NULL;
  read_fasta(fasta_path, &count, &sequences);
  CHECK(count == 1);
  CHECK(!strcmp(sequences[0].name, "chr1"));
  CHECK(!memcmp(sequences[0].bases, "ATGAAATAA", 9));
  CHECK(sequences[0].length == 9);
  read_gff3(gff_path, sequences, count);
  CHECK(sequences[0].labels != NULL);
  for (int i = 0; i < 9; i++)
    CHECK(sequences[0].labels[i] == PLUS_START_SINGLE_EXON);
  free_sequences(sequences, count);
  unlink(fasta_path);
  unlink(gff_path);
}

static void test_features_and_orf(void) {
  const char bases[] = "ACGTN";
  float coefficients[5 * CWT_CHANNELS * COMPLEX_COMPONENTS];
  float power[5 * CWT_CHANNELS];
  cwt(bases, 5, 0, 5, coefficients);
  for (size_t i = 0; i < sizeof(coefficients) / sizeof(*coefficients); i++)
    CHECK(isfinite(coefficients[i]));
  cwt_power(bases, 5, 0, 5, power);
  for (size_t i = 0; i < sizeof(power) / sizeof(*power); i++)
    CHECK(isfinite(power[i]) && power[i] >= 0.0f);

  float features[INPUT_CHANNELS];
  for (int p = 0; p < 5; p++) {
    compute_features(bases[p], power + p * CWT_CHANNELS, features);
    for (int n = 0; n < NUCLEOTIDES; n++)
      CHECK(features[n] == (float)(base_index(bases[p]) == n));
    for (int c = 0; c < CWT_CHANNELS; c++)
      CHECK(isfinite(features[NUCLEOTIDES + c]) &&
            features[NUCLEOTIDES + c] == power[p * CWT_CHANNELS + c]);
  }

  CHECK(is_stop("TAA") && is_stop("TAG") && is_stop("TGA"));
  CHECK(!is_stop("TAC") && !is_stop("ATG"));
  int64_t start = -1, end = -1;
  CHECK(longest_orf("CCCATGAAATAAGGG", 15, &start, &end));
  CHECK(start == 3 && end == 12);
  CHECK(!longest_orf("CCCCCCCC", 8, &start, &end));
}

static void test_cwt_window_independence(void) {
  enum { LENGTH = 2400 };
  char bases[LENGTH];
  float whole[LENGTH * CWT_CHANNELS], split[LENGTH * CWT_CHANNELS];
  for (int i = 0; i < LENGTH; i++)
    bases[i] = "ACGTN"[(i * 13 + i / 7) % 5];
  cwt_power(bases, LENGTH, 0, LENGTH, whole);
  cwt_power(bases, LENGTH, 0, 331, split);
  cwt_power(bases, LENGTH, 331, 1400, split + 331 * CWT_CHANNELS);
  cwt_power(bases, LENGTH, 1400, LENGTH, split + 1400 * CWT_CHANNELS);
  for (int i = 0; i < LENGTH * CWT_CHANNELS; i++)
    CHECK(fabsf(whole[i] - split[i]) < 1e-5f);
}

static void test_intervals_and_gff_output(void) {
  IntervalList list = {0};
  push_interval(&list, 4, 8);
  push_interval(&list, 1, 3);
  CHECK(list.count == 2 && list.capacity >= 2);
  qsort(list.items, (size_t)list.count, sizeof(*list.items), compare_interval);
  CHECK(list.items[0].start == 1 && list.items[1].start == 4);
  CHECK(compare_interval(&list.items[0], &list.items[1]) < 0);
  free(list.items);

  char bases[] = "CCCATGAAATAAGGG";
  uint8_t path[15];
  for (int i = 0; i < 15; i++)
    path[i] = PLUS_START_SINGLE_EXON;
  Sequence sequence = {(char *)"chr1", bases, 15, NULL};
  FILE *output = tmpfile();
  unsigned long gene = 0;
  CHECK(output != NULL);
  CHECK(write_gff(output, &sequence, path, &gene) == 1);
  CHECK(gene == 1);
  rewind(output);
  char text[2048] = {0};
  CHECK(fread(text, 1, sizeof(text) - 1, output) > 0);
  CHECK(strstr(text, "gene1") != NULL);
  CHECK(strstr(text, "\tCDS\t") != NULL);
  fclose(output);
}

static void test_multiexon_gff(void) {
  const char *bases[] = {"NCCATGANNAAACCNNGGGTAACTN",
                         "NAGTTACCCNNGGTTTNNTCATGGN"};
  const int starts[][3] = {{1, 9, 16}, {1, 11, 18}};
  const int ends[][3] = {{7, 14, 24}, {9, 16, 24}};
  const char *expected[] = {
      "chr1\tanno\tgene\t2\t24\t.\t+\t.\tID=gene1\n"
      "chr1\tanno\tmRNA\t2\t24\t.\t+\t.\tID=gene1.t1;Parent=gene1\n"
      "chr1\tanno\texon\t2\t7\t.\t+\t.\tID=gene1.t1.exon1;Parent=gene1.t1\n"
      "chr1\tanno\texon\t10\t14\t.\t+\t.\tID=gene1.t1.exon2;Parent=gene1.t1\n"
      "chr1\tanno\texon\t17\t24\t.\t+\t.\tID=gene1.t1.exon3;Parent=gene1.t1\n"
      "chr1\tanno\tfive_prime_UTR\t2\t3\t.\t+\t.\tID=gene1.t1.utr5p1;Parent="
      "gene1.t1\n"
      "chr1\tanno\tCDS\t4\t7\t.\t+\t0\tID=gene1.t1.cds1;Parent=gene1.t1\n"
      "chr1\tanno\tCDS\t10\t14\t.\t+\t2\tID=gene1.t1.cds2;Parent=gene1.t1\n"
      "chr1\tanno\tCDS\t17\t22\t.\t+\t0\tID=gene1.t1.cds3;Parent=gene1.t1\n"
      "chr1\tanno\tthree_prime_UTR\t23\t24\t.\t+\t.\tID=gene1.t1.utr3p1;Parent="
      "gene1.t1\n",
      "chr1\tanno\tgene\t2\t24\t.\t-\t.\tID=gene1\n"
      "chr1\tanno\tmRNA\t2\t24\t.\t-\t.\tID=gene1.t1;Parent=gene1\n"
      "chr1\tanno\texon\t2\t9\t.\t-\t.\tID=gene1.t1.exon1;Parent=gene1.t1\n"
      "chr1\tanno\texon\t12\t16\t.\t-\t.\tID=gene1.t1.exon2;Parent=gene1.t1\n"
      "chr1\tanno\texon\t19\t24\t.\t-\t.\tID=gene1.t1.exon3;Parent=gene1.t1\n"
      "chr1\tanno\tfive_prime_UTR\t23\t24\t.\t-\t.\tID=gene1.t1.utr5p1;Parent="
      "gene1.t1\n"
      "chr1\tanno\tCDS\t4\t9\t.\t-\t0\tID=gene1.t1.cds1;Parent=gene1.t1\n"
      "chr1\tanno\tCDS\t12\t16\t.\t-\t2\tID=gene1.t1.cds2;Parent=gene1.t1\n"
      "chr1\tanno\tCDS\t19\t22\t.\t-\t0\tID=gene1.t1.cds3;Parent=gene1.t1\n"
      "chr1\tanno\tthree_prime_UTR\t2\t3\t.\t-\t.\tID=gene1.t1.utr3p1;Parent="
      "gene1.t1\n"};
  for (int strand = 0; strand < 2; strand++) {
    uint8_t path[25] = {0};
    for (int e = 0; e < 3; e++) {
      for (int p = starts[strand][e]; p < ends[strand][e]; p++)
        path[p] = strand ? MINUS_START_SINGLE_EXON : PLUS_START_SINGLE_EXON;
      if (e < 2)
        for (int p = ends[strand][e]; p < starts[strand][e + 1]; p++)
          path[p] = strand ? MINUS_INTRON : PLUS_INTRON;
    }
    Sequence sequence = {(char *)"chr1", (char *)bases[strand], 25, NULL};
    FILE *output = tmpfile();
    CHECK(output != NULL);
    if (!output)
      continue;
    unsigned long gene = 0;
    CHECK(write_gff(output, &sequence, path, &gene) == 1 && gene == 1);
    rewind(output);
    char text[2048] = {0};
    CHECK(fread(text, 1, sizeof(text) - 1, output) > 0);
    CHECK(!strcmp(text, expected[strand]));
    fclose(output);

    char file_path[32];
    FILE *file = make_temp(file_path, "");
    CHECK(file != NULL);
    if (file) {
      fclose(file);
      file = fopen(file_path, "r");
      CHECK(file != NULL);
      if (file) {
        gene = 0;
        CHECK(write_gff(file, &sequence, path, &gene) == -1 && gene == 0);
        fclose(file);
      }
      unlink(file_path);
    }
  }
}

static void test_multiexon_gff_input(void) {
  char path[32];
  FILE *gff =
      make_temp(path, "plus\tanno\texon\t2\t7\t.\t+\t.\tParent=tx1\n"
                      "plus\tanno\texon\t10\t14\t.\t+\t.\tParent=tx1\n"
                      "plus\tanno\texon\t17\t24\t.\t+\t.\tParent=tx1\n"
                      "minus\tanno\texon\t2\t9\t.\t-\t.\tParent=tx2\n"
                      "minus\tanno\texon\t12\t16\t.\t-\t.\tParent=tx2\n"
                      "minus\tanno\texon\t19\t24\t.\t-\t.\tParent=tx2\n");
  CHECK(gff != NULL);
  if (!gff)
    return;
  fclose(gff);
  Sequence sequences[] = {
      {(char *)"plus", (char *)"NCCATGANNAAACCNNGGGTAACTN", 25, NULL},
      {(char *)"minus", (char *)"NAGTTACCCNNGGTTTNNTCATGGN", 25, NULL}};
  read_gff3(path, sequences, 2);
  CHECK(sequences[0].labels[1] == PLUS_START_SINGLE_EXON);
  CHECK(sequences[0].labels[7] == PLUS_INTRON);
  CHECK(sequences[0].labels[9] == PLUS_INTERNAL_EXON);
  CHECK(sequences[0].labels[16] == PLUS_END_EXON);
  CHECK(sequences[1].labels[1] == MINUS_END_EXON);
  CHECK(sequences[1].labels[9] == MINUS_INTRON);
  CHECK(sequences[1].labels[11] == MINUS_INTERNAL_EXON);
  CHECK(sequences[1].labels[18] == MINUS_START_SINGLE_EXON);
  free(sequences[0].labels);
  free(sequences[1].labels);
  unlink(path);
}

static void test_decoder_and_grammar(void) {
  Grammar grammar;
  build_grammar(&grammar);
  CHECK(grammar.count ==
        1 + STRANDS * (3 * MIN_EXON_LENGTH + MIN_INTRON_LENGTH));
  CHECK(grammar.cls[0] == INTERGENIC && grammar.backpointers[0] == 5);

  Decoder *decoder = decoder_create(1);
  CHECK(decoder != NULL);
  double emissions[CLASSES] = {0};
  uint8_t path[1] = {255};
  decoder_step(decoder, emissions, NULL);
  CHECK(decoder_finish(decoder, path) == 0 && path[0] == INTERGENIC);
  decoder_destroy(decoder);
  CHECK(decoder_create(0) == NULL);

  decoder = decoder_create(8);
  CHECK(decoder != NULL);
  uint8_t gene_path[8] = {255};
  for (int position = 0; position < 8; position++) {
    for (int class = 0; class < CLASSES; class ++)
      emissions[class] =
          class == (position == 0 || position == 7 ? INTERGENIC
                                                   : PLUS_START_SINGLE_EXON)
              ? 100.0
              : 0.0;
    decoder_step(decoder, emissions, NULL);
  }
  CHECK(decoder_finish(decoder, gene_path) == 0);
  CHECK(gene_path[0] == INTERGENIC && gene_path[7] == INTERGENIC);
  for (int position = 1; position < 7; position++)
    CHECK(gene_path[position] == PLUS_START_SINGLE_EXON);
  decoder_destroy(decoder);
}

static void test_crf_and_lbfgs(void) {
  Model model = {0};
  float features[INPUT_CHANNELS] = {0};
  double scores[CLASSES];
  model.emission_weights[2 * (INPUT_CHANNELS + 1) + 0] = 2.0;
  model.emission_weights[2 * (INPUT_CHANNELS + 1) + INPUT_CHANNELS] = 1.0;
  features[0] = 3.0f;
  compute_emissions(&model, features, scores);
  CHECK(close_enough(scores[2], 7.0));
  CHECK(close_enough(scores[0], 0.0));
  CHECK(is_model_finite(&model));
  model.transition_weights[0] = -2.0;
  CHECK(is_model_finite(&model));
  model.transition_weights[0] = 0.0;

  Grammar grammar;
  build_grammar(&grammar);
  CHECK(grammar.cls[0] == INTERGENIC && grammar.backpointers[0] == 5);
  WorkspacePlan plan;
  CHECK(workspace_plan(5, 2, 2, &plan));
  CHECK(plan.block_length == 2 && plan.blocks == 3);
  CHECK(!workspace_plan(0, 2, 0, &plan));

  uint8_t bases[1] = {0}, labels[1] = {INTERGENIC};
  float cwt[3] = {0};
  Sample sample = {1, bases, cwt, labels};
  ModelObjectiveContext context = {&grammar, &sample, 1, 1, 1, 0};
  double parameters[PARAMETERS] = {0};
  double gradient[PARAMETERS];
  double objective = model_objective(parameters, gradient, &context);
  CHECK(isfinite(objective) && close_enough(objective, 0.0));
  CHECK(fabs(gradient[0]) < 1e-12);
  model.emission_weights[INPUT_CHANNELS] = NAN;
  CHECK(!is_model_finite(&model));
  model.emission_weights[INPUT_CHANNELS] = 0.0;
  CHECK(train(&model, &grammar, &sample, 1, 1, 1) >= 0);
  Grammar invalid = grammar;
  invalid.count = GRAMMAR_MAX_STATES + 1;
  CHECK(train(&model, &invalid, &sample, 1, 1, 1) == LBFGS_INVALID_ARGUMENT);
  invalid = grammar;
  invalid.backpointers[0] = GRAMMAR_MAX_BACKPOINTERS + 1;
  CHECK(train(&model, &invalid, &sample, 1, 1, 1) == LBFGS_INVALID_ARGUMENT);
  invalid = grammar;
  invalid.from[0][0] = invalid.count;
  CHECK(train(&model, &invalid, &sample, 1, 1, 1) == LBFGS_INVALID_ARGUMENT);

  double x[1] = {0};
  int completed = -1;
  LbfgsStatus status =
      lbfgs_minimize(1, x, quadratic, NULL, NULL, 20, &completed);
  CHECK(status == LBFGS_CONVERGED || status == LBFGS_STOPPING_CRITERION ||
        status == LBFGS_MAX_ITERATIONS);
  CHECK(fabs(x[0] - 3.0) < 1e-4);
  CHECK(completed >= 0);
  CHECK(lbfgs_minimize(0, x, quadratic, NULL, NULL, 1, &completed) ==
        LBFGS_INVALID_ARGUMENT);

  char model_path[32];
  FILE *model_file = make_temp(model_path, "");
  CHECK(model_file != NULL);
  if (model_file)
    fclose(model_file);
  CHECK(save_model(model_path, &model) == 0);
  Model loaded = {0};
  CHECK(load_model(model_path, &loaded) == 0);
  CHECK(!memcmp(&model, &loaded, sizeof(model)));
  model_file = fopen(model_path, "r+b");
  CHECK(model_file != NULL);
  if (model_file) {
    CHECK(fwrite("SURF", 1, 4, model_file) == 4);
    CHECK(fclose(model_file) == 0);
    CHECK(load_model(model_path, &loaded) == EINVAL);
    CHECK(!memcmp(&model, &loaded, sizeof(model)));
  }
  unlink(model_path);
}

static void test_crf_gradient(void) {
  Grammar grammar = {0};
  grammar.count = 2;
  grammar.cls[1] = PLUS_START_SINGLE_EXON;
  for (int state = 0; state < grammar.count; state++) {
    grammar.backpointers[state] = 2;
    grammar.from[state][0] = 0;
    grammar.from[state][1] = 1;
  }
  uint8_t bases[] = {0, 1, 2, 3, 0};
  uint8_t labels[] = {0, 1, 1, 0, 0};
  float cwt[5 * CWT_CHANNELS];
  for (int i = 0; i < 5 * CWT_CHANNELS; i++)
    cwt[i] = (float)(i - 7) / 10;
  Sample sample = {5, bases, cwt, labels};
  ModelObjectiveContext context = {&grammar, &sample, 1, 5, 1, 0};
  double parameters[PARAMETERS], gradient[PARAMETERS], scratch[PARAMETERS];
  for (int i = 0; i < PARAMETERS; i++)
    parameters[i] = (double)(i % 7 - 3) / 10;
  CHECK(isfinite(model_objective(parameters, gradient, &context)));
  for (int i = 0; i < PARAMETERS; i++) {
    double original = parameters[i];
    parameters[i] = original + 1e-6;
    double plus = model_objective(parameters, scratch, &context);
    parameters[i] = original - 1e-6;
    double minus = model_objective(parameters, scratch, &context);
    parameters[i] = original;
    CHECK(isfinite(plus) && isfinite(minus) &&
          fabs(gradient[i] - (plus - minus) / 2e-6) < 1e-7);
  }
}

static void test_prediction(void) {
  Model model = {0};
  char bases[] = "A";
  Sequence sequence = {(char *)"tiny", bases, 1, NULL};
  FILE *output = tmpfile();
  CHECK(output != NULL);
  if (!output)
    return;
  CHECK(predict(&model, &sequence, 1, output, NULL) == 0);
  fclose(output);
}

int main(void) {
  test_utils();
  test_fasta_and_gff_input();
  test_features_and_orf();
  test_cwt_window_independence();
  test_intervals_and_gff_output();
  test_multiexon_gff();
  test_multiexon_gff_input();
  test_decoder_and_grammar();
  test_crf_and_lbfgs();
  test_crf_gradient();
  test_prediction();
  fprintf(stderr, "%d checks, %d failures\n", checks, failures);
  return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
