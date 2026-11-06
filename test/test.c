/**
 * @file test.c
 * @brief Small, hand-checkable unit tests for every public surf function.
 */

#include "../src/include/crf.h"
#include "../src/include/lbfgs.h"
#include "../src/include/surf.h"
#include "../src/include/utils.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

typedef struct {
  int64_t start;
  int64_t end;
} Interval;

typedef struct {
  Interval *items;
  int count;
  int capacity;
} IntervalList;

void push_interval(IntervalList *, int64_t, int64_t);
int compare_interval(const void *, const void *);
int is_stop(const char *);
LongGrammar long2short(const ShortGrammar *, uint8_t *, uint8_t *, int *);
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
  FILE *gff = make_temp(
      gff_path,
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
  Sequence sequence = {(char *)"tiny", (char *)bases, 5, NULL};
  cwt(bases, 5, 0, 5, coefficients);
  for (size_t i = 0; i < sizeof(coefficients) / sizeof(*coefficients); i++)
    CHECK(isfinite(coefficients[i]));
  cwt_power(bases, 5, 0, 5, power);
  for (size_t i = 0; i < sizeof(power) / sizeof(*power); i++)
    CHECK(isfinite(power[i]) && power[i] >= 0.0f);

  float features[INPUT_CHANNELS * 5];
  memset(features, 0, sizeof(features));
  compute_features(&sequence, power, 0, 0, 5, features);
  CHECK(features[0] == 1.0f && features[1] == 0.0f);
  CHECK(features[3 * 5 + 4] == 0.0f);
  for (int i = NUCLEOTIDES * 5; i < INPUT_CHANNELS * 5; i++)
    CHECK(isfinite(features[i]));

  CHECK(is_stop("TAA") && is_stop("TAG") && is_stop("TGA"));
  CHECK(!is_stop("TAC") && !is_stop("ATG"));
  int64_t start = -1, end = -1;
  CHECK(longest_orf("CCCATGAAATAAGGG", 15, &start, &end));
  CHECK(start == 3 && end == 12);
  CHECK(!longest_orf("CCCCCCCC", 8, &start, &end));
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

static void test_decoder_and_grammar(void) {
  ShortGrammar grammar;
  build_grammar(&grammar);
  CHECK(grammar.count == 1 + STRANDS *
                                     (3 * MIN_EXON_LENGTH + MIN_INTRON_LENGTH));
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
    for (int class = 0; class < CLASSES; class++)
      emissions[class] = class == (position == 0 || position == 7
                                       ? INTERGENIC
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

  ShortGrammar short_grammar;
  build_grammar(&short_grammar);
  uint8_t classes[SHORTGRAMMER_MAX_STATES];
  uint8_t counts[SHORTGRAMMER_MAX_STATES];
  int backpointers[SHORTGRAMMER_MAX_STATES * SHORTGRAMMAR_MAX_BACKPOINTERS];
  LongGrammar grammar = long2short(&short_grammar, classes, counts, backpointers);
  CHECK(grammar.states == short_grammar.count);
  CHECK(grammar.state_class[0] == INTERGENIC && grammar.backpointer_count[0] == 5);
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
  unlink(model_path);
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
  test_intervals_and_gff_output();
  test_decoder_and_grammar();
  test_crf_and_lbfgs();
  test_prediction();
  fprintf(stderr, "%d checks, %d failures\n", checks, failures);
  return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
