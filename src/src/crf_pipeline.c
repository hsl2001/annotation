/**
 * @file crf_pipeline.c
 * @brief Labeled-window preparation, CRF orchestration, and model I/O.
 */
#include "anno_internal.h"
#include "crf_internal.h"

#include <errno.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/**
 * @brief Extract compact labeled examples for every fixed training window.
 * @param[in] sequences Input sequences with per-base labels.
 * @param[in] sequence_count Number of input sequences.
 * @param[out] example_count Receives the number of extracted windows.
 * @return Allocated examples and their owned feature arrays.
 */
static CrfExample *crf_training_examples(const AnnoSequence *sequences,
                                         int sequence_count,
                                         size_t *example_count) {
  size_t count = 0;
  for (int s = 0; s < sequence_count; ++s) {
    if (!sequences[s].labels || sequences[s].length < 1)
      anno_fail("sequence %s has no training labels", sequences[s].name);
    size_t windows = (size_t)((sequences[s].length - 1) / CRF_TRAIN_WINDOW + 1);
    if (windows > SIZE_MAX - count)
      anno_fail("too many CRF training windows");
    count += windows;
  }
  if (!count)
    anno_fail("no CRF training windows");

  CrfExample *examples = anno_xcalloc(count, sizeof(*examples));
  float *input = anno_xcalloc((size_t)ANNO_INPUT_CHANNELS * ANNO_FEATURE_WINDOW,
                         sizeof(*input));
  size_t at = 0;
  for (int s = 0; s < sequence_count; ++s) {
    const AnnoSequence *sequence = &sequences[s];
    for (int64_t start = 0; start < sequence->length;) {
      int64_t length = sequence->length - start;
      if (length > CRF_TRAIN_WINDOW)
        length = CRF_TRAIN_WINDOW;
      size_t positions = (size_t)length;
      size_t cwt_count = anno_mul_size(positions, ANNO_CWT_CHANNELS);
      float *power = anno_xcalloc(cwt_count, sizeof(*power));
      anno_cwt_power(sequence->bases, sequence->length, start, start + length, power);
      memset(input, 0, (size_t)ANNO_INPUT_CHANNELS * ANNO_FEATURE_WINDOW *
                          sizeof(*input));
      anno_features(sequence, power, start, start, length, input);

      CrfExample *example = &examples[at++];
      example->length = length;
      uint8_t *bases = anno_xcalloc(positions, sizeof(*bases));
      float *cwt = anno_xcalloc(cwt_count, sizeof(*cwt));
      uint8_t *labels = anno_xcalloc(positions, sizeof(*labels));
      for (size_t p = 0; p < positions; ++p) {
        int base = anno_base_index(sequence->bases[start + (int64_t)p]);
        bases[p] = (uint8_t)(base < 0 ? ANNO_NUCLEOTIDES : base);
        labels[p] = sequence->labels[start + (int64_t)p];
        for (int c = 0; c < ANNO_CWT_CHANNELS; ++c)
          cwt[p * ANNO_CWT_CHANNELS + (size_t)c] =
              input[(size_t)(ANNO_NUCLEOTIDES + c) * ANNO_FEATURE_WINDOW + p];
      }
      example->bases = bases;
      example->cwt = cwt;
      example->labels = labels;
      free(power);
      start += length;
    }
  }
  free(input);
  *example_count = count;
  return examples;
}

/**
 * @brief Release a set of training examples and all owned feature buffers.
 * @param[in,out] examples Examples to destroy; NULL is accepted.
 * @param[in] count Number of examples.
 * @return None.
 */
static void crf_training_examples_destroy(CrfExample *examples, size_t count) {
  if (!examples)
    return;
  for (size_t i = 0; i < count; ++i) {
    free((void *)examples[i].bases);
    free((void *)examples[i].cwt);
    free((void *)examples[i].labels);
  }
  free(examples);
}

/**
 * @brief Expose the decoder grammar in the CRF training representation.
 * @param[in] grammar Duration-expanded decoder grammar.
 * @param[out] classes Receives the observed class for each grammar state.
 * @param[out] counts Receives each state's predecessor count.
 * @param[out] predecessors Receives flattened predecessor indices.
 * @return Grammar view referencing the caller-owned output arrays.
 */
CrfGrammar crf_grammar_view(const Grammar *grammar, uint8_t *classes,
                                   uint8_t *counts, int *predecessors) {
  for (int s = 0; s < grammar->count; ++s) {
    uint8_t cls = grammar->cls[s];
    classes[s] = cls == ANNO_PLUS_INTERNAL_EXON || cls == ANNO_MINUS_INTERNAL_EXON
                     ? ANNO_INTERGENIC : cls;
    counts[s] = (uint8_t)grammar->predecessors[s];
    for (int k = 0; k < grammar->predecessors[s]; ++k)
      predecessors[(size_t)s * GRAMMAR_MAX_PREDECESSORS + (size_t)k] =
          grammar->from[s][k];
  }
  return (CrfGrammar){grammar->count, classes, GRAMMAR_MAX_PREDECESSORS,
                      counts, predecessors};
}

/**
 * @brief Train a CRF from labeled sequences using L-BFGS.
 * @param[in] sequences Sequences with per-base training labels.
 * @param[in] count Number of sequences.
 * @param[in] iterations Maximum optimizer iterations.
 * @param[in] threads Number of parallel objective workers.
 * @return Newly allocated trained CRF model.
 */
CrfModel *train_crf(const AnnoSequence *sequences, int count,
                           int iterations, int threads) {
  Grammar grammar;
  grammar_build(&grammar);
  uint8_t classes[GRAMMAR_MAX_STATES], predecessor_count[GRAMMAR_MAX_STATES];
  int predecessors[GRAMMAR_MAX_STATES * GRAMMAR_MAX_PREDECESSORS];
  CrfGrammar view = crf_grammar_view(&grammar, classes, predecessor_count,
                                    predecessors);
  size_t example_count = 0;
  CrfExample *examples = crf_training_examples(sequences, count, &example_count);
  CrfModel *model = anno_xcalloc(1, sizeof(*model));
  double started = anno_now_seconds();
  anno_log("training CRF with L-BFGS: %zu windows, %d iterations", example_count,
      iterations);
  int status = crf_train(model, &view, examples, example_count, iterations,
                         threads);
  crf_training_examples_destroy(examples, example_count);
  if (status != LBFGS_MAX_ITERATIONS && status != LBFGS_CONVERGED &&
      status != LBFGS_STOPPING_CRITERION)
    anno_fail("CRF L-BFGS training failed (status %d)", status);
  anno_log("CRF L-BFGS status %d after %.1fs", status, anno_now_seconds() - started);
  return model;
}

/**
 * @brief Write a versioned CRF model file.
 * @param[in] path Destination model path.
 * @param[in] model Finite model parameters to serialize.
 * @return Zero on success, or an errno-style error code.
 */
int crf_save_model(const char *path, const CrfModel *model) {
  static const char magic[8] = {'A', 'N', 'N', 'O', 'C', 'R', 'F', '1'};
  if (!path || !model || !crf_model_is_finite(model))
    return EINVAL;
  FILE *output = fopen(path, "wb");
  if (!output)
    return errno;
  uint32_t count = CRF_PARAMETER_COUNT;
  int ok = fwrite(magic, sizeof(magic), 1, output) == 1 &&
           fwrite(&count, sizeof(count), 1, output) == 1 &&
           fwrite(model, sizeof(*model), 1, output) == 1;
  int closed = fclose(output) == 0;
  return ok && closed ? 0 : EIO;
}

/**
 * @brief Load and validate a versioned CRF model file.
 * @param[in] path Source model path.
 * @param[out] model Receives the loaded model on success.
 * @return Zero on success, or an errno-style error code.
 */
int crf_load_model(const char *path, CrfModel *model) {
  static const char magic[8] = {'A', 'N', 'N', 'O', 'C', 'R', 'F', '1'};
  if (!path || !model)
    return EINVAL;
  FILE *input = fopen(path, "rb");
  if (!input)
    return errno;
  char found[sizeof(magic)];
  uint32_t count = 0;
  int ok = fread(found, sizeof(found), 1, input) == 1 &&
           fread(&count, sizeof(count), 1, input) == 1 &&
           !memcmp(found, magic, sizeof(magic)) && count == CRF_PARAMETER_COUNT &&
           fread(model, sizeof(*model), 1, input) == 1;
  int extra = fgetc(input);
  if (extra != EOF || ferror(input) || !crf_model_is_finite(model))
    ok = 0;
  int closed = fclose(input) == 0;
  return ok && closed ? 0 : EINVAL;
}

