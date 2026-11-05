/**
 * @file prediction.c
 * @brief Genome inference and optional class bedGraph output.
 */
#include "include/anno_internal.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/** Class names; also used for bedGraph file names. */
static const char *const class_names[ANNO_CLASSES] = {
    "intergenic", "plus_start_single_exon", "plus_end_exon", "plus_intron",
    "minus_start_single_exon", "minus_end_exon", "minus_intron",
    "plus_internal_exon", "minus_internal_exon"};

/* ------------------------------------------------------------------------ */
/* Prediction pipeline                                                      */
/* ------------------------------------------------------------------------ */

/** Shared context of one batch of inference windows. */

typedef struct {
  FILE *files[ANNO_CLASSES];
  int value[ANNO_CLASSES];
  int64_t start[ANNO_CLASSES];
} Bedgraph;

/**
 * @brief Flush one pending bedGraph run.
 * @param[in,out] b BedGraph state containing output files and current runs.
 * @param[in] name Sequence name.
 * @param[in] j Class index.
 * @param[in] position Exclusive genomic end position.
 * @return None.
 */
static void bedgraph_emit(Bedgraph *b, const char *name, int j,
                          int64_t position) {
  if (b->value[j] >= 0 && position > b->start[j] &&
      fprintf(b->files[j], "%s\t%lld\t%lld\t%d.%02d\n", name,
              (long long)b->start[j], (long long)position, b->value[j] / 100,
              b->value[j] % 100) < 0)
    anno_fail("cannot write bedGraph");
}

/**
 * @brief Append a probability row to the class bedGraph runs.
 * @param[in,out] b BedGraph state to update.
 * @param[in] name Sequence name.
 * @param[in] position Genomic position.
 * @param[in] row Class probability row.
 * @return None.
 */
static void bedgraph_add(Bedgraph *b, const char *name, int64_t position,
                         const float *row) {
  if (!b->files[0])
    return;
  for (int j = 0; j < ANNO_CLASSES; j++) {
    int value = (int)lrintf(row[j] * 100.0f);
    if (value == b->value[j])
      continue;
    bedgraph_emit(b, name, j, position);
    b->value[j] = value;
    b->start[j] = position;
  }
}

/**
 * @brief Predict, decode, and write GFF3 and optional bedGraph files.
 * @param[in] model Trained CRF.
 * @param[in] sequences Genome sequences to predict.
 * @param[in] count Number of sequences.
 * @param[in,out] gff Destination GFF3 stream.
 * @param[in] prefix Optional bedGraph output prefix.
 * @return Number of genes written.
 */
unsigned long predict_genome(const CrfModel *model,
                                    const AnnoSequence *sequences, int count,
                                    FILE *gff, const char *prefix) {
  Bedgraph bedgraph = {0};
  for (int j = 0; prefix && j < ANNO_CLASSES; j++) {
    char path[ANNO_PATH_BUFFER];
    if (snprintf(path, sizeof(path), "%s.%s.bedgraph", prefix,
                 class_names[j]) >= (int)sizeof(path))
      anno_fail("output prefix too long");
    if (!(bedgraph.files[j] = fopen(path, "w")))
      anno_fail("cannot create %s", path);
    fprintf(bedgraph.files[j], "track type=bedGraph name=\"%s\"\n",
            class_names[j]);
  }
  unsigned long genes = 0;
  for (int i = 0; i < count; i++) {
    const AnnoSequence *s = &sequences[i];
    Decoder *decoder = decoder_create(s->length);
    for (int j = 0; j < ANNO_CLASSES; j++)
      bedgraph.value[j] = -1;
    int64_t window_count =
        (s->length + CRF_TRAIN_WINDOW - 1) / CRF_TRAIN_WINDOW;
    float *input =
        anno_xcalloc((size_t)ANNO_INPUT_CHANNELS * ANNO_FEATURE_WINDOW, sizeof(float));
    float *power =
        anno_xcalloc((size_t)ANNO_FEATURE_WINDOW * ANNO_CWT_CHANNELS, sizeof(float));

    for (int64_t window = 0; window < window_count; window++) {
      int64_t start = window * CRF_TRAIN_WINDOW;
      int64_t valid = s->length - start < CRF_TRAIN_WINDOW
                          ? s->length - start : CRF_TRAIN_WINDOW;
      memset(input, 0,
             (size_t)ANNO_INPUT_CHANNELS * ANNO_FEATURE_WINDOW * sizeof(float));
      memset(power, 0,
             (size_t)ANNO_FEATURE_WINDOW * ANNO_CWT_CHANNELS * sizeof(float));
      anno_cwt_power(s->bases, s->length, start, start + valid, power);
      anno_features(s, power, start, start, valid, input);
      for (int64_t p = 0; p < valid; p++) {
        float row_features[ANNO_INPUT_CHANNELS], probability[ANNO_CLASSES];
        double scores[ANNO_CLASSES], top = -INFINITY, total = 0.0;
        for (int c = 0; c < ANNO_INPUT_CHANNELS; ++c)
          row_features[c] = input[(size_t)c * ANNO_FEATURE_WINDOW + p];
        crf_emissions(model, row_features, scores);
        for (int c = 0; c < ANNO_CLASSES; ++c)
          if (scores[c] > top)
            top = scores[c];
        for (int c = 0; c < ANNO_CLASSES; ++c) {
          probability[c] = (float)exp(scores[c] - top);
          total += probability[c];
        }
        for (int c = 0; c < ANNO_CLASSES; ++c)
          probability[c] = (float)(probability[c] / total);
        bedgraph_add(&bedgraph, s->name, start + p, probability);
        decoder_step(decoder, scores, model->transition_weights);
      }
    }
    free(input);
    free(power);

    for (int j = 0; j < ANNO_CLASSES && bedgraph.files[0]; j++)
      bedgraph_emit(&bedgraph, s->name, j, s->length);
    uint8_t *path = anno_xcalloc((size_t)s->length, 1);
    if (decoder_finish(decoder, path))
      anno_fail("grammar decoding failed for %s", s->name);
    decoder_destroy(decoder);
    long written = anno_write_gff(gff, s, path, &genes);
    free(path);
    if (written < 0)
      anno_fail("cannot write GFF3");
    anno_log("%s: %ld genes", s->name, written);
  }
  for (int j = 0; j < ANNO_CLASSES; j++)
    if (bedgraph.files[j] && fclose(bedgraph.files[j]))
      anno_fail("cannot close bedGraph");
  return genes;
}
