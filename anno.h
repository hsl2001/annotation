#ifndef ANNO_H
#define ANNO_H

#include <stdint.h>
#include <stdio.h>

#define ANNO_CLASSES 7
#define ANNO_DECODED_STATES 9
#define ANNO_CWT_CHANNELS 3
#define ANNO_CWT_WINDOW 1024
#define ANNO_DEFAULT_ITERATIONS 300

typedef enum {
  ANNO_INTERGENIC = 0,
  ANNO_PLUS_INTRON = 1,
  ANNO_PLUS_5UTR = 2,
  ANNO_PLUS_3UTR = 3,
  ANNO_MINUS_INTRON = 4,
  ANNO_MINUS_5UTR = 5,
  ANNO_MINUS_3UTR = 6
} AnnoClass;

typedef enum {
  ANNO_STATE_INTERGENIC = 0,
  ANNO_STATE_PLUS_5UTR = 1,
  ANNO_STATE_PLUS_EXON = 2,
  ANNO_STATE_PLUS_INTRON = 3,
  ANNO_STATE_PLUS_3UTR = 4,
  ANNO_STATE_MINUS_5UTR = 5,
  ANNO_STATE_MINUS_EXON = 6,
  ANNO_STATE_MINUS_INTRON = 7,
  ANNO_STATE_MINUS_3UTR = 8
} AnnoDecodedState;

typedef struct {
  int input_channels;
  int classes;
  int base_channels;
  int depth;
  int kernel_size;
} AnnoUnetConfig;

typedef struct {
  float feature_mean[ANNO_CWT_CHANNELS + 4];
  float feature_std[ANNO_CWT_CHANNELS + 4];
  uint32_t window_size;
  uint32_t cwt_scales[ANNO_CWT_CHANNELS];
  char class_names[ANNO_CLASSES][24];
} AnnoModelMetadata;

typedef struct AnnoModel AnnoModel;

/* Input: decoded state. Output: matching class, or -1 for an invalid state. */
int anno_class_for_state(uint8_t state);
/* Input: adjacent decoded states. Output: nonzero when the transition is valid. */
int anno_grammar_transition(uint8_t previous, uint8_t current);
/* Input: class probabilities and minimum run lengths. Output: decoded path, or 0. */
int anno_grammar_decode(const float *probability, int64_t length,
                        int min_exon, int min_intron, uint8_t *path);

/* Input: valid network configuration. Output: initialized model, or NULL. */
AnnoModel *anno_model_create(const AnnoUnetConfig *config);
/* Input: model. Output: none; releases all model storage. */
void anno_model_destroy(AnnoModel *model);
/* Input: model and feature matrix. Output: class probabilities, or an error code. */
int anno_model_predict(const AnnoModel *model, const float *input,
                        int64_t length, float *probability);
/* Input: model, labeled feature matrix, and training parameters. Output: 0 or an error code. */
int anno_model_train(AnnoModel *model, const float *input, const uint8_t *label,
                     int64_t length, int iterations, int batch_size,
                     float learning_rate, uint64_t seed);
/* Input: model and destination path. Output: 0 on success, otherwise an error code. */
int anno_model_save(const AnnoModel *model, const char *path);
/* Input: model path. Output: allocated model through model, or an error code. */
int anno_model_load(AnnoModel **model, const char *path);
/* Input: model. Output: read-only metadata, or NULL for a null model. */
const AnnoModelMetadata *anno_model_metadata(const AnnoModel *model);

#endif
