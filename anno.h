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

/**
 * @brief Executes anno_class_for_state.
 * @param state Input value.
 * @return Function result.
 */
int anno_class_for_state(uint8_t state);
/**
 * @brief Executes anno_grammar_transition.
 * @param previous Input value.
 * @param current Input value.
 * @return Function result.
 */
int anno_grammar_transition(uint8_t previous, uint8_t current);
/**
 * @brief Executes anno_grammar_decode.
 * @param probability Input value.
 * @param length Input value.
 * @param min_exon Input value.
 * @param min_intron Input value.
 * @param path Input value.
 * @return Function result.
 */
int anno_grammar_decode(const float *probability, int64_t length,
                        int min_exon, int min_intron, uint8_t *path);

/**
 * @brief Executes anno_model_create.
 * @param config Input value.
 * @return Function result.
 */
AnnoModel *anno_model_create(const AnnoUnetConfig *config);
/**
 * @brief Executes anno_model_destroy.
 * @param model Input value.
 * @return No value.
 */
void anno_model_destroy(AnnoModel *model);
/**
 * @brief Executes anno_model_predict.
 * @param model Input value.
 * @param input Input value.
 * @param length Input value.
 * @param probability Input value.
 * @return Function result.
 */
int anno_model_predict(const AnnoModel *model, const float *input,
                        int64_t length, float *probability);
/**
 * @brief Executes anno_model_train.
 * @param model Input value.
 * @param input Input value.
 * @param label Input value.
 * @param length Input value.
 * @param iterations Input value.
 * @param batch_size Input value.
 * @param learning_rate Input value.
 * @param seed Input value.
 * @return Function result.
 */
int anno_model_train(AnnoModel *model, const float *input, const uint8_t *label,
                     int64_t length, int iterations, int batch_size,
                     float learning_rate, uint64_t seed);
/**
 * @brief Executes anno_model_save.
 * @param model Input value.
 * @param path Input value.
 * @return Function result.
 */
int anno_model_save(const AnnoModel *model, const char *path);
/**
 * @brief Executes anno_model_load.
 * @param model Input value.
 * @param path Input value.
 * @return Function result.
 */
int anno_model_load(AnnoModel **model, const char *path);
/**
 * @brief Executes anno_model_metadata.
 * @param model Input value.
 * @return Function result.
 */
const AnnoModelMetadata *anno_model_metadata(const AnnoModel *model);

#endif
