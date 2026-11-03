/**
 * @file unet.h
 * @brief 1D UNet runtime.
 */
#ifndef ANNO_ANN_H
#define ANNO_ANN_H

#include "anno.h"
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Return the default UNet configuration.
 * @return Default configuration for the annotation model.
 */
UnetConfig unet_default_config(void);
/**
 * @brief Create an initialized UNet model.
 * @param config Model configuration.
 * @return Newly allocated model, or NULL on invalid configuration/allocation
 * failure.
 */
Unet *unet_create(const UnetConfig *config);
/**
 * @brief Release a model created by unet_create() or unet_load().
 * @param unet Model to release; NULL is accepted.
 */
void unet_destroy(Unet *unet);
/**
 * @brief Evaluate interleaved input features.
 * @param unet Model to evaluate.
 * @param input Interleaved input [position][channel].
 * @param length Number of input positions.
 * @param probability Length x ANNO_CLASSES output probabilities.
 * @return 0 on success, or an errno value otherwise.
 */
int unet_predict(const Unet *unet, const float *input, int64_t length,
                 float *probability);
/**
 * @brief Evaluate channel-major input.
 * @param unet Model to evaluate.
 * @param input Channel-major input [channel][ANNO_UNET_WINDOW].
 * @param length Number of input positions.
 * @param probability Length x ANNO_CLASSES output.
 * @return 0 on success, or an errno value otherwise.
 */
int unet_predict_channels(const Unet *unet, const float *input, int64_t length,
                          float *probability);
/**
 * @brief Train one bounded batch using the masked Soft F1 objective.
 *
 * The objective is the mean per-class Soft F1 loss, evaluated only at
 * positions enabled by mask.
 * @param unet Model to update.
 * @param input Channel-major input batch.
 * @param truth One-hot channel-major target batch.
 * @param mask Position mask batch.
 * @param batch Number of batch elements.
 * @param learning_rate Adam learning rate.
 * @param loss Receives the mean objective value; may be NULL.
 * @return 0 on success, or an errno value otherwise.
 */
int unet_train_batch(Unet *unet, const float *input, const float *truth,
                     const float *mask, size_t batch, float learning_rate,
                     float *loss);
/**
 * @brief Write a UNet model to a binary file.
 * @param path Model file path.
 * @param unet Model to write.
 * @return 0 on success, or an errno value otherwise.
 */
int unet_save(const char *path, const Unet *unet);
/**
 * @brief Read a UNet model from a binary file.
 * @param path Model file path.
 * @param out Receives the loaded model.
 * @return 0 on success, or an errno value otherwise.
 */
int unet_load(const char *path, Unet **out);

#ifdef __cplusplus
}
#endif
#endif
