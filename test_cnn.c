#include "anno_cnn.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int main(void) {
  srand(7);
  CNN *network = cnn_new();
  CNNCache *cache = cnn_cache_new(41);
  CNNGradient *gradient = cnn_gradient_new(network, 0);
  double features[41 * CWT_CHANNELS];
  int targets[41];
  for (int position = 0; position < 41; position++)
    targets[position] = position >= 18 && position < 22 ? position - 18 : -1;
  for (int index = 0; index < 41 * CWT_CHANNELS; index++)
    features[index] = sin(index * 0.3);
  cnn_forward(network, cache, features, 0, 41);
  double before = cnn_loss(network, cache, targets, 0, gradient);
  const int indices[] = {0, 5, 30, 90};
  for (int layer = 0; layer <= CNN_DEPTH; layer++) {
    genann *mapping = layer == CNN_DEPTH ? network->heads[0] : network->layers[layer];
    double *derivative = layer == CNN_DEPTH ? gradient->head : gradient->layers[layer];
    for (int sample = 0; sample < 4; sample++) {
      int weight = indices[sample];
      double original = mapping->weight[weight];
      mapping->weight[weight] = original + 1e-5;
      cnn_forward(network, cache, features, 0, 41);
      double upper = cnn_loss(network, cache, targets, 0, NULL);
      mapping->weight[weight] = original - 1e-5;
      cnn_forward(network, cache, features, 0, 41);
      double lower = cnn_loss(network, cache, targets, 0, NULL);
      mapping->weight[weight] = original;
      assert(fabs((upper - lower) / 2e-5 - derivative[weight]) < 1e-6);
    }
  }
  cnn_update(network, gradient, 0, 0.01);
  cnn_forward(network, cache, features, 0, 41);
  assert(cnn_loss(network, cache, targets, 0, NULL) < before);
  CNNOptimizer *optimizer = cnn_optimizer_new(network, 0);
  cnn_loss(network, cache, targets, 0, gradient);
  CNNGradient *batch = cnn_gradient_new(network, 0);
  cnn_gradient_add(network, batch, gradient, 0, 3.0);
  cnn_gradient_add(network, batch, gradient, 0, 7.0);
  cnn_gradient_scale(network, batch, 0, 0.1);
  for (int layer = 0; layer < CNN_DEPTH; layer++)
    for (int weight = 0; weight < network->layers[layer]->total_weights; weight++)
      assert(fabs(batch->layers[layer][weight] - gradient->layers[layer][weight]) < 1e-12);
  cnn_gradient_free(batch);
  double original_bias = network->heads[0]->weight[0];
  double derivative = gradient->head[0];
  cnn_adam_update(network, gradient, optimizer, 0, 0.0001);
  assert(fabs(network->heads[0]->weight[0] -
              (original_bias - 0.0001 * derivative / (fabs(derivative) + 1e-8))) < 1e-9);
  cnn_optimizer_free(optimizer);
  CNNCache *full = cnn_cache_new(240);
  CNNCache *tile = cnn_cache_new(TILE_SIZE + 2 * CNN_HALO);
  double full_features[240 * CWT_CHANNELS];
  double tile_features[(TILE_SIZE + 2 * CNN_HALO) * CWT_CHANNELS];
  for (int index = 0; index < 240 * CWT_CHANNELS; index++)
    full_features[index] = cos(index * 0.07);
  cnn_forward(network, full, full_features, 0, 240);
  for (int start = 0; start < 240; start += TILE_SIZE) {
    memset(tile_features, 0, sizeof(tile_features));
    for (int offset = 0; offset < tile->count; offset++) {
      int genomic = start - CNN_HALO + offset;
      if (genomic >= 0 && genomic < 240)
        memcpy(tile_features + offset * CWT_CHANNELS,
               full_features + genomic * CWT_CHANNELS, CWT_CHANNELS * sizeof(double));
    }
    cnn_forward(network, tile, tile_features, start - CNN_HALO, 240);
    for (int offset = 0; offset < TILE_SIZE && start + offset < 240; offset++) {
      double full_output[2], tile_output[2];
      cnn_probabilities(network, full, start + offset, 1, full_output);
      cnn_probabilities(network, tile, CNN_HALO + offset, 1, tile_output);
      assert(fabs(full_output[1] - tile_output[1]) < 1e-12);
    }
  }
  cnn_cache_free(full);
  cnn_cache_free(tile);
  FILE *saved = tmpfile();
  assert(saved);
  cnn_write(network, saved);
  rewind(saved);
  CNN *copy = cnn_new();
  cnn_read(copy, saved);
  for (int layer = 0; layer < CNN_DEPTH; layer++)
    assert(memcmp(network->layers[layer]->weight, copy->layers[layer]->weight,
                  network->layers[layer]->total_weights * sizeof(double)) == 0);
  fclose(saved);
  cnn_free(copy);
  cnn_gradient_free(gradient);
  cnn_cache_free(cache);
  cnn_free(network);
  puts("CNN gradient, loss, tile equivalence and serialization tests passed");
  return 0;
}