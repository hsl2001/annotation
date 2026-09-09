#include "anno_cnn.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

static const int dilations[CNN_DEPTH] = {1, 2, 4, 8, 16, 1, 2, 4};

static size_t cell(int position, int scale, int channel) {
  return ((size_t)position * WAVE_COUNT + scale) * CNN_WIDTH + channel;
}

static genann *linear_new(int inputs, int outputs) {
  genann *layer = genann_init(inputs, 0, 0, outputs);
  if (!layer)
    anno_fail("Cannot allocate genann layer");
  layer->activation_output = genann_act_linear;
  double limit = sqrt(6.0 / (inputs + outputs));
  for (int output = 0; output < outputs; output++) {
    int base = output * (inputs + 1);
    layer->weight[base] = 0.0;
    for (int input = 0; input < inputs; input++)
      layer->weight[base + input + 1] =
          (2.0 * rand() / RAND_MAX - 1.0) * limit;
  }
  return layer;
}

CNN *cnn_new(void) {
  CNN *network = anno_alloc(1, sizeof(*network));
  for (int layer = 0; layer < CNN_DEPTH; layer++)
    network->layers[layer] = linear_new(9 * (layer ? CNN_WIDTH : 2), CNN_WIDTH);
  network->heads[0] = linear_new(WAVE_COUNT * CNN_WIDTH, 4);
  network->heads[1] = linear_new(WAVE_COUNT * CNN_WIDTH, 2);
  return network;
}

void cnn_free(CNN *network) {
  for (int layer = 0; layer < CNN_DEPTH; layer++)
    genann_free(network->layers[layer]);
  genann_free(network->heads[0]);
  genann_free(network->heads[1]);
  free(network);
}

CNNCache *cnn_cache_new(int count) {
  CNNCache *cache = anno_alloc(1, sizeof(*cache));
  cache->count = count;
  for (int layer = 0; layer <= CNN_DEPTH; layer++)
    cache->values[layer] = anno_alloc((size_t)count * WAVE_COUNT * CNN_WIDTH,
                                      sizeof(double));
  return cache;
}

void cnn_cache_free(CNNCache *cache) {
  for (int layer = 0; layer <= CNN_DEPTH; layer++)
    free(cache->values[layer]);
  free(cache);
}

static int valid_position(const CNNCache *cache, int position) {
  long genomic = (long)cache->origin + position;
  return position >= 0 && position < cache->count && genomic >= 0 &&
         genomic < cache->sequence_length;
}

static void patch_get(const CNNCache *cache, int layer, int position, int scale,
                      double *patch) {
  int channels = layer ? CNN_WIDTH : 2;
  int input = 0;
  for (int shift = -1; shift <= 1; shift++) {
    int neighbor = position + shift * dilations[layer];
    for (int scale_shift = -1; scale_shift <= 1; scale_shift++) {
      int neighbor_scale = scale + scale_shift;
      for (int channel = 0; channel < channels; channel++)
        patch[input++] = valid_position(cache, neighbor) && neighbor_scale >= 0 &&
                                 neighbor_scale < WAVE_COUNT
                             ? cache->values[layer][cell(neighbor, neighbor_scale, channel)]
                             : 0.0;
    }
  }
}

void cnn_forward(CNN *network, CNNCache *cache, const double *features,
                 int origin, int sequence_length) {
  cache->origin = origin;
  cache->sequence_length = sequence_length;
  size_t cells = (size_t)cache->count * WAVE_COUNT * CNN_WIDTH;
  memset(cache->values[0], 0, cells * sizeof(double));
  for (int position = 0; position < cache->count; position++)
    for (int scale = 0; scale < WAVE_COUNT; scale++)
      for (int part = 0; part < 2; part++)
        cache->values[0][cell(position, scale, part)] =
            features[(size_t)position * CWT_CHANNELS + 2 * scale + part];
  for (int layer = 0; layer < CNN_DEPTH; layer++) {
    memset(cache->values[layer + 1], 0, cells * sizeof(double));
    for (int position = 0; position < cache->count; position++) {
      if (!valid_position(cache, position))
        continue;
      for (int scale = 0; scale < WAVE_COUNT; scale++) {
        double patch[9 * CNN_WIDTH];
        patch_get(cache, layer, position, scale, patch);
        const double *output = genann_run(network->layers[layer], patch);
        for (int channel = 0; channel < CNN_WIDTH; channel++) {
          size_t index = cell(position, scale, channel);
          double residual = layer ? cache->values[layer][index] : 0.0;
          cache->values[layer + 1][index] = tanh(output[channel] + residual);
        }
      }
    }
  }
}

void cnn_probabilities(CNN *network, CNNCache *cache, int position, int task,
                       double *probabilities) {
  genann *head = network->heads[task];
  const double *logits = genann_run(head, cache->values[CNN_DEPTH] + cell(position, 0, 0));
  double maximum = logits[0];
  for (int output = 1; output < head->outputs; output++)
    maximum = fmax(maximum, logits[output]);
  double sum = 0.0;
  for (int output = 0; output < head->outputs; output++) {
    probabilities[output] = exp(logits[output] - maximum);
    sum += probabilities[output];
  }
  for (int output = 0; output < head->outputs; output++)
    probabilities[output] /= sum;
}

CNNGradient *cnn_gradient_new(CNN *network, int task) {
  CNNGradient *gradient = anno_alloc(1, sizeof(*gradient));
  for (int layer = 0; layer < CNN_DEPTH; layer++)
    gradient->layers[layer] = anno_alloc(network->layers[layer]->total_weights,
                                          sizeof(double));
  gradient->head = anno_alloc(network->heads[task]->total_weights, sizeof(double));
  return gradient;
}

void cnn_gradient_free(CNNGradient *gradient) {
  for (int layer = 0; layer < CNN_DEPTH; layer++)
    free(gradient->layers[layer]);
  free(gradient->head);
  free(gradient);
}

void cnn_gradient_add(CNN *network, CNNGradient *destination,
                       const CNNGradient *source, int task, double factor) {
  for (int layer = 0; layer <= CNN_DEPTH; layer++) {
    genann *mapping = layer == CNN_DEPTH ? network->heads[task] : network->layers[layer];
    double *output = layer == CNN_DEPTH ? destination->head : destination->layers[layer];
    const double *input = layer == CNN_DEPTH ? source->head : source->layers[layer];
    for (int weight = 0; weight < mapping->total_weights; weight++)
      output[weight] += factor * input[weight];
  }
}

void cnn_gradient_scale(CNN *network, CNNGradient *gradient, int task, double factor) {
  for (int layer = 0; layer <= CNN_DEPTH; layer++) {
    genann *mapping = layer == CNN_DEPTH ? network->heads[task] : network->layers[layer];
    double *values = layer == CNN_DEPTH ? gradient->head : gradient->layers[layer];
    for (int weight = 0; weight < mapping->total_weights; weight++) values[weight] *= factor;
  }
}

static void linear_backward(const genann *layer, const double *input,
                            const double *delta, double *weights, double *dx) {
  for (int output = 0; output < layer->outputs; output++) {
    int base = output * (layer->inputs + 1);
    weights[base] -= delta[output];
    for (int channel = 0; channel < layer->inputs; channel++) {
      weights[base + channel + 1] += delta[output] * input[channel];
      dx[channel] += delta[output] * layer->weight[base + channel + 1];
    }
  }
}

double cnn_loss(CNN *network, CNNCache *cache, const int *targets, int task,
                CNNGradient *gradient) {
  int examples = 0;
  for (int position = 0; position < cache->count; position++)
    if (targets[position] >= 0 && valid_position(cache, position))
      examples++;
  if (!examples)
    return 0.0;
  size_t cells = (size_t)cache->count * WAVE_COUNT * CNN_WIDTH;
  double *current = gradient ? anno_alloc(cells, sizeof(double)) : NULL;
  double *previous = gradient ? anno_alloc(cells, sizeof(double)) : NULL;
  genann *head = network->heads[task];
  if (gradient) {
    memset(gradient->head, 0, head->total_weights * sizeof(double));
    for (int layer = 0; layer < CNN_DEPTH; layer++)
      memset(gradient->layers[layer], 0,
             network->layers[layer]->total_weights * sizeof(double));
  }
  double loss = 0.0;
  for (int position = 0; position < cache->count; position++) {
    int target = targets[position];
    if (target < 0 || !valid_position(cache, position))
      continue;
    if (target >= head->outputs)
      anno_fail("Invalid CNN target");
    double probabilities[4];
    cnn_probabilities(network, cache, position, task, probabilities);
    loss -= log(fmax(probabilities[target], 1e-300)) / examples;
    if (gradient) {
      double delta[4];
      for (int output = 0; output < head->outputs; output++)
        delta[output] = (probabilities[output] - (output == target)) / examples;
      linear_backward(head, cache->values[CNN_DEPTH] + cell(position, 0, 0),
                      delta, gradient->head, current + cell(position, 0, 0));
    }
  }
  if (!gradient)
    return loss;
  for (int layer = CNN_DEPTH - 1; layer >= 0; layer--) {
    memset(previous, 0, cells * sizeof(double));
    int channels = layer ? CNN_WIDTH : 2;
    for (int position = 0; position < cache->count; position++) {
      if (!valid_position(cache, position))
        continue;
      for (int scale = 0; scale < WAVE_COUNT; scale++) {
        double patch[9 * CNN_WIDTH], dx[9 * CNN_WIDTH] = {0};
        double delta[CNN_WIDTH];
        patch_get(cache, layer, position, scale, patch);
        for (int channel = 0; channel < CNN_WIDTH; channel++) {
          size_t index = cell(position, scale, channel);
          double value = cache->values[layer + 1][index];
          delta[channel] = current[index] * (1.0 - value * value);
          if (layer)
            previous[index] += delta[channel];
        }
        linear_backward(network->layers[layer], patch, delta,
                        gradient->layers[layer], dx);
        int input = 0;
        for (int shift = -1; shift <= 1; shift++) {
          int neighbor = position + shift * dilations[layer];
          for (int scale_shift = -1; scale_shift <= 1; scale_shift++) {
            int neighbor_scale = scale + scale_shift;
            for (int channel = 0; channel < channels; channel++) {
              if (valid_position(cache, neighbor) && neighbor_scale >= 0 &&
                  neighbor_scale < WAVE_COUNT)
                previous[cell(neighbor, neighbor_scale, channel)] += dx[input];
              input++;
            }
          }
        }
      }
    }
    double *swap = current;
    current = previous;
    previous = swap;
  }
  free(current);
  free(previous);
  return loss;
}

void cnn_update(CNN *network, CNNGradient *gradient, int task, double rate) {
  double norm = 0.0;
  for (int layer = 0; layer <= CNN_DEPTH; layer++) {
    genann *mapping = layer == CNN_DEPTH ? network->heads[task] : network->layers[layer];
    double *derivative = layer == CNN_DEPTH ? gradient->head : gradient->layers[layer];
    for (int weight = 0; weight < mapping->total_weights; weight++) {
      if (!isfinite(derivative[weight]))
        anno_fail("Non-finite CNN gradient; reduce learning rate");
      norm += derivative[weight] * derivative[weight];
    }
  }
  double step = rate / fmax(1.0, sqrt(norm) / 5.0);
  for (int layer = 0; layer <= CNN_DEPTH; layer++) {
    genann *mapping = layer == CNN_DEPTH ? network->heads[task] : network->layers[layer];
    double *derivative = layer == CNN_DEPTH ? gradient->head : gradient->layers[layer];
    for (int weight = 0; weight < mapping->total_weights; weight++)
      mapping->weight[weight] -= step * derivative[weight];
  }
}

CNNOptimizer *cnn_optimizer_new(CNN *network, int task) {
  CNNOptimizer *optimizer = anno_alloc(1, sizeof(*optimizer));
  optimizer->first = cnn_gradient_new(network, task);
  optimizer->second = cnn_gradient_new(network, task);
  return optimizer;
}

void cnn_optimizer_free(CNNOptimizer *optimizer) {
  cnn_gradient_free(optimizer->first);
  cnn_gradient_free(optimizer->second);
  free(optimizer);
}

void cnn_adam_update(CNN *network, CNNGradient *gradient, CNNOptimizer *optimizer,
                     int task, double rate) {
  double norm = 0.0;
  for (int layer = 0; layer <= CNN_DEPTH; layer++) {
    genann *mapping = layer == CNN_DEPTH ? network->heads[task] : network->layers[layer];
    double *derivative = layer == CNN_DEPTH ? gradient->head : gradient->layers[layer];
    for (int weight = 0; weight < mapping->total_weights; weight++) {
      if (!isfinite(derivative[weight])) anno_fail("Non-finite Adam gradient");
      norm += derivative[weight] * derivative[weight];
    }
  }
  double clipping = 1.0 / fmax(1.0, sqrt(norm) / 5.0);
  optimizer->step++;
  double first_correction = 1.0 - pow(0.9, optimizer->step);
  double second_correction = 1.0 - pow(0.999, optimizer->step);
  for (int layer = 0; layer <= CNN_DEPTH; layer++) {
    genann *mapping = layer == CNN_DEPTH ? network->heads[task] : network->layers[layer];
    double *derivative = layer == CNN_DEPTH ? gradient->head : gradient->layers[layer];
    double *first = layer == CNN_DEPTH ? optimizer->first->head : optimizer->first->layers[layer];
    double *second = layer == CNN_DEPTH ? optimizer->second->head : optimizer->second->layers[layer];
    for (int weight = 0; weight < mapping->total_weights; weight++) {
      double value = derivative[weight] * clipping;
      first[weight] = 0.9 * first[weight] + 0.1 * value;
      second[weight] = 0.999 * second[weight] + 0.001 * value * value;
      mapping->weight[weight] -= rate * (first[weight] / first_correction) /
          (sqrt(second[weight] / second_correction) + 1e-8);
    }
  }
}

void cnn_write(CNN *network, FILE *stream) {
  for (int layer = 0; layer < CNN_DEPTH + 2; layer++) {
    genann *mapping = layer < CNN_DEPTH ? network->layers[layer] : network->heads[layer - CNN_DEPTH];
    fprintf(stream, "%d %d\n", mapping->inputs, mapping->outputs);
    for (int weight = 0; weight < mapping->total_weights; weight++)
      fprintf(stream, "%.17g%c", mapping->weight[weight],
              weight + 1 == mapping->total_weights ? '\n' : ' ');
  }
}

void cnn_read(CNN *network, FILE *stream) {
  for (int layer = 0; layer < CNN_DEPTH + 2; layer++) {
    genann *mapping = layer < CNN_DEPTH ? network->layers[layer] : network->heads[layer - CNN_DEPTH];
    int inputs, outputs;
    if (fscanf(stream, "%d %d", &inputs, &outputs) != 2 ||
        inputs != mapping->inputs || outputs != mapping->outputs)
      anno_fail("Incompatible or corrupt CNN dimensions");
    for (int weight = 0; weight < mapping->total_weights; weight++)
      if (fscanf(stream, "%lf", &mapping->weight[weight]) != 1 ||
          !isfinite(mapping->weight[weight]))
        anno_fail("Truncated or non-finite CNN weights");
  }
}