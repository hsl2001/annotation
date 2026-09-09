#ifndef ANNO_CNN_H
#define ANNO_CNN_H

#include "anno_signal.h"
#include "genann.h"

#define CNN_DEPTH 8
#define CNN_WIDTH 8
#define CNN_HALO 38
#define TILE_SIZE 128

typedef struct {
  genann *layers[CNN_DEPTH];
  genann *heads[2];
} CNN;

typedef struct {
  int count, origin, sequence_length;
  double *values[CNN_DEPTH + 1];
} CNNCache;

typedef struct {
  double *layers[CNN_DEPTH];
  double *head;
} CNNGradient;

typedef struct {
  CNNGradient *first, *second;
  unsigned long step;
} CNNOptimizer;

CNN *cnn_new(void);
void cnn_free(CNN *network);
CNNCache *cnn_cache_new(int count);
void cnn_cache_free(CNNCache *cache);
void cnn_forward(CNN *network, CNNCache *cache, const double *features,
                 int origin, int sequence_length);
void cnn_probabilities(CNN *network, CNNCache *cache, int position, int task,
                       double *probabilities);
CNNGradient *cnn_gradient_new(CNN *network, int task);
void cnn_gradient_free(CNNGradient *gradient);
void cnn_gradient_add(CNN *network, CNNGradient *destination,
                       const CNNGradient *source, int task, double factor);
void cnn_gradient_scale(CNN *network, CNNGradient *gradient, int task, double factor);
double cnn_loss(CNN *network, CNNCache *cache, const int *targets, int task,
                CNNGradient *gradient);
void cnn_update(CNN *network, CNNGradient *gradient, int task, double rate);
CNNOptimizer *cnn_optimizer_new(CNN *network, int task);
void cnn_optimizer_free(CNNOptimizer *optimizer);
void cnn_adam_update(CNN *network, CNNGradient *gradient, CNNOptimizer *optimizer,
                     int task, double rate);
void cnn_write(CNN *network, FILE *stream);
void cnn_read(CNN *network, FILE *stream);

#endif