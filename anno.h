#ifndef ANNO_H
#define ANNO_H

#include <stddef.h>
#include <stdint.h>

/* CWT */
#ifndef CWT_SCALES
#define CWT_SCALES 4, 5, 6, 7, 8, 9, 16, 32, 64, 128
#endif
#define CWT_WINDOW_SIZE 1024
#define MORLET_OMEGA0 6.0
#define MORLET_CUTOFF 8.0

/* CRF features and decoding */
#define SPLICE_RADIUS 12
#define POWER_SLOPE_EPSILON 1e-12
#define MIN_CDS 30

/* Training */
#define TRAINING_WINDOW 65536
#define DEFAULT_EPOCHS 72
#define ADAM_LEARNING_RATE 0.03
#define ADAM_BETA1 0.9
#define ADAM_BETA2 0.999
#define ADAM_EPSILON 1e-8

/* Model file */
#define DEFAULT_MODEL "anno.model"
#define MODEL_MAGIC "ANNOCRFE"

#define CWT_COUNT_IMPL(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, \
                       _17, _18, _19, _20, _21, _22, _23, _24, _25, _26, _27, _28, _29, _30, \
                       _31, _32, count, ...) count
#define CWT_COUNT(...) CWT_COUNT_IMPL(__VA_ARGS__, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, \
                                      21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, \
                                      6, 5, 4, 3, 2, 1)
#define WAVE_COUNT CWT_COUNT(CWT_SCALES)

typedef struct {
   char *name, *seq;
   int length, genome_length, offset;
   float *cwt_power;
   uint8_t *label;
} Contig;

void anno_fail(const char *format, ...);
void *anno_alloc(size_t count, size_t size);
int base_index(char base);

#endif
