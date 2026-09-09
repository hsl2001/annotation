#ifndef ANNO_HMM_H
#define ANNO_HMM_H

#include "anno_signal.h"
#include <stdio.h>

#define HMM_STATES 19
#define HMM_BINS 8
#define HMM_BACKGROUND 0
#define HMM_C0 1
#define HMM_DONOR(phase) (4 + 5 * (phase))
#define HMM_UNKNOWN 255

typedef struct {
  double initial[HMM_STATES];
  double transition[HMM_STATES][HMM_STATES];
  double nucleotide[HMM_STATES][4];
  double score[HMM_STATES][HMM_BINS];
} HMM;

typedef struct {
  const char *sequence;
  int length;
  const float *probability;
  const unsigned char *states;
} HMMObservation;

int hmm_is_coding(int state);
int hmm_edge(int previous, int next);
void hmm_seed(HMM *hmm, const HMMObservation *observations, int count);
double hmm_expectation(const HMM *hmm, const HMMObservation *observation,
                        HMM *counts);
void hmm_fit(HMM *hmm, const HMMObservation *observations, int count, int iterations);
unsigned char *hmm_decode(const HMM *hmm, const HMMObservation *observation);
void hmm_write(const HMM *hmm, FILE *stream);
void hmm_read(HMM *hmm, FILE *stream);
void hmm_gff(FILE *stream, const char *name, int length, int strand,
             const unsigned char *path, unsigned long *gene_number);

#endif