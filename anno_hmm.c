#include "anno_hmm.h"

#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  int count[HMM_STATES];
  int next[HMM_STATES][4];
} Edges;

int hmm_is_coding(int state) {
  return state >= HMM_C0 && state < HMM_C0 + 3;
}

int hmm_edge(int previous, int next) {
  if (previous == HMM_BACKGROUND)
    return next == HMM_BACKGROUND || hmm_is_coding(next);
  if (hmm_is_coding(previous)) {
    int phase = (previous - HMM_C0 + 1) % 3;
    return next == HMM_BACKGROUND || next == HMM_C0 + phase || next == HMM_DONOR(phase);
  }
  if (previous < 4 || previous >= HMM_STATES) return 0;
  int phase = (previous - 4) / 5;
  int kind = (previous - 4) % 5;
  if (kind == 2) return next == previous || next == previous + 1;
  if (kind == 4) return next == HMM_C0 + phase;
  return next == previous + 1;
}

static Edges edges_new(void) {
  Edges edges = {0};
  for (int previous = 0; previous < HMM_STATES; previous++)
    for (int next = 0; next < HMM_STATES; next++)
      if (hmm_edge(previous, next))
        edges.next[previous][edges.count[previous]++] = next;
  return edges;
}

static int score_bin(float probability) {
  if (!isfinite(probability) || probability < 0.0f || probability > 1.0f)
    anno_fail("Invalid CNN probability for HMM");
  int bin = (int)(probability * HMM_BINS);
  return bin == HMM_BINS ? HMM_BINS - 1 : bin;
}

static void emissions(const HMM *hmm, const HMMObservation *observation,
                       int position, double *values) {
  int base = base_index(observation->sequence[position]);
  int bin = score_bin(observation->probability[position]);
  for (int state = 0; state < HMM_STATES; state++)
    values[state] = fmax(hmm->score[state][bin] *
                        (base < 0 ? 1.0 : hmm->nucleotide[state][base]), 1e-300);
}

static void training_emissions(const HMM *hmm, const HMMObservation *observation,
                                int position, double *values) {
  emissions(hmm, observation, position, values);
  int label = observation->states ? observation->states[position] : HMM_UNKNOWN;
  if (label >= HMM_STATES) return;
  for (int state = 0; state < HMM_STATES; state++) {
    int compatible = label == 0 ? state == 0 : hmm_is_coding(label) ?
                     hmm_is_coding(state) : state >= 4 && (state - 4) % 5 == (label - 4) % 5;
    if (!compatible) values[state] = 0.0;
  }
}

static void normalize(double *values, int count) {
  double sum = 0.0;
  for (int index = 0; index < count; index++) sum += values[index];
  if (!(sum > 0.0) || !isfinite(sum)) anno_fail("Invalid HMM normalization");
  for (int index = 0; index < count; index++) values[index] /= sum;
}

static void model_normalize(HMM *hmm) {
  normalize(hmm->initial, HMM_STATES);
  for (int state = 0; state < HMM_STATES; state++) {
    normalize(hmm->transition[state], HMM_STATES);
    normalize(hmm->nucleotide[state], 4);
    normalize(hmm->score[state], HMM_BINS);
  }
}

void hmm_seed(HMM *hmm, const HMMObservation *observations, int count) {
  memset(hmm, 0, sizeof(*hmm));
  for (int state = 0; state < HMM_STATES; state++) {
    hmm->initial[state] = state < 4 ? 0.1 : 0.0;
    for (int next = 0; next < HMM_STATES; next++)
      hmm->transition[state][next] = hmm_edge(state, next) ? 0.1 : 0.0;
    for (int base = 0; base < 4; base++) hmm->nucleotide[state][base] = 0.1;
    for (int bin = 0; bin < HMM_BINS; bin++) hmm->score[state][bin] = 0.1;
  }
  for (int record = 0; record < count; record++) {
    const HMMObservation *observation = &observations[record];
    for (int position = 0; position < observation->length; position++) {
      int state = observation->states[position];
      if (state >= HMM_STATES) continue;
      if (!position && state < 4) hmm->initial[state]++;
      int base = base_index(observation->sequence[position]);
      if (base >= 0) hmm->nucleotide[state][base]++;
      hmm->score[state][score_bin(observation->probability[position])]++;
      if (position) {
        int previous = observation->states[position - 1];
        if (previous < HMM_STATES && hmm_edge(previous, state))
          hmm->transition[previous][state]++;
      }
    }
  }
  model_normalize(hmm);
}

double hmm_expectation(const HMM *hmm, const HMMObservation *observation,
                        HMM *counts) {
  int length = observation->length;
  if (length <= 0) return 0.0;
  Edges edges = edges_new();
  double *alpha = anno_alloc((size_t)length * HMM_STATES, sizeof(double));
  double *scaling = anno_alloc(length, sizeof(double));
  double likelihood = 0.0;
  for (int position = 0; position < length; position++) {
    double emission[HMM_STATES];
    training_emissions(hmm, observation, position, emission);
    double *row = alpha + (size_t)position * HMM_STATES;
    if (!position) {
      for (int state = 0; state < HMM_STATES; state++) row[state] = hmm->initial[state];
    } else {
      double *previous_row = row - HMM_STATES;
      for (int previous = 0; previous < HMM_STATES; previous++)
        for (int edge = 0; edge < edges.count[previous]; edge++) {
          int next = edges.next[previous][edge];
          row[next] += previous_row[previous] * hmm->transition[previous][next];
        }
    }
    for (int state = 0; state < HMM_STATES; state++) {
      row[state] *= emission[state];
      scaling[position] += row[state];
    }
    if (!(scaling[position] > 0.0)) anno_fail("HMM has no valid path");
    for (int state = 0; state < HMM_STATES; state++) row[state] /= scaling[position];
    likelihood += log(scaling[position]);
  }
  double terminal = 0.0;
  for (int state = 0; state < 4; state++)
    terminal += alpha[(size_t)(length - 1) * HMM_STATES + state];
  if (!(terminal > 0.0)) anno_fail("HMM has no complete terminal state");
  likelihood += log(terminal);
  if (counts) {
    double beta[HMM_STATES] = {0};
    for (int state = 0; state < 4; state++) beta[state] = 1.0 / terminal;
    for (int position = length - 1; position >= 0; position--) {
      const double *row = alpha + (size_t)position * HMM_STATES;
      int base = base_index(observation->sequence[position]);
      int bin = score_bin(observation->probability[position]);
      for (int state = 0; state < HMM_STATES; state++) {
        double posterior = row[state] * beta[state];
        if (!position) counts->initial[state] += posterior;
        if (base >= 0) counts->nucleotide[state][base] += posterior;
        counts->score[state][bin] += posterior;
      }
      if (position) {
        double emission[HMM_STATES], previous_beta[HMM_STATES] = {0};
        training_emissions(hmm, observation, position, emission);
        for (int previous = 0; previous < HMM_STATES; previous++)
          for (int edge = 0; edge < edges.count[previous]; edge++) {
            int next = edges.next[previous][edge];
            double contribution = hmm->transition[previous][next] * emission[next] *
                                  beta[next] / scaling[position];
            previous_beta[previous] += contribution;
            counts->transition[previous][next] += row[previous - HMM_STATES] * contribution;
          }
        memcpy(beta, previous_beta, sizeof(beta));
      }
    }
  }
  free(alpha);
  free(scaling);
  return likelihood;
}

static void maximize_row(double *destination, const double *counts, int length) {
  double sum = 0.0;
  for (int index = 0; index < length; index++) sum += counts[index];
  if (sum > 1e-200)
    for (int index = 0; index < length; index++) destination[index] = counts[index] / sum;
}

void hmm_fit(HMM *hmm, const HMMObservation *observations, int count, int iterations) {
  for (int iteration = 0; iteration < iterations; iteration++) {
    HMM counts = {0};
    double likelihood = 0.0;
    for (int record = 0; record < count; record++)
      likelihood += hmm_expectation(hmm, &observations[record], &counts);
    maximize_row(hmm->initial, counts.initial, HMM_STATES);
    for (int state = 0; state < HMM_STATES; state++) {
      maximize_row(hmm->transition[state], counts.transition[state], HMM_STATES);
      maximize_row(hmm->nucleotide[state], counts.nucleotide[state], 4);
      maximize_row(hmm->score[state], counts.score[state], HMM_BINS);
    }
    fprintf(stderr, "HMM EM %d/%d log-likelihood %.6f\n", iteration + 1, iterations, likelihood);
  }
}

unsigned char *hmm_decode(const HMM *hmm, const HMMObservation *observation) {
  int length = observation->length;
  unsigned char *path = anno_alloc(length, 1);
  if (length <= 0) return path;
  unsigned char *back = anno_alloc((size_t)length * HMM_STATES, 1);
  Edges edges = edges_new();
  double previous[HMM_STATES], transition[HMM_STATES][HMM_STATES];
  for (int state = 0; state < HMM_STATES; state++) {
    previous[state] = log(hmm->initial[state]);
    for (int next = 0; next < HMM_STATES; next++)
      transition[state][next] = log(hmm->transition[state][next]);
  }
  for (int position = 0; position < length; position++) {
    double emission[HMM_STATES], current[HMM_STATES];
    emissions(hmm, observation, position, emission);
    for (int state = 0; state < HMM_STATES; state++) current[state] = -INFINITY;
    if (!position) {
      memcpy(current, previous, sizeof(current));
    } else {
      for (int state = 0; state < HMM_STATES; state++)
        for (int edge = 0; edge < edges.count[state]; edge++) {
          int next = edges.next[state][edge];
          double value = previous[state] + transition[state][next];
          if (value > current[next]) {
            current[next] = value;
            back[(size_t)position * HMM_STATES + next] = (unsigned char)state;
          }
        }
    }
    double maximum = -INFINITY;
    for (int state = 0; state < HMM_STATES; state++) {
      current[state] += log(emission[state]);
      maximum = fmax(maximum, current[state]);
    }
    if (!isfinite(maximum)) anno_fail("No finite Viterbi path");
    for (int state = 0; state < HMM_STATES; state++) previous[state] = current[state] - maximum;
  }
  int final_state = 0;
  for (int state = 1; state < 4; state++)
    if (previous[state] > previous[final_state]) final_state = state;
  if (!isfinite(previous[final_state])) anno_fail("No complete Viterbi path");
  path[length - 1] = (unsigned char)final_state;
  for (int position = length - 1; position > 0; position--)
    path[position - 1] = back[(size_t)position * HMM_STATES + path[position]];
  free(back);
  return path;
}

static void write_row(const double *row, int count, FILE *stream) {
  for (int index = 0; index < count; index++)
    fprintf(stream, "%.17g%c", row[index], index + 1 == count ? '\n' : ' ');
}

void hmm_write(const HMM *hmm, FILE *stream) {
  write_row(hmm->initial, HMM_STATES, stream);
  for (int state = 0; state < HMM_STATES; state++) {
    write_row(hmm->transition[state], HMM_STATES, stream);
    write_row(hmm->nucleotide[state], 4, stream);
    write_row(hmm->score[state], HMM_BINS, stream);
  }
}

static void read_row(double *row, int count, FILE *stream) {
  double sum = 0.0;
  for (int index = 0; index < count; index++) {
    if (fscanf(stream, "%lf", &row[index]) != 1 || !isfinite(row[index]) || row[index] < 0.0)
      anno_fail("Invalid or truncated HMM parameters");
    sum += row[index];
  }
  if (fabs(sum - 1.0) > 1e-6) anno_fail("Unnormalized HMM parameters");
}

void hmm_read(HMM *hmm, FILE *stream) {
  read_row(hmm->initial, HMM_STATES, stream);
  for (int state = 0; state < HMM_STATES; state++) {
    read_row(hmm->transition[state], HMM_STATES, stream);
    read_row(hmm->nucleotide[state], 4, stream);
    read_row(hmm->score[state], HMM_BINS, stream);
    if (state >= 4 && hmm->initial[state] != 0.0) anno_fail("Invalid HMM initial topology");
    for (int next = 0; next < HMM_STATES; next++)
      if (!hmm_edge(state, next) && hmm->transition[state][next] != 0.0)
        anno_fail("Invalid HMM splice topology");
  }
}

static char *gff_identifier(const char *name) {
  char *escaped = anno_alloc(strlen(name) * 3 + 1, 1);
  char *destination = escaped;
  for (const unsigned char *source = (const unsigned char *)name; *source; source++) {
    if (isalnum(*source) || strchr(".:^*$@!+_?-|", *source)) *destination++ = (char)*source;
    else { snprintf(destination, 4, "%%%02X", (unsigned)*source); destination += 3; }
  }
  return escaped;
}

void hmm_gff(FILE *stream, const char *name, int length, int strand,
             const unsigned char *path, unsigned long *gene_number) {
  char *identifier = gff_identifier(name);
  if (!strand) fprintf(stream, "##sequence-region %s 1 %d\n", identifier, length);
  int position = 0;
  while (position < length) {
    if (!hmm_is_coding(path[position])) { position++; continue; }
    int start = position, last_coding = position;
    while (position < length && path[position] != HMM_BACKGROUND) {
      if (hmm_is_coding(path[position])) last_coding = position;
      position++;
    }
    int end = last_coding + 1;
    int genomic_start = strand ? length - end + 1 : start + 1;
    int genomic_end = strand ? length - start : end;
    char direction = strand ? '-' : '+';
    unsigned long gene = ++*gene_number;
    fprintf(stream, "%s\tMorletCNN\tgene\t%d\t%d\t.\t%c\t.\tID=gene%lu;gene_biotype=protein_coding;prediction_status=unvalidated\n",
            identifier, genomic_start, genomic_end, direction, gene);
    fprintf(stream, "%s\tMorletCNN\tmRNA\t%d\t%d\t.\t%c\t.\tID=tx%lu;Parent=gene%lu\n",
            identifier, genomic_start, genomic_end, direction, gene, gene);
    int block = 0;
    for (int current = start; current < end;) {
      if (!hmm_is_coding(path[current])) { current++; continue; }
      int cds_start = current;
      int phase = (3 - (path[current] - HMM_C0)) % 3;
      while (current < end && hmm_is_coding(path[current])) current++;
      int left = strand ? length - current + 1 : cds_start + 1;
      int right = strand ? length - cds_start : current;
      block++;
      fprintf(stream, "%s\tMorletCNN\texon\t%d\t%d\t.\t%c\t.\tID=exon%lu_%d;Parent=tx%lu\n",
              identifier, left, right, direction, gene, block, gene);
      fprintf(stream, "%s\tMorletCNN\tCDS\t%d\t%d\t.\t%c\t%d\tID=cds%lu_%d;Parent=tx%lu\n",
              identifier, left, right, direction, phase, gene, block, gene);
    }
  }
  free(identifier);
}