#include "anno_hmm.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static double brute_sum, brute_best;
static unsigned char brute_path[8], candidate[8];

static void enumerate(const HMM *hmm, const HMMObservation *observation,
                       int position, int previous, double probability) {
  if (position == observation->length) {
    if (previous >= 4) return;
    brute_sum += probability;
    if (probability > brute_best) {
      brute_best = probability;
      memcpy(brute_path, candidate, sizeof(brute_path));
    }
    return;
  }
  int bin = (int)(observation->probability[position] * HMM_BINS);
  for (int state = 0; state < HMM_STATES; state++) {
    if (observation->states && observation->states[position] < HMM_STATES) {
      int label = observation->states[position];
      int allowed = label == 0 ? state == 0 : hmm_is_coding(label) ?
                    hmm_is_coding(state) : state >= 4 && (state - 4) % 5 == (label - 4) % 5;
      if (!allowed) continue;
    }
    double transition = position ? hmm->transition[previous][state] : hmm->initial[state];
    if (!transition) continue;
    candidate[position] = state;
    double emission = hmm->score[state][bin] * hmm->nucleotide[state][base_index(observation->sequence[position])];
    enumerate(hmm, observation, position + 1, state, probability * transition * emission);
  }
}

int main(void) {
  HMM hmm;
  const char *sequence = "ACGTACGT";
  float probability[8] = {.1, .9, .9, .1, .1, .1, .9, .1};
  unsigned char labels[8] = {0, 1, 2, 3, 0, 0, 0, 0};
  HMMObservation observation = {sequence, 8, probability, labels};
  hmm_seed(&hmm, &observation, 1);
  observation.states = NULL;
  enumerate(&hmm, &observation, 0, 0, 1.0);
  HMM counts = {0};
  double likelihood = hmm_expectation(&hmm, &observation, &counts);
  assert(fabs(likelihood - log(brute_sum)) < 1e-10);
  double emissions = 0.0, transitions = 0.0, initial = 0.0;
  for (int state = 0; state < HMM_STATES; state++) {
    initial += counts.initial[state];
    for (int bin = 0; bin < HMM_BINS; bin++) emissions += counts.score[state][bin];
    for (int next = 0; next < HMM_STATES; next++) transitions += counts.transition[state][next];
  }
  assert(fabs(initial - 1.0) < 1e-10);
  assert(fabs(emissions - 8.0) < 1e-10 && fabs(transitions - 7.0) < 1e-10);
  unsigned char *path = hmm_decode(&hmm, &observation);
  assert(!memcmp(path, brute_path, 8));
  free(path);
  observation.states = labels;
  brute_sum = brute_best = 0.0;
  enumerate(&hmm, &observation, 0, 0, 1.0);
  memset(&counts, 0, sizeof(counts));
  likelihood = hmm_expectation(&hmm, &observation, &counts);
  assert(fabs(likelihood - log(brute_sum)) < 1e-10);
  double coding_count = 0.0;
  for (int state = HMM_C0; state < HMM_C0 + 3; state++)
    for (int bin = 0; bin < HMM_BINS; bin++) coding_count += counts.score[state][bin];
  assert(fabs(coding_count - 3.0) < 1e-10);
  for (int iteration = 0; iteration < 3; iteration++) {
    hmm_fit(&hmm, &observation, 1, 1);
    double updated = hmm_expectation(&hmm, &observation, NULL);
    assert(updated + 1e-8 >= likelihood);
    likelihood = updated;
  }
  FILE *saved = tmpfile();
  assert(saved);
  hmm_write(&hmm, saved);
  rewind(saved);
  HMM restored;
  hmm_read(&restored, saved);
  assert(!memcmp(&hmm, &restored, sizeof(hmm)));
  fclose(saved);
  unsigned char spliced[20] = {0, 1, 2, 3, 1, 2, 3, 1, 9, 10, 11, 12, 13, 2, 3, 1, 2, 3, 0, 0};
  for (int position = 1; position < 20; position++) assert(hmm_edge(spliced[position - 1], spliced[position]));
  char splice_sequence[21] = "AAAAAAAAAAAAAAAAAAAA";
  splice_sequence[8] = 'G'; splice_sequence[9] = 'T';
  splice_sequence[11] = 'A'; splice_sequence[12] = 'G';
  float splice_probability[20];
  for (int position = 0; position < 20; position++)
    splice_probability[position] = hmm_is_coding(spliced[position]) ? .95f : .05f;
  HMMObservation repeats[100];
  for (int record = 0; record < 100; record++)
    repeats[record] = (HMMObservation){splice_sequence, 20, splice_probability, spliced};
  hmm_seed(&hmm, repeats, 100);
  assert(hmm.nucleotide[HMM_DONOR(1)][2] > .99);
  assert(hmm.nucleotide[HMM_DONOR(1) + 1][3] > .99);
  assert(hmm.nucleotide[HMM_DONOR(1) + 3][0] > .99);
  assert(hmm.nucleotide[HMM_DONOR(1) + 4][2] > .99);
  hmm_fit(&hmm, repeats, 100, 2);
  repeats[0].states = NULL;
  path = hmm_decode(&hmm, &repeats[0]);
  assert(!memcmp(path, spliced, 20));
  free(path);
  saved = tmpfile();
  unsigned long gene = 0;
  hmm_gff(saved, "chr;%", 20, 1, spliced, &gene);
  rewind(saved);
  char output[4096] = {0};
  assert(fread(output, 1, sizeof(output) - 1, saved) > 0);
  assert(strstr(output, "chr%3B%25\tMorletCNN"));
  assert(strstr(output, "gene\t3\t19\t.\t-\t."));
  assert(strstr(output, "CDS\t13\t19\t.\t-\t0"));
  assert(strstr(output, "CDS\t3\t7\t.\t-\t2"));
  assert(strstr(output, "gene_biotype=protein_coding"));
  fclose(saved);
  puts("HMM brute-force, EM, topology, serialization and GFF tests passed");
  return 0;
}