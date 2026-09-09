#define main anno_cli_main
#include "anno.c"
#undef main
#include <assert.h>

typedef struct {
  unsigned long positive, predicted, correct, total;
  unsigned long reference_blocks, predicted_blocks, matched_blocks;
  unsigned long reference_introns, predicted_introns, matched_introns;
  double loss;
} Metrics;

static int interval_state(int state, int coding) {
  return coding ? hmm_is_coding(state) : state >= 4 && state < HMM_STATES;
}

static void intervals(const unsigned char *truth, const unsigned char *prediction,
                       int length, int coding, unsigned long *reference_count,
                       unsigned long *prediction_count, unsigned long *matched) {
  for (int position = 0; position < length; position++) {
    if (interval_state(truth[position], coding) &&
        (!position || !interval_state(truth[position - 1], coding)))
      (*reference_count)++;
    if (!interval_state(prediction[position], coding) ||
        (position && interval_state(prediction[position - 1], coding))) continue;
    (*prediction_count)++;
    int end = position + 1;
    while (end < length && interval_state(prediction[end], coding)) end++;
    int exact = (!position || !interval_state(truth[position - 1], coding)) &&
                (end == length || !interval_state(truth[end], coding));
    for (int inside = position; inside < end && exact; inside++)
      exact = interval_state(truth[inside], coding);
    if (exact) (*matched)++;
  }
}

static void measure(Metrics *metrics, const Sequence *sequence, int strand,
                     const float *probability, const unsigned char *prediction) {
  for (int position = 0; position < sequence->length; position++) {
    int target = sequence->coding[strand][position];
    int decision = prediction ? hmm_is_coding(prediction[position]) : probability[position] >= 0.5f;
    metrics->positive += target;
    metrics->predicted += decision;
    metrics->correct += target && decision;
    metrics->total++;
    double score = target ? probability[position] : 1.0 - probability[position];
    metrics->loss -= log(fmax(score, 1e-30));
  }
  if (prediction) {
    intervals(sequence->states[strand], prediction, sequence->length, 1,
              &metrics->reference_blocks, &metrics->predicted_blocks, &metrics->matched_blocks);
    intervals(sequence->states[strand], prediction, sequence->length, 0,
              &metrics->reference_introns, &metrics->predicted_introns, &metrics->matched_introns);
  }
}

static void report(const char *name, const Metrics *metrics) {
  double precision = metrics->predicted ? (double)metrics->correct / metrics->predicted : 0.0;
  double recall = metrics->positive ? (double)metrics->correct / metrics->positive : 0.0;
  double f1 = precision + recall ? 2.0 * precision * recall / (precision + recall) : 0.0;
  printf("%s: CE=%.6f CDS precision=%.6f recall=%.6f F1=%.6f\n",
         name, metrics->loss / metrics->total, precision, recall, f1);
  if (!metrics->reference_blocks && !metrics->predicted_blocks) return;
  printf("%s: exact CDS=%lu reference=%lu predicted=%lu; exact introns=%lu reference=%lu predicted=%lu\n",
         name, metrics->matched_blocks, metrics->reference_blocks, metrics->predicted_blocks,
         metrics->matched_introns, metrics->reference_introns, metrics->predicted_introns);
}

int main(int argc, char **argv) {
  if (argc == 2 && !strcmp(argv[1], "--self-test")) {
    unsigned char truth[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 0};
    unsigned char shifted[] = {0, 1, 2, 3, 1, 4, 5, 6, 7, 8, 2, 3, 0};
    unsigned long reference_count = 0, prediction_count = 0, matched = 0;
    intervals(truth, truth, 13, 1, &reference_count, &prediction_count, &matched);
    assert(reference_count == 2 && prediction_count == 2 && matched == 2);
    reference_count = prediction_count = matched = 0;
    intervals(truth, shifted, 13, 1, &reference_count, &prediction_count, &matched);
    assert(reference_count == 2 && prediction_count == 2 && matched == 0);
    reference_count = prediction_count = matched = 0;
    intervals(truth, truth, 13, 0, &reference_count, &prediction_count, &matched);
    assert(reference_count == 1 && prediction_count == 1 && matched == 1);
    reference_count = prediction_count = matched = 0;
    intervals(truth, shifted, 13, 0, &reference_count, &prediction_count, &matched);
    assert(reference_count == 1 && prediction_count == 1 && matched == 0);
    puts("exact interval evaluation tests passed");
    return 0;
  }
  int check = argc > 1 && !strcmp(argv[1], "--check");
  if (check) { argc--; argv++; }
  if (argc < 4 || argc > 6) anno_fail("Usage: test_overfit [--check] FASTA GFF MODEL_PREFIX [REFERENCE_GFF [PREDICTED_GFF]]");
  Dataset *dataset = dataset_read(argv[1]);
  FILE *annotation = fopen(argv[2], "r");
  if (!annotation) anno_fail("Cannot read annotation");
  annotation_read(dataset, annotation);
  fclose(annotation);
  Model model = model_load(argv[3]);
  if (model.stage != 2) anno_fail("Evaluation requires a fine-tuned model");
  Wavelets *wavelets = anno_alloc(1, sizeof(*wavelets));
  wavelets_init(wavelets);
  HMMObservation *observations = anno_alloc(dataset->count * 2, sizeof(*observations));
  for (int record = 0; record < dataset->count; record++) {
    Sequence *sequence = &dataset->sequences[record];
    for (int strand = 0; strand < 2; strand++) {
      HMMObservation *observation = &observations[2 * record + strand];
      observation->sequence = strand ? reverse_complement(sequence->bases, sequence->length) : sequence->bases;
      observation->length = sequence->length;
      observation->states = sequence->states[strand];
      observation->probability = network_predict(model.network, wavelets, observation->sequence, observation->length);
    }
  }
  HMM seeded;
  hmm_seed(&seeded, observations, dataset->count * 2);
  HMM fitted = seeded;
  hmm_fit(&fitted, observations, dataset->count * 2, 10);
  FILE *reference = argc > 4 ? fopen(argv[4], "w") : NULL;
  FILE *prediction = argc > 5 ? fopen(argv[5], "w") : NULL;
  if ((argc > 4 && !reference) || (argc > 5 && !prediction)) anno_fail("Cannot write evaluation GFF");
  if (reference) fputs("##gff-version 3\n", reference);
  if (prediction) fputs("##gff-version 3\n", prediction);
  unsigned long reference_genes = 0, predicted_genes = 0;
  Metrics raw = {0}, saved = {0}, initialized = {0}, constrained = {0};
  for (int record = 0; record < dataset->count; record++) {
    Sequence *sequence = &dataset->sequences[record];
    for (int strand = 0; strand < 2; strand++) {
      HMMObservation *observation = &observations[2 * record + strand];
      HMMObservation unlabeled = *observation;
      unlabeled.states = NULL;
      unsigned char *fitted_path = hmm_decode(&fitted, &unlabeled);
      unsigned char *path = hmm_decode(&model.hmm, &unlabeled);
      unsigned char *seeded_path = hmm_decode(&seeded, &unlabeled);
      measure(&raw, sequence, strand, observation->probability, NULL);
      measure(&saved, sequence, strand, observation->probability, path);
      measure(&initialized, sequence, strand, observation->probability, seeded_path);
      measure(&constrained, sequence, strand, observation->probability, fitted_path);
      if (reference) hmm_gff(reference, sequence->name, sequence->length, strand, observation->states, &reference_genes);
      if (prediction) hmm_gff(prediction, sequence->name, sequence->length, strand, path, &predicted_genes);
      free(path); free(seeded_path); free(fitted_path);
      free((void *)observation->probability);
      if (strand) free((void *)observation->sequence);
    }
  }
  if (reference) fclose(reference);
  if (prediction) fclose(prediction);
  double prior = (double)raw.positive / raw.total;
    double baseline = prior > 0.0 && prior < 1.0 ?
            -prior * log(prior) - (1.0 - prior) * log1p(-prior) : 0.0;
  printf("labels: positive=%lu total=%lu prior=%.6f baseline_CE=%.6f\n",
      raw.positive, raw.total, prior, baseline);
  report("CNN", &raw);
  report("saved HMM", &saved);
  report("seeded HMM", &initialized);
  report("constrained EM", &constrained);
  free(observations);
  free(wavelets);
  cnn_free(model.network);
  dataset_free(dataset);
    if (!check) return 0;
    double base_f1 = saved.positive + saved.predicted ?
      2.0 * saved.correct / (saved.positive + saved.predicted) : 0.0;
    double cds_f1 = saved.reference_blocks + saved.predicted_blocks ?
      2.0 * saved.matched_blocks / (saved.reference_blocks + saved.predicted_blocks) : 0.0;
    double intron_f1 = saved.reference_introns + saved.predicted_introns ?
      2.0 * saved.matched_introns / (saved.reference_introns + saved.predicted_introns) : 0.0;
    int passed = base_f1 >= 0.99 && cds_f1 >= 0.80 && intron_f1 >= 0.90;
    printf("Overfit gate %s: base F1 %.6f >= .99; exact CDS F1 %.6f >= .80; exact intron F1 %.6f >= .90\n",
       passed ? "PASS" : "FAIL", base_f1, cds_f1, intron_f1);
    return passed ? 0 : 1;
}