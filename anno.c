#define _POSIX_C_SOURCE 200809L
#include "anno_cnn.h"
#include "anno_data.h"
#include "anno_hmm.h"

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#define MODEL_VERSION 1

typedef struct {
  CNN *network;
  HMM hmm;
  int stage;
  unsigned long pretrain_steps, train_steps;
  unsigned seed;
} Model;

static char *model_path(const char *prefix, const char *suffix) {
  size_t size = strlen(prefix) + strlen(suffix) + 1;
  char *path = anno_alloc(size, 1);
  snprintf(path, size, "%s%s", prefix, suffix);
  return path;
}

static void model_save(Model *model, const char *prefix) {
  char *path = model_path(prefix, ".bin");
  char *temporary = model_path(prefix, ".bin.tmp.XXXXXX");
  int descriptor = mkstemp(temporary);
  if (descriptor < 0) anno_fail("Cannot create model: %s", path);
  FILE *stream = fdopen(descriptor, "w");
  if (!stream) anno_fail("Cannot open model stream: %s", path);
  fprintf(stream, "ANNO_MORLET_CNN %d\n%d %d %d %d\n", MODEL_VERSION,
          CNN_DEPTH, CNN_WIDTH, WAVE_COUNT, HMM_STATES);
  for (int scale = 0; scale < WAVE_COUNT; scale++)
    fprintf(stream, "%d%c", wave_sizes[scale], scale + 1 == WAVE_COUNT ? '\n' : ' ');
  fprintf(stream, "%d %lu %lu %u\n", model->stage, model->pretrain_steps,
          model->train_steps, model->seed);
  cnn_write(model->network, stream);
  if (model->stage == 2) hmm_write(&model->hmm, stream);
  int failed = ferror(stream);
  if (fclose(stream)) failed = 1;
  if (failed || rename(temporary, path)) {
    unlink(temporary);
    anno_fail("Failed to save model: %s", path);
  }
  fprintf(stderr, "Saved %s\n", path);
  free(path);
  free(temporary);
}

static Model model_load(const char *prefix) {
  char *path = model_path(prefix, ".bin");
  FILE *stream = fopen(path, "r");
  if (!stream) anno_fail("Cannot open model: %s", path);
  char magic[64];
  int version, depth, width, scales, states;
  if (fscanf(stream, "%63s %d", magic, &version) != 2 ||
      strcmp(magic, "ANNO_MORLET_CNN") || version != MODEL_VERSION)
    anno_fail("Unsupported model (old MLP-Mixer models cannot be loaded): %s", path);
  if (fscanf(stream, "%d %d %d %d", &depth, &width, &scales, &states) != 4 ||
      depth != CNN_DEPTH || width != CNN_WIDTH || scales != WAVE_COUNT || states != HMM_STATES)
    anno_fail("Incompatible model architecture: %s", path);
  for (int scale = 0; scale < WAVE_COUNT; scale++) {
    int size;
    if (fscanf(stream, "%d", &size) != 1 || size != wave_sizes[scale])
      anno_fail("Incompatible Morlet kernel sizes: %s", path);
  }
  Model model = {0};
  if (fscanf(stream, "%d %lu %lu %u", &model.stage, &model.pretrain_steps,
             &model.train_steps, &model.seed) != 4 ||
      (model.stage != 1 && model.stage != 2) || !model.pretrain_steps ||
      (model.stage == 2 && !model.train_steps))
    anno_fail("Invalid model training metadata: %s", path);
  model.network = cnn_new();
  cnn_read(model.network, stream);
  if (model.stage == 2) hmm_read(&model.hmm, stream);
  int trailing;
  do { trailing = fgetc(stream); } while (trailing != EOF && isspace((unsigned char)trailing));
  if (trailing != EOF || ferror(stream)) anno_fail("Unexpected model payload: %s", path);
  fclose(stream);
  free(path);
  return model;
}

static unsigned char *make_mask(const char *sequence, int length) {
  unsigned char *mask = anno_alloc(length, 1);
  int first_known = -1, masked = 0;
  for (int position = 0; position < length; position++) {
    if (base_index(sequence[position]) < 0) continue;
    if (first_known < 0) first_known = position;
    if ((double)rand() / ((double)RAND_MAX + 1.0) < 0.15) {
      mask[position] = 1;
      masked++;
    }
  }
  if (!masked && first_known >= 0) mask[first_known] = 1;
  return mask;
}

static void network_train(Model *model, Dataset *dataset, const Wavelets *wavelets,
                           int task, int epochs, double rate) {
  int count = TILE_SIZE + 2 * CNN_HALO;
  CNNCache *cache = cnn_cache_new(count);
  CNNGradient *gradient = cnn_gradient_new(model->network, task);
  CNNOptimizer *optimizer = task ? cnn_optimizer_new(model->network, task) : NULL;
  CNNGradient *batch = task ? cnn_gradient_new(model->network, task) : NULL;
  double **feature_cache = anno_alloc((size_t)dataset->count * 2, sizeof(*feature_cache));
  size_t total_bases = 0;
  for (int record = 0; record < dataset->count; record++) total_bases += dataset->sequences[record].length;
  if (task) {
    size_t positives = 0;
    for (int record = 0; record < dataset->count; record++)
      for (int strand = 0; strand < 2; strand++)
        for (int position = 0; position < dataset->sequences[record].length; position++)
          positives += dataset->sequences[record].coding[strand][position];
    double prior = (double)positives / (2.0 * total_bases);
    double baseline = prior > 0.0 && prior < 1.0 ?
                      -prior * log(prior) - (1.0 - prior) * log1p(-prior) : 0.0;
    fprintf(stderr, "Training labels: %zu positive / %zu strand-positions; prior %.6f; baseline CE %.6f\n",
            positives, total_bases * 2, prior, baseline);
  }
  int cache_features = task && total_bases <= 2000000;
  double *features = anno_alloc((size_t)count * CWT_CHANNELS, sizeof(double));
  int *targets = anno_alloc(count, sizeof(int));
  for (int epoch = 0; epoch < epochs; epoch++) {
    clock_t started = clock();
    double total_loss = 0.0;
    unsigned long examples = 0;
    double epoch_rate = task && epochs > 1 ? rate * pow(0.1, (double)epoch / (epochs - 1)) : rate;
    for (int record = 0; record < dataset->count; record++) {
      Sequence *sequence = &dataset->sequences[record];
      int tiles = (sequence->length + TILE_SIZE - 1) / TILE_SIZE;
      int *order = anno_alloc((size_t)tiles * 2, sizeof(int));
      char *reverse = reverse_complement(sequence->bases, sequence->length);
      const char *oriented[2] = {sequence->bases, reverse};
      unsigned char *masks[2] = {task ? NULL : make_mask(oriented[0], sequence->length),
                                 task ? NULL : make_mask(oriented[1], sequence->length)};
      if (cache_features && !epoch) {
        for (int strand = 0; strand < 2; strand++) {
          feature_cache[2 * record + strand] = anno_alloc((size_t)sequence->length * CWT_CHANNELS, sizeof(double));
          cwt_extract(wavelets, oriented[strand], sequence->length, NULL, 0, sequence->length,
                      feature_cache[2 * record + strand]);
        }
      }
        for (int tile = 0; tile < 2 * tiles; tile++) order[tile] = tile;
        for (int tile = 2 * tiles - 1; tile > 0; tile--) {
          int other = rand() % (tile + 1);
          int swap = order[tile]; order[tile] = order[other]; order[other] = swap;
        }
        fprintf(stderr, "%s epoch %d/%d: %s (both strands), %d bp\n",
                task ? "Train" : "Pretrain", epoch + 1, epochs,
          sequence->name, sequence->length);
              int batch_targets = 0, batch_tiles = 0;
        for (int tile = 0; tile < 2 * tiles; tile++) {
          int strand = order[tile] / tiles;
          int start = (order[tile] % tiles) * TILE_SIZE;
          const char *bases = oriented[strand];
          unsigned char *mask = masks[strand];
          int labeled = 0;
          for (int position = 0; position < count; position++) targets[position] = -1;
          for (int offset = 0; offset < TILE_SIZE && start + offset < sequence->length; offset++) {
            int position = start + offset;
            if (task || mask[position]) {
              targets[CNN_HALO + offset] = task ? sequence->coding[strand][position] : base_index(bases[position]);
              labeled++;
            }
          }
          if (!labeled) continue;
          if (cache_features) {
            int origin = start - CNN_HALO;
            int first = origin < 0 ? -origin : 0;
            int end = origin + count > sequence->length ? sequence->length - origin : count;
            memset(features, 0, (size_t)count * CWT_CHANNELS * sizeof(double));
            memcpy(features + (size_t)first * CWT_CHANNELS,
                   feature_cache[2 * record + strand] + (size_t)(origin + first) * CWT_CHANNELS,
                   (size_t)(end - first) * CWT_CHANNELS * sizeof(double));
          } else {
            cwt_extract(wavelets, bases, sequence->length, mask, start - CNN_HALO, count, features);
          }
          cnn_forward(model->network, cache, features, start - CNN_HALO, sequence->length);
          double loss = cnn_loss(model->network, cache, targets, task, gradient);
          if (!isfinite(loss)) anno_fail("Non-finite training loss");
          if (task) {
            cnn_gradient_add(model->network, batch, gradient, task, labeled);
            batch_targets += labeled;
            batch_tiles++;
            if (batch_tiles == 8 || tile + 1 == 2 * tiles) {
              cnn_gradient_scale(model->network, batch, task, 1.0 / batch_targets);
              cnn_adam_update(model->network, batch, optimizer, task, epoch_rate);
              cnn_gradient_scale(model->network, batch, task, 0.0);
              batch_tiles = batch_targets = 0;
              model->train_steps++;
            }
          } else {
            cnn_update(model->network, gradient, task, rate);
            model->pretrain_steps++;
          }
          total_loss += loss * labeled;
          examples += labeled;
        }
      free(masks[0]); free(masks[1]);
      free(reverse);
      free(order);
    }
    if (!examples) anno_fail("No usable training targets (pretraining needs A/C/G/T)");
    double elapsed = (double)(clock() - started) / CLOCKS_PER_SEC;
        fprintf(stderr, "%s epoch %d/%d mean cross-entropy %.6f; %lu targets; lr %.6g; %.2f CPU seconds\n",
            task ? "Train" : "Pretrain", epoch + 1, epochs,
          total_loss / examples, examples, epoch_rate, elapsed);
  }
  free(targets);
  free(features);
  if (optimizer) cnn_optimizer_free(optimizer);
  if (batch) cnn_gradient_free(batch);
  for (int index = 0; index < 2 * dataset->count; index++) free(feature_cache[index]);
  free(feature_cache);
  cnn_gradient_free(gradient);
  cnn_cache_free(cache);
}

static float *network_predict(CNN *network, const Wavelets *wavelets,
                               const char *sequence, int length) {
  int count = TILE_SIZE + 2 * CNN_HALO;
  CNNCache *cache = cnn_cache_new(count);
  double *features = anno_alloc((size_t)count * CWT_CHANNELS, sizeof(double));
  float *probability = anno_alloc(length, sizeof(float));
  for (int start = 0; start < length; start += TILE_SIZE) {
    cwt_extract(wavelets, sequence, length, NULL, start - CNN_HALO, count, features);
    cnn_forward(network, cache, features, start - CNN_HALO, length);
    for (int offset = 0; offset < TILE_SIZE && start + offset < length; offset++) {
      double output[2];
      cnn_probabilities(network, cache, CNN_HALO + offset, 1, output);
      probability[start + offset] = (float)output[1];
    }
  }
  free(features);
  cnn_cache_free(cache);
  return probability;
}

static void train_hmm(Model *model, Dataset *dataset, const Wavelets *wavelets,
                       int iterations) {
  HMMObservation *observations = anno_alloc((size_t)dataset->count * 2, sizeof(*observations));
  unsigned long donors = 0;
  for (int record = 0; record < dataset->count; record++) {
    Sequence *sequence = &dataset->sequences[record];
    for (int strand = 0; strand < 2; strand++) {
      HMMObservation *observation = &observations[record * 2 + strand];
      observation->sequence = strand ? reverse_complement(sequence->bases, sequence->length) : sequence->bases;
      observation->length = sequence->length;
      observation->states = sequence->states[strand];
      fprintf(stderr, "HMM observations: %s (%c)\n", sequence->name, strand ? '-' : '+');
      observation->probability = network_predict(model->network, wavelets, observation->sequence, observation->length);
      for (int position = 0; position < observation->length; position++)
        for (int phase = 0; phase < 3; phase++)
          if (observation->states[position] == HMM_DONOR(phase)) donors++;
    }
  }
  if (!donors) fprintf(stderr, "Warning: no selected multi-CDS splice junctions; splice training is underdetermined\n");
  else fprintf(stderr, "HMM initialization: %lu annotated splice junctions\n", donors);
  hmm_seed(&model->hmm, observations, dataset->count * 2);
  hmm_fit(&model->hmm, observations, dataset->count * 2, iterations);
  for (int record = 0; record < dataset->count * 2; record++) {
    free((void *)observations[record].probability);
    if (record % 2) free((void *)observations[record].sequence);
  }
  free(observations);
}

static int integer_option(const char *text, const char *name, int minimum) {
  char *end;
  errno = 0;
  long value = strtol(text, &end, 10);
  if (errno || !*text || *end || value < minimum || value > INT_MAX)
    anno_fail("Invalid %s: %s", name, text);
  return (int)value;
}

static double rate_option(const char *text) {
  char *end;
  errno = 0;
  double rate = strtod(text, &end);
  if (errno || !*text || *end || !isfinite(rate) || rate <= 0.0 || rate > 1.0)
    anno_fail("Learning rate must be finite and in (0, 1]: %s", text);
  return rate;
}

static void usage(const char *program) {
  fprintf(stderr,
          "Complex Morlet CWT + genann 8-layer CNN + splice HMM\n"
          "Usage:\n"
          "  %s pretrain <genome.fasta> <out_prefix> [epochs=3] [lr=0.01] [seed=1]\n"
          "  %s train <genome.fasta> <annotation.gff> <pretrained_prefix> <out_prefix> [epochs=3] [lr=0.001] [em_iterations=5] [seed=1]\n"
          "  %s predict <genome.fasta> <trained_prefix>\n",
          program, program, program);
}

int main(int argc, char **argv) {
  if (argc == 2 && (!strcmp(argv[1], "--help") || !strcmp(argv[1], "-h"))) {
    usage(argv[0]);
    return 0;
  }
  if (argc < 2) { usage(argv[0]); return 1; }
  int pretrain = !strcmp(argv[1], "pretrain");
  int train = !strcmp(argv[1], "train");
  int predict = !strcmp(argv[1], "predict");
  if ((!pretrain && !train && !predict) || (pretrain && (argc < 4 || argc > 7)) ||
      (train && (argc < 6 || argc > 10)) || (predict && argc != 4)) {
    usage(argv[0]);
    return 1;
  }
  int option_start = pretrain ? 4 : 6;
  int epochs = !predict && argc > option_start ? integer_option(argv[option_start], "epochs", 1) : 3;
  double rate = !predict && argc > option_start + 1 ? rate_option(argv[option_start + 1]) : (train ? 0.001 : 0.01);
  int em_iterations = train && argc > 8 ? integer_option(argv[8], "EM iterations", 1) : 5;
  int seed_index = pretrain ? 6 : 9;
  unsigned seed = !predict && argc > seed_index ? (unsigned)integer_option(argv[seed_index], "seed", 0) : 1;
  srand(seed);
  Model model = {0};
  if (pretrain) model.network = cnn_new();
  else model = model_load(argv[train ? 4 : 3]);
  if (predict && model.stage != 2)
    anno_fail("Prediction requires a fine-tuned CNN and fitted HMM; run train first");
  if (!predict) model.seed = seed;
  Dataset *dataset = dataset_read(argv[2]);
  Wavelets *wavelets = anno_alloc(1, sizeof(*wavelets));
  wavelets_init(wavelets);
  if (pretrain) {
    network_train(&model, dataset, wavelets, 0, epochs, rate);
    model.stage = 1;
    model_save(&model, argv[3]);
  } else if (train) {
    FILE *annotation = fopen(argv[3], "r");
    if (!annotation) anno_fail("Cannot open annotation: %s", argv[3]);
    annotation_read(dataset, annotation);
    fclose(annotation);
    network_train(&model, dataset, wavelets, 1, epochs, rate);
    train_hmm(&model, dataset, wavelets, em_iterations);
    model.stage = 2;
    model_save(&model, argv[5]);
  } else {
    printf("##gff-version 3\n");
    unsigned long gene_number = 0;
    for (int record = 0; record < dataset->count; record++) {
      Sequence *sequence = &dataset->sequences[record];
      for (int strand = 0; strand < 2; strand++) {
        char *reverse = strand ? reverse_complement(sequence->bases, sequence->length) : NULL;
        const char *bases = strand ? reverse : sequence->bases;
        fprintf(stderr, "Predict: %s (%c), %d bp\n", sequence->name, strand ? '-' : '+', sequence->length);
        float *probability = network_predict(model.network, wavelets, bases, sequence->length);
        HMMObservation observation = {bases, sequence->length, probability, NULL};
        unsigned char *path = hmm_decode(&model.hmm, &observation);
        hmm_gff(stdout, sequence->name, sequence->length, strand, path, &gene_number);
        free(path);
        free(probability);
        free(reverse);
      }
    }
    if (fflush(stdout) || ferror(stdout)) anno_fail("Cannot write GFF output");
    fprintf(stderr, "Predicted %lu protein-coding candidates\n", gene_number);
  }
  free(wavelets);
  dataset_free(dataset);
  cnn_free(model.network);
  return 0;
}