#include "include/model.h"
#include "include/lbfgs.h"
#include "include/utils.h"

#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

const char *const classifier_names[CLASSIFIERS] = {
    "GS", "IS", "IE", "GE", "IS-IE-pair", "GS-GE-pair"};

typedef struct {
  char *parent;
  int sequence, exon;
  char strand;
  Interval span;
} Feature;

static int feature_compare(const void *a, const void *b) {
  const Feature *x = a, *y = b;
  if (x->sequence != y->sequence)
    return x->sequence < y->sequence ? -1 : 1;
  if (x->strand != y->strand)
    return x->strand < y->strand ? -1 : 1;
  int order = strcmp(x->parent, y->parent);
  return order ? order : compare_interval(&x->span, &y->span);
}

static void unique_intervals(IntervalList *list) {
  if (list->count > 1)
    qsort(list->items, list->count, sizeof(*list->items), compare_interval);
  int count = 0;
  for (int i = 0; i < list->count; ++i)
    if (!count || compare_interval(&list->items[i], &list->items[count - 1]))
      list->items[count++] = list->items[i];
  list->count = count;
}

static void observe_pair(Sequence *s, int strand, int type, int64_t a,
                         int64_t b) {
  int start = type == JUNCTION ? IS : GS, end = type == JUNCTION ? IE : GE;
  s->labels[strand ? b : a] |= 1U << (strand * SITES + start);
  s->labels[strand ? a : b] |= 1U << (strand * SITES + end);
  IntervalList *pairs = type == JUNCTION ? s->introns : s->spans;
  push_interval(&pairs[strand], strand ? s->length - b : a,
                strand ? s->length - a : b);
  int64_t *maximum = type == JUNCTION ? &s->max_intron : &s->max_span;
  *maximum = b - a > *maximum ? b - a : *maximum;
}

void read_gff3(const char *path, Sequence *sequences, int count) {
  gzFile input = gzopen(path, "rb");
  if (!input)
    bye("cannot open annotation %s", path);
  for (int i = 0; i < count; ++i) {
    Sequence *s = &sequences[i];
    free(s->labels);
    s->labels = scalloc((size_t)s->length + 1, 1);
    s->max_intron = s->max_span = 0;
    for (int strand = 0; strand < 2; ++strand) {
      free(s->introns[strand].items);
      free(s->spans[strand].items);
      s->introns[strand] = s->spans[strand] = (IntervalList){0};
    }
  }
  Feature *features = NULL;
  size_t used = 0, capacity = 0;
  char line[65536];
  long line_number = 0;
  while (gzgets(input, line, sizeof(line))) {
    ++line_number;
    if (!strncmp(line, "##FASTA", 7))
      break;
    if (*line == '#' || *line == '\n' || *line == '\r')
      continue;
    char *fields[9], *cursor = line;
    int columns = 0;
    while (columns < 9 && cursor) {
      fields[columns++] = cursor;
      cursor = strchr(cursor, '\t');
      if (cursor)
        *cursor++ = '\0';
    }
    if (columns != 9 || cursor)
      bye("%s:%ld: expected 9 columns", path, line_number);
    if (strcmp(fields[2], "exon") && strcmp(fields[2], "CDS") &&
        strcmp(fields[2], "five_prime_UTR") &&
        strcmp(fields[2], "three_prime_UTR"))
      continue;
    int sequence = -1;
    for (int i = 0; i < count; ++i)
      if (!strcmp(fields[0], sequences[i].name))
        sequence = i;
    if (sequence < 0)
      bye("%s:%ld: sequence absent from FASTA", path, line_number);
    char *end;
    errno = 0;
    int64_t a = strtoll(fields[3], &end, 10);
    if (errno || end == fields[3] || *end || a < 1)
      bye("%s:%ld: invalid exon/CDS start", path, line_number);
    int64_t b = strtoll(fields[4], &end, 10);
    char strand = fields[6][0];
    if (errno || end == fields[4] || *end || b < a ||
        b > sequences[sequence].length || fields[6][1] ||
        (strand != '+' && strand != '-'))
      bye("%s:%ld: invalid exon/CDS end or strand", path, line_number);
    char *parent = NULL, *save;
    fields[8][strcspn(fields[8], "\r\n")] = 0;
    for (char *field = strtok_r(fields[8], ";", &save); field;
         field = strtok_r(NULL, ";", &save))
      if (!strncmp(field, "Parent=", 7))
        parent = field + 7;
    if (!parent || !*parent)
      bye("%s:%ld: exon/CDS without Parent", path, line_number);
    for (char *id = strtok_r(parent, ",", &save); id;
         id = strtok_r(NULL, ",", &save)) {
      int type = fields[2][0];
      if (used == capacity) {
        capacity = capacity ? capacity * 2 : 256;
        features = srealloc(features, capacity, sizeof(*features));
      }
      features[used++] =
          (Feature){sstrdup(id),
                    sequence,
                    (type == 'e') | (type == 'f') * 2 | (type == 't') * 4,
                    strand,
                    {a - 1, b}};
    }
  }
  if (gzclose(input) != Z_OK)
    bye("cannot read annotation %s", path);
  if (!used)
    bye("no exon or CDS features in %s", path);
  qsort(features, used, sizeof(*features), feature_compare);
  for (size_t first = 0, last; first < used; first = last) {
    Feature *f = &features[first];
    for (last = first + 1;
         last < used && features[last].sequence == f->sequence &&
         features[last].strand == f->strand &&
         !strcmp(features[last].parent, f->parent);
         ++last)
      ;
    Sequence *s = &sequences[f->sequence];
    IntervalList exons = {0};
    int flags = 0;
    for (size_t i = first; i < last; ++i)
      flags |= features[i].exon;
    for (size_t i = first; i < last; ++i) {
      if (flags & 1 ? !(features[i].exon & 1) : features[i].exon & 6)
        continue;
      Interval span = features[i].span;
      if (exons.count && span.start <= exons.items[exons.count - 1].end) {
        if (span.end > exons.items[exons.count - 1].end)
          exons.items[exons.count - 1].end = span.end;
      } else
        push_interval(&exons, span.start, span.end);
    }
    int reverse = f->strand == '-';
    observe_pair(s, reverse, SPAN, exons.items[0].start,
                 exons.items[exons.count - 1].end);
    for (int e = 1; e < exons.count; ++e)
      observe_pair(s, reverse, JUNCTION, exons.items[e - 1].end,
                   exons.items[e].start);

    free(exons.items);
  }
  for (size_t i = 0; i < used; ++i)
    free(features[i].parent);
  free(features);
  for (int i = 0; i < count; ++i)
    for (int strand = 0; strand < 2; ++strand) {
      unique_intervals(&sequences[i].introns[strand]);
      unique_intervals(&sequences[i].spans[strand]);
    }
}

static uint8_t nucleotide(const Sequence *s, int64_t p) {
  int base = p >= 0 && p < s->length ? base_index(s->bases[p]) : -1;
  return (uint8_t)(base < 0 ? NUCLEOTIDES : base);
}

void sample_site(const Sequence *s, const float *power, int64_t origin,
                 int64_t p, int type, Sample *sample) {
  *sample = (Sample){.type = (uint8_t)type, .weight = 1};
  double mean = 0;
  for (int k = 0; k < SITE_WIDTH; ++k) {
    int64_t at = p + k - SITE_RADIUS;
    sample->bases[k] = nucleotide(s, at);

    double value = 0;
    if (at >= 0 && at < s->length)
      for (int c = 0; c < CWT_CHANNELS; ++c)
        value += power[(at - origin) * CWT_CHANNELS + c];
    sample->power[k] = (float)log1p(value / CWT_CHANNELS);
    mean += sample->power[k] / SITE_WIDTH;
  }
  for (int k = 0; k < SITE_WIDTH; ++k)
    sample->power[k] -= (float)mean;
}

void sample_pair(const Sequence *s, int64_t start, int64_t end, int type,
                 Sample *sample) {
  *sample = (Sample){.type = (uint8_t)type, .weight = 1};
  for (int side = 0; side < 2; ++side)
    for (int k = 0; k < SITE_WIDTH; ++k)
      sample->ends[side][k] =
          nucleotide(s, (side ? end : start) + k - SITE_RADIUS);
  sample->distance[0] = (float)(log1p(end - start));
  sample->distance[1] = sample->distance[0] * sample->distance[0];
}

static double sample_logit(const double *weights, const Sample *sample,
                           double *gradient) {
  int indices[4 * SITE_WIDTH + 3], count = 0;
  int pair = sample->type >= SITES;
  double values[4 * SITE_WIDTH + 3];
  for (int k = 0; k < SITE_WIDTH; ++k) {
    int at = k * (pair ? 2 * NUCLEOTIDES : NUCLEOTIDES + 1);
    for (int end = 0; end <= pair; ++end) {
      int a = pair ? sample->ends[end][k] : sample->bases[k];
      if (a < NUCLEOTIDES) {
        indices[count] = at + end * NUCLEOTIDES + a;
        values[count++] = 1;
      }
    }
    if (!pair) {
      indices[count] = at + NUCLEOTIDES;
      values[count++] = sample->power[k];
    }
  }

  if (pair)
    for (int k = 0; k < 2; ++k) {
      indices[count] = PAIR_FEATURES - 3 + k;
      values[count++] = sample->distance[k];
    }
  indices[count] = (pair ? PAIR_FEATURES : SITE_FEATURES) - 1;
  values[count++] = 1;
  double score = 0;
  for (int k = 0; k < count; ++k)
    score += weights[indices[k]] * values[k];
  if (gradient) {
    double residual =
        (score >= 0 ? 1 / (1 + exp(-score)) : exp(score) / (1 + exp(score))) -
        sample->target;
    for (int k = 0; k < count; ++k)
      gradient[indices[k]] += sample->weight * residual * values[k];
  }
  return score;
}

double site_logit(const Model *model, const Sample *sample) {
  return sample_logit(model->weights[sample->type], sample, NULL);
}

static int valid_pair(const Sequence *s, int64_t start, int64_t end, int type,
                      int64_t maximum) {
  return (type == JUNCTION || type == SPAN) && start >= 0 && end <= s->length &&
         end >= start && end - start >= (type == JUNCTION ? 20 : 1) &&
         end - start <= maximum;
}

double pair_logit(const Model *model, const Sequence *s, int64_t start,
                  int64_t end, int type) {
  int64_t maximum = type == JUNCTION ? model->max_intron : model->max_span;
  if (!valid_pair(s, start, end, type, maximum))
    return -INFINITY;
  Sample sample;
  sample_pair(s, start, end, type, &sample);
  return site_logit(model, &sample);
}

typedef struct {
  Sample *items;
  size_t count, capacity;
  double population;
} Samples;

static Sample *append_sample(Samples *samples) {
  if (samples->count == samples->capacity) {
    samples->capacity = samples->capacity ? samples->capacity * 2 : 1024;
    samples->items =
        srealloc(samples->items, samples->capacity, sizeof(*samples->items));
  }
  return &samples->items[samples->count++];
}

static void collect_pairs(Samples *samples, const Sequence *s,
                          const IntervalList *observed, int type,
                          int64_t maximum) {
  IntervalList pairs = {0};
  for (int i = 0; i < observed->count; ++i) {
    Interval pair = observed->items[i];
    for (int j = i > 4 ? i - 4 : 0; j < observed->count && j <= i + 4; ++j)
      push_interval(&pairs, pair.start, observed->items[j].end);
    for (int offset = -SITE_RADIUS; offset <= SITE_RADIUS; ++offset) {
      push_interval(&pairs, pair.start + offset, pair.end);
      push_interval(&pairs, pair.start, pair.end + offset);
    }
  }
  unique_intervals(&pairs);
  for (int i = 0; i < pairs.count; ++i) {
    Interval pair = pairs.items[i];
    if (!valid_pair(s, pair.start, pair.end, type, maximum))
      continue;
    Sample *sample = append_sample(samples);
    sample_pair(s, pair.start, pair.end, type, sample);
    // Only an observed endpoint pair is positive, never its recombinations.
    sample->target = bsearch(&pair, observed->items, observed->count,
                             sizeof(pair), compare_interval) != NULL;
  }
  free(pairs.items);
}

static double objective(const double *weights, double *gradient,
                        void *context) {
  Samples *samples = context;
  int type = samples->items[0].type;
  int features = type < SITES ? SITE_FEATURES : PAIR_FEATURES;
  memset(gradient, 0, features * sizeof(*gradient));
  double loss = 0;
  for (size_t i = 0; i < samples->count; ++i) {
    const Sample *sample = &samples->items[i];
    double score = sample_logit(weights, sample, gradient);
    loss += sample->weight * (fmax(score, 0) - sample->target * score +
                              log1p(exp(-fabs(score))));
  }
  // L2 regularization to prevent overfitting.
  double lambda = 1e-5;
  for (int k = 0; k < features - 1; ++k) {
    loss += 0.5 * lambda * samples->population * weights[k] * weights[k];
    gradient[k] += lambda * samples->population * weights[k];
  }
  return loss;
}

Model *train_model(const Sequence *sequences, int count, int iterations,
                   int threads) {
  Samples samples[CLASSIFIERS] = {0};
  int stride = 1;
  double length = 0;
  for (int s = 0; s < count; ++s)
    length += sequences[s].length;
  // Bound negative storage, preserving its population weight and all positives.
  while (2 * length / stride > STRIDE)
    stride *= 2;
  cry("site negative sampling stride: %d; all positives retained", stride);
  Model *model = scalloc(1, sizeof(*model));
  for (int s = 0; s < count; ++s) {
    model->max_intron = fmax(model->max_intron, sequences[s].max_intron);
    model->max_span = fmax(model->max_span, sequences[s].max_span);
  }
  model->max_intron = (int64_t)fmin(STRIDE, fmax(2000, model->max_intron));
  model->max_span = (int64_t)fmin(STRIDE, fmax(2000, model->max_span));
  for (int s = 0; s < count; ++s) {
    const Sequence *sequence = &sequences[s];
    char *reverse = scalloc((size_t)sequence->length + 1, 1);
    for (int64_t p = 0; p < sequence->length; ++p)
      reverse[p] = complement(sequence->bases[sequence->length - 1 - p]);
    for (int strand = 0; strand < 2; ++strand) {
      Sequence oriented = {.name = sequence->name,
                           .bases = strand ? reverse : sequence->bases,
                           .length = sequence->length};
      for (int type = JUNCTION; type < CLASSIFIERS; ++type)
        collect_pairs(
            &samples[type], &oriented,
            &(type == JUNCTION ? sequence->introns : sequence->spans)[strand],
            type, type == JUNCTION ? model->max_intron : model->max_span);
      enum { TILE = 16384 };
      float *power =
          scalloc((TILE + 2 * SITE_RADIUS) * CWT_CHANNELS, sizeof(float));
      for (int64_t begin = 0; begin <= oriented.length; begin += TILE) {
        int64_t end = begin + TILE;
        if (end > oriented.length)
          end = oriented.length;
        int64_t origin = begin > SITE_RADIUS ? begin - SITE_RADIUS : 0;
        int64_t limit = end + SITE_RADIUS;
        if (limit > oriented.length)
          limit = oriented.length;
        cwt_power(oriented.bases, oriented.length, origin, limit, power);
        for (int64_t p = begin; p < end + (end == oriented.length); ++p) {
          unsigned label = sequence->labels[strand ? oriented.length - p : p] >>
                           (strand * SITES);
          for (int type = 0; type < SITES; ++type) {
            int target = !!(label & (1U << type));
            uint64_t key = (uint64_t)p + (uint64_t)s * 0x9e3779b97f4a7c15ULL +
                           (unsigned)(strand * SITES + type);
            key = (key ^ (key >> 30)) * 0xbf58476d1ce4e5b9ULL;
            key = (key ^ (key >> 27)) * 0x94d049bb133111ebULL;
            key ^= key >> 31;
            if (!target && (key & (unsigned)(stride - 1)))
              continue;
            Sample *sample = append_sample(&samples[type]);
            sample_site(&oriented, power, origin, p, type, sample);
            sample->target = (uint8_t)target;
            sample->weight = target ? 1 : (float)stride;
          }
        }
        if (end == oriented.length)
          break;
      }
      free(power);
    }
    free(reverse);
  }

#pragma omp parallel for num_threads(                                          \
        threads < CLASSIFIERS ? threads : CLASSIFIERS) schedule(static)
  for (int type = 0; type < CLASSIFIERS; ++type) {
    size_t positives = 0;
    for (size_t i = 0; i < samples[type].count; ++i) {
      positives += samples[type].items[i].target;
      samples[type].population += samples[type].items[i].weight;
    }
    double negatives = samples[type].population - positives;
    if (!positives || positives == samples[type].count)
      bye("training annotation lacks positive or negative %s examples",
          classifier_names[type]);
    // Equalize effective positive and negative weight for site classifiers.
    double positive_weight = type < SITES ? negatives / positives : 1;
    for (size_t i = 0; i < samples[type].count; ++i)
      if (samples[type].items[i].target)
        samples[type].items[i].weight = (float)positive_weight;
    int features = type < SITES ? SITE_FEATURES : PAIR_FEATURES;
    model->weights[type][features - 1] =
        type < SITES ? 0 : log((double)positives / negatives);
    int completed;
    int status = lbfgs_minimize(features, model->weights[type], objective,
                                &samples[type], iterations, &completed);
    if (status < 0)
      bye("%s optimization failed: status %d", classifier_names[type], status);
    cry("%s: %zu positive, %zu sampled negative (weighted %.0f), %d iterations",
        classifier_names[type], positives, samples[type].count - positives,
        negatives, completed);
    free(samples[type].items);
  }
  return model;
}

int save_model(const char *path, const Model *model) {
  FILE *output = fopen(path, "wb");
  if (!output)
    return errno;
  uint32_t features[] = {SITE_FEATURES, PAIR_FEATURES};
  int ok = fwrite("SURFFING", 8, 1, output) == 1 &&
           fwrite(features, sizeof(features), 1, output) == 1 &&
           fwrite(model, sizeof(*model), 1, output) == 1;
  int closed = fclose(output) == 0;
  return ok && closed ? 0 : EIO;
}

int load_model(const char *path, Model *model) {
  FILE *input = fopen(path, "rb");
  if (!input)
    return errno;
  char magic[8];
  uint32_t features[2];
  Model candidate;
  int ok = fread(magic, 8, 1, input) == 1 && !memcmp(magic, "SURFFING", 8) &&
           fread(features, sizeof(features), 1, input) == 1 &&
           features[0] == SITE_FEATURES && features[1] == PAIR_FEATURES &&
           fread(&candidate, sizeof(candidate), 1, input) == 1 &&
           fgetc(input) == EOF && !ferror(input);
  int closed = fclose(input) == 0;
  if (!ok || !closed)
    return EINVAL;
  *model = candidate;
  return 0;
}
