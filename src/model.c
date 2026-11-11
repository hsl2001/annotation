#include "include/model.h"
#include "include/lbfgs.h"
#include "include/utils.h"

#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

const char *const classifier_names[CLASSIFIERS] = {"TSS", "IS", "IE", "TTS",
                                                   "IS-IE pair"};

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

void read_gff3(const char *path, Sequence *sequences, int count) {
  gzFile input = gzopen(path, "rb");
  if (!input)
    bye("cannot open annotation %s", path);
  for (int i = 0; i < count; ++i) {
    free(sequences[i].labels);
    sequences[i].labels = scalloc((size_t)sequences[i].length + 1, 1);
    sequences[i].max_intron = sequences[i].max_exon = 0;
    for (int strand = 0; strand < 2; ++strand) {
      free(sequences[i].introns[strand].items);
      sequences[i].introns[strand] = (IntervalList){0};
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
    if (strcmp(fields[2], "exon") && strcmp(fields[2], "CDS"))
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
      if (used == capacity) {
        capacity = capacity ? capacity * 2 : 256;
        features = srealloc(features, capacity, sizeof(*features));
      }
      features[used++] = (Feature){sstrdup(id),
                                   sequence,
                                   !strcmp(fields[2], "exon"),
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
    int has_exons = 0;
    for (size_t i = first; i < last; ++i)
      has_exons |= features[i].exon;
    for (size_t i = first; i < last; ++i) {
      if (has_exons && !features[i].exon)
        continue;
      Interval span = features[i].span;
      if (exons.count && span.start <= exons.items[exons.count - 1].end) {
        if (span.end > exons.items[exons.count - 1].end)
          exons.items[exons.count - 1].end = span.end;
      } else
        push_interval(&exons, span.start, span.end);
    }
    for (int e = 0; e < exons.count; ++e)
      if (exons.items[e].end - exons.items[e].start > s->max_exon)
        s->max_exon = exons.items[e].end - exons.items[e].start;
    int reverse = f->strand == '-';
    int shift = reverse * SITES;
    int64_t tss =
        reverse ? exons.items[exons.count - 1].end : exons.items[0].start;
    int64_t tts =
        reverse ? exons.items[0].start : exons.items[exons.count - 1].end;
    s->labels[tss] |= 1U << (shift + TSS);
    s->labels[tts] |= 1U << (shift + TTS);
    for (int e = 1; e < exons.count; ++e) {
      int64_t a = exons.items[e - 1].end, b = exons.items[e].start;
      s->labels[reverse ? b : a] |= 1U << (shift + IS);
      s->labels[reverse ? a : b] |= 1U << (shift + IE);
      push_interval(&s->introns[reverse], reverse ? s->length - b : a,
                    reverse ? s->length - a : b);
      if (b - a > s->max_intron)
        s->max_intron = b - a;
    }

    free(exons.items);
  }
  for (size_t i = 0; i < used; ++i)
    free(features[i].parent);
  free(features);
  for (int i = 0; i < count; ++i)
    for (int strand = 0; strand < 2; ++strand)
      unique_intervals(&sequences[i].introns[strand]);
}

int site_candidate(const char *b, int64_t length, int64_t p, int type) {
  if (p < 0 || p > length)
    return 0;
  if (type == TSS || type == TTS)
    return type == TSS ? p < length && base_index(b[p]) >= 0
                       : p > 0 && base_index(b[p - 1]) >= 0;
  if (type == IS)
    return p + 2 <= length &&
           ((b[p] == 'G' && (b[p + 1] == 'T' || b[p + 1] == 'C')) ||
            ((b[p] == 'A' || b[p] == 'T') && b[p + 1] == 'T'));
  return p >= 2 && ((b[p - 2] == 'A' && (b[p - 1] == 'G' || b[p - 1] == 'C')) ||
                    (b[p - 2] == 'G' && b[p - 1] == 'T'));
}

static uint8_t nucleotide(const Sequence *s, int64_t p) {
  int base = p >= 0 && p < s->length ? base_index(s->bases[p]) : -1;
  return (uint8_t)(base < 0 ? NUCLEOTIDES : base);
}

void sample_site(const Sequence *s, const float *power, int64_t origin,
                 int64_t p, int type, Sample *sample) {
  sample->type = (uint8_t)type;
  sample->target = 0;
  sample->weight = 1;
  for (int k = 0; k < SITE_WIDTH; ++k) {
    int64_t at = p + k - SITE_RADIUS;
    sample->bases[k] = nucleotide(s, at);

    double value = 0;
    if (at >= 0 && at < s->length)
      for (int c = 0; c < CWT_CHANNELS; ++c)
        value += power[(at - origin) * CWT_CHANNELS + c];
    sample->power[k] = (float)log1p(value / CWT_CHANNELS);
  }
  double mean = 0;
  for (int k = 0; k < SITE_WIDTH; ++k)
    mean += sample->power[k] / SITE_WIDTH;
  for (int k = 0; k < SITE_WIDTH; ++k)
    sample->power[k] -= (float)mean;
}

void sample_junction(const Sequence *s, int64_t donor, int64_t acceptor,
                     Sample *sample) {
  sample->type = JUNCTION;
  sample->target = 0;
  sample->weight = 1;
  for (int k = 0; k < SITE_WIDTH; ++k) {
    sample->donor[k] = nucleotide(s, donor + k - SITE_RADIUS);
    sample->acceptor[k] = nucleotide(s, acceptor + k - SITE_RADIUS);
  }
  sample->distance[0] = (float)(log1p(acceptor - donor) / 10);
  sample->distance[1] = sample->distance[0] * sample->distance[0];
}

static double sample_logit(const double *weights, const Sample *sample,
                           double *gradient) {
  int indices[4 * SITE_WIDTH + 3], count = 0;
  int junction = sample->type == JUNCTION;
  double values[4 * SITE_WIDTH + 3];
  for (int k = 0; k < SITE_WIDTH; ++k) {
    int at = k * ((junction + 1) * NUCLEOTIDES + 1);
    for (int end = 0; end <= junction; ++end) {
      int a = junction ? (end ? sample->acceptor[k] : sample->donor[k])
                       : sample->bases[k];
      if (a < NUCLEOTIDES) {
        indices[count] = at + end * NUCLEOTIDES + a;
        values[count++] = 1;
      }
    }
    indices[count] = at + (junction + 1) * NUCLEOTIDES;
    values[count++] = junction ? sample->donor[k] < NUCLEOTIDES &&
                                     sample->donor[k] == sample->acceptor[k]
                               : sample->power[k];
  }
  if (!junction)
    for (int k = 0; k + 2 < SITE_WIDTH; ++k)
      if (sample->bases[k] < 4 && sample->bases[k + 1] < 4 &&
          sample->bases[k + 2] < 4) {
        indices[count] = SITE_WIDTH * (NUCLEOTIDES + 1) + k * 64 +
                         sample->bases[k] * 16 + sample->bases[k + 1] * 4 +
                         sample->bases[k + 2];
        values[count++] = 1;
      }
  if (junction)
    for (int k = 0; k < 2; ++k) {
      indices[count] = JUNCTION_FEATURES - 3 + k;
      values[count++] = sample->distance[k];
    }
  indices[count] = (junction ? JUNCTION_FEATURES : SITE_FEATURES) - 1;
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

static int valid_junction(const Sequence *s, int64_t donor, int64_t acceptor,
                          int64_t max_intron) {
  return donor >= 0 && acceptor <= s->length && acceptor >= donor &&
         acceptor - donor >= 20 && acceptor - donor <= max_intron;
}

double junction_logit(const Model *model, const Sequence *s, int64_t donor,
                      int64_t acceptor) {
  if (!valid_junction(s, donor, acceptor, model->max_intron))
    return -INFINITY;
  Sample sample;
  sample_junction(s, donor, acceptor, &sample);
  return sample_logit(model->weights[JUNCTION], &sample, NULL);
}

typedef struct {
  Sample *items;
  size_t count, capacity;
} Samples;

static Sample *append_sample(Samples *samples) {
  if (samples->count == samples->capacity) {
    samples->capacity = samples->capacity ? samples->capacity * 2 : 1024;
    samples->items =
        srealloc(samples->items, samples->capacity, sizeof(*samples->items));
  }
  return &samples->items[samples->count++];
}

static void propose_junction(IntervalList *pairs, const Sequence *s,
                             int64_t donor, int64_t acceptor,
                             int64_t max_intron) {
  if (valid_junction(s, donor, acceptor, max_intron))
    push_interval(pairs, donor, acceptor);
}

static void collect_junctions(Samples *samples, const Sequence *s,
                              const IntervalList *observed,
                              int64_t max_intron) {
  IntervalList pairs = {0};
  for (int i = 0; i < observed->count; ++i) {
    Interval pair = observed->items[i];
    push_interval(&pairs, pair.start, pair.end);
    for (int j = i > 4 ? i - 4 : 0; j < observed->count && j <= i + 4; ++j)
      propose_junction(&pairs, s, pair.start, observed->items[j].end,
                       max_intron);
    for (int offset = -SITE_RADIUS; offset <= SITE_RADIUS; ++offset) {
      if (site_candidate(s->bases, s->length, pair.start + offset, IS))
        propose_junction(&pairs, s, pair.start + offset, pair.end, max_intron);
      if (site_candidate(s->bases, s->length, pair.end + offset, IE))
        propose_junction(&pairs, s, pair.start, pair.end + offset, max_intron);
    }
  }
  unique_intervals(&pairs);
  for (int i = 0; i < pairs.count; ++i) {
    Interval pair = pairs.items[i];
    Sample *sample = append_sample(samples);
    sample_junction(s, pair.start, pair.end, sample);
    sample->target =
        bsearch(&pair, observed->items, observed->count,
                sizeof(*observed->items), compare_interval) != NULL;
  }
  free(pairs.items);
}

static double objective(const double *weights, double *gradient,
                        void *context) {
  Samples *samples = context;
  int type = samples->items[0].type;
  int features = type == JUNCTION ? JUNCTION_FEATURES : SITE_FEATURES;
  memset(gradient, 0, features * sizeof(*gradient));
  double loss = 0;
  for (size_t i = 0; i < samples->count; ++i) {
    const Sample *sample = &samples->items[i];
    double score = sample_logit(weights, sample, gradient);
    loss += sample->weight * (fmax(score, 0) - sample->target * score +
                              log1p(exp(-fabs(score))));
  }
  double lambda = type == JUNCTION ? 1e-4 : 1e-5;
  for (int k = 0; k < features - 1; ++k) {
    loss += 0.5 * lambda * samples->count * weights[k] * weights[k];
    gradient[k] += lambda * samples->count * weights[k];
  }
  return loss;
}

Model *train_model(const Sequence *sequences, int count, int iterations,
                   int threads) {
  Samples samples[CLASSIFIERS] = {0};
  Model *model = scalloc(1, sizeof(*model));
  for (int s = 0; s < count; ++s)
    if (sequences[s].max_intron > model->max_intron)
      model->max_intron = sequences[s].max_intron;
  model->max_intron = (int64_t)fmin(1000000, 5 * fmax(2000, model->max_intron));
  for (int s = 0; s < count; ++s)
    if (sequences[s].max_exon > model->max_exon)
      model->max_exon = sequences[s].max_exon;
  model->max_exon = (int64_t)fmin(1000000, 5 * fmax(1000, model->max_exon));
  for (int s = 0; s < count; ++s) {
    const Sequence *sequence = &sequences[s];
    char *reverse = scalloc((size_t)sequence->length + 1, 1);
    for (int64_t p = 0; p < sequence->length; ++p)
      reverse[p] = complement(sequence->bases[sequence->length - 1 - p]);
    for (int strand = 0; strand < 2; ++strand) {
      Sequence oriented = {.name = sequence->name,
                           .bases = strand ? reverse : sequence->bases,
                           .length = sequence->length};
      collect_junctions(&samples[JUNCTION], &oriented,
                        &sequence->introns[strand], model->max_intron);
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
            if (!target &&
                !site_candidate(oriented.bases, oriented.length, p, type))
              continue;
            Sample *sample = append_sample(&samples[type]);
            sample_site(&oriented, power, origin, p, type, sample);
            sample->target = (uint8_t)target;
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
    for (size_t i = 0; i < samples[type].count; ++i)
      positives += samples[type].items[i].target;
    if (!positives || positives == samples[type].count)
      bye("training annotation lacks positive or negative %s examples",
          classifier_names[type]);
    // Sublinear balancing improves rare-site recall without equalizing priors.
    double positive_weight =
        type == JUNCTION
            ? 1
            : sqrt((double)(samples[type].count - positives) / positives);
    for (size_t i = 0; i < samples[type].count; ++i)
      if (samples[type].items[i].target)
        samples[type].items[i].weight = (float)positive_weight;
    int features = type == JUNCTION ? JUNCTION_FEATURES : SITE_FEATURES;
    model->weights[type][features - 1] =
        log((double)positives / (samples[type].count - positives));
    int completed;
    int status = lbfgs_minimize(features, model->weights[type], objective,
                                &samples[type], iterations, &completed);
    if (status < 0)
      bye("%s optimization failed: status %d", classifier_names[type], status);
    cry("%s: %zu positive, %zu negative, %d iterations", classifier_names[type],
        positives, samples[type].count - positives, completed);
    free(samples[type].items);
  }
  return model;
}

static int model_valid(const Model *model) {
  for (int t = 0; t < CLASSIFIERS; ++t)
    for (int k = 0; k < MODEL_FEATURES; ++k)
      if (!isfinite(model->weights[t][k]))
        return 0;
  return model->max_intron >= 20 && model->max_intron <= 1000000 &&
         model->max_exon >= 1 && model->max_exon <= 1000000;
}

int save_model(const char *path, const Model *model) {
  if (!model_valid(model))
    return EINVAL;
  FILE *output = fopen(path, "wb");
  if (!output)
    return errno;
  uint32_t features[] = {SITE_FEATURES, JUNCTION_FEATURES};
  int ok = fwrite("SURFTX01", 8, 1, output) == 1 &&
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
  int ok = fread(magic, 8, 1, input) == 1 && !memcmp(magic, "SURFTX01", 8) &&
           fread(features, sizeof(features), 1, input) == 1 &&
           features[0] == SITE_FEATURES && features[1] == JUNCTION_FEATURES &&
           fread(&candidate, sizeof(candidate), 1, input) == 1 &&
           model_valid(&candidate) && fgetc(input) == EOF && !ferror(input);
  int closed = fclose(input) == 0;
  if (!ok || !closed)
    return EINVAL;
  *model = candidate;
  return 0;
}
