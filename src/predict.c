#include "include/model.h"
#include "include/utils.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

enum { TILE = 16384 };

typedef struct {
  int64_t position;
  float logit;
  unsigned kind;
} Point;

typedef struct {
  Point *points;
  int count;
} Candidates;

static Candidates detect_sites(const Model *model, const Sequence *s,
                               int64_t begin, int64_t end, unsigned mask) {
  Candidates c = {0};
  int capacity = 0;
  float *power =
      scalloc((TILE + 2 * SITE_RADIUS) * CWT_CHANNELS, sizeof(float));
  float (*logits)[SITES] = scalloc(TILE, sizeof(*logits));
  for (; begin <= end; begin += TILE) {
    int64_t stop = begin + TILE < end + 1 ? begin + TILE : end + 1;
    int64_t origin = begin > SITE_RADIUS ? begin - SITE_RADIUS : 0;
    int64_t limit =
        stop + SITE_RADIUS < s->length ? stop + SITE_RADIUS : s->length;
    cwt_power(s->bases, s->length, origin, limit, power);
#pragma omp parallel for schedule(static)
    for (int i = 0; i < stop - begin; ++i)
      for (int type = 0; type < SITES; ++type) {
        double z = -INFINITY;
        if (mask & (1U << type)) {
          Sample sample;
          sample_site(s, power, origin, begin + i, type, &sample);
          z = site_logit(model, &sample);
          if (!isfinite(z))
            bye("non-finite %s prediction for %s", classifier_names[type],
                s->name);
        }
        logits[i][type] = (float)z;
      }
    for (int i = 0; i < stop - begin; ++i)
      // GE precedes GS at shared boundaries, allowing adjacent spans.
      for (int type = SITES - 1; type >= 0; --type)
        if (logits[i][type] >= 0) {
          if (c.count == capacity) {
            if (capacity > INT_MAX / 2)
              bye("point graph allocation exceeds addressable entries");
            capacity = capacity ? capacity * 2 : 256;
            c.points = srealloc(c.points, capacity, sizeof(*c.points));
          }
          c.points[c.count++] =
              (Point){begin + i, logits[i][type], (unsigned)type};
        }
  }
  free(logits);
  free(power);
  return c;
}

static IntervalList find_spans(const Model *model, const Sequence *s,
                               const Candidates *c) {
  typedef struct {
    double score;
    int start, previous;
  } Node;
  Node *nodes = scalloc(c->count, sizeof(*nodes));
  int best = -1;
  for (int p = 0; p < c->count; ++p) {
    Point end = c->points[p];
    nodes[p] = (Node){best < 0 ? 0 : nodes[best].score, -1, best};
    if (end.kind != GE)
      continue;
    for (int j = p - 1;
         j >= 0 && end.position - c->points[j].position <= model->max_span;
         --j) {
      Point start = c->points[j];
      if (start.kind != GS)
        continue;
      double pair = pair_logit(model, s, start.position, end.position, SPAN);
      double score = nodes[j].score + start.logit + end.logit + pair;
      if (isfinite(pair) && pair >= 0 && score > nodes[p].score)
        nodes[p] = (Node){score, j, nodes[j].previous};
    }
    if (nodes[p].start >= 0)
      best = p;
  }
  IntervalList spans = {0};
  for (int p = best; p >= 0; p = nodes[p].previous)
    push_interval(&spans, c->points[nodes[p].start].position,
                  c->points[p].position);
  free(nodes);
  return spans;
}

typedef struct {
  Interval span;
  double score;
} Intron;

typedef struct {
  Intron *items;
  int count, capacity;
} Introns;

static void push_intron(Introns *list, Interval span, double score) {
  if (list->count == list->capacity) {
    list->capacity = list->capacity ? list->capacity * 2 : 8;
    list->items =
        srealloc(list->items, (size_t)list->capacity, sizeof(*list->items));
  }
  list->items[list->count++] = (Intron){span, score};
}

static Introns find_introns(const Model *model, const Sequence *s,
                            const Candidates *c) {
  Introns introns = {0};
  for (int i = 0; i < c->count; ++i)
    if (c->points[i].kind == IS)
      for (int j = i + 1;
           j < c->count &&
           c->points[j].position - c->points[i].position <= model->max_intron;
           ++j)
        if (c->points[j].kind == IE) {
          double pair = pair_logit(model, s, c->points[i].position,
                                   c->points[j].position, JUNCTION);
          double score = c->points[i].logit + c->points[j].logit + pair;
          if (score > 0)
            push_intron(
                &introns,
                (Interval){c->points[i].position, c->points[j].position},
                score);
        }
  return introns;
}

static IntervalList refine_chain(const Sequence *s, char strand,
                                 IntervalList exons) {
  for (;;) {
    Interval best = longest_orf(s, &exons, strand);
    int64_t length = best.end - best.start;
    int remove = -1;
    for (int i = 0; i + 1 < exons.count; ++i) {
      IntervalList merged = {0};
      for (int j = 0; j < exons.count; ++j) {
        if (j == i) {
          push_interval(&merged, exons.items[j].start, exons.items[j + 1].end);
          ++j;
        } else
          push_interval(&merged, exons.items[j].start, exons.items[j].end);
      }
      Interval orf = longest_orf(s, &merged, strand);
      if (orf.end - orf.start > length) {
        length = orf.end - orf.start;
        remove = i;
      }
      free(merged.items);
    }
    if (remove < 0)
      break;
    exons.items[remove].end = exons.items[remove + 1].end;
    memmove(&exons.items[remove + 1], &exons.items[remove + 2],
            (size_t)(exons.count - remove - 2) * sizeof(*exons.items));
    --exons.count;
  }
  return exons;
}

static IntervalList best_chain(Interval span, const Introns *introns) {
  int count = introns->count;
  double *score = scalloc((size_t)count, sizeof(*score));
  int *previous = scalloc((size_t)count, sizeof(*previous));
  int best = -1;
  for (int j = 0; j < count; ++j) {
    double value = introns->items[j].score;
    score[j] = value;
    previous[j] = -1;
    for (int i = 0; i < j; ++i)
      if (introns->items[i].span.end < introns->items[j].span.start &&
          score[i] + value > score[j]) {
        score[j] = score[i] + value;
        previous[j] = i;
      }
    if (best < 0 || score[j] > score[best])
      best = j;
  }
  IntervalList exons = {0};
  int *chain = scalloc((size_t)count, sizeof(*chain));
  int length = 0;
  for (int j = best; j >= 0; j = previous[j])
    chain[length++] = j;
  int64_t start = span.start;
  for (int k = length - 1; k >= 0; --k) {
    Interval intron = introns->items[chain[k]].span;
    push_interval(&exons, start, intron.start);
    start = intron.end;
  }
  push_interval(&exons, start, span.end);
  free(chain);
  free(previous);
  free(score);
  return exons;
}

static int emit_forms(FILE *output, const Sequence *s, char strand,
                      Interval span, const Introns *introns,
                      unsigned long gene) {
  IntervalList exons = best_chain(span, introns);
  exons = refine_chain(s, strand, exons);
  int status = write_isoform(output, s, &exons, strand, gene, 1);
  free(exons.items);
  if (status < 0)
    bye("cannot write transcript for %s", s->name);
  return status > 0;
}

static void write_points(FILE *const *files, const Sequence *s,
                         const Candidates *c, int reverse) {
  for (int i = 0; i < c->count; ++i) {
    Point point = c->points[reverse ? c->count - 1 - i : i];
    int type = (int)point.kind;
    int64_t p = point.position - (type == GE);
    if (reverse)
      p = s->length - 1 - p;
    if (files[type] && p >= 0 && p < s->length &&
        fprintf(files[type], "%s\t%lld\t%lld\t%.6g\n", s->name, (long long)p,
                (long long)p + 1, 1 / (1 + exp(-point.logit))) < 0)
      bye("cannot write %s bedGraph", classifier_names[type]);
  }
}

unsigned long predict(const Model *model, const Sequence *sequences, int count,
                      FILE *gff, const char *prefix) {
  FILE *files[2 * SITES] = {0};
  for (int i = 0; prefix && i < 2 * SITES; ++i) {
    char path[PATH_BUFFER];
    if (snprintf(path, sizeof(path), "%s.%s.%s.bedgraph", prefix,
                 classifier_names[i % SITES],
                 i < SITES ? "plus" : "minus") >= (int)sizeof(path))
      bye("output prefix too long");
    if (!(files[i] = fopen(path, "w")))
      bye("cannot create %s", path);
    fprintf(files[i], "track type=bedGraph name=\"%s %c\"\n",
            classifier_names[i % SITES], i < SITES ? '+' : '-');
  }
  unsigned long genes = 0;
  for (int i = 0; i < count; ++i) {
    const Sequence *s = &sequences[i];
    char *reverse = scalloc((size_t)s->length + 1, 1);
    for (int64_t p = 0; p < s->length; ++p)
      reverse[p] = complement(s->bases[s->length - 1 - p]);
    unsigned long previous = genes;
    for (int strand = 0; strand < 2; ++strand) {
      Sequence oriented = *s;
      oriented.bases = strand ? reverse : s->bases;
      Candidates ends =
          detect_sites(model, &oriented, 0, s->length, (1U << GS) | (1U << GE));
      write_points(files + strand * SITES, s, &ends, strand);
      IntervalList spans = find_spans(model, &oriented, &ends);
      free(ends.points);
      for (int j = 0; j < spans.count; ++j) {
        Interval span = spans.items[strand ? j : spans.count - 1 - j];
        Candidates sites = detect_sites(model, &oriented, span.start + 1,
                                        span.end - 1, (1U << IS) | (1U << IE));
        write_points(files + strand * SITES, s, &sites, strand);
        Introns introns = find_introns(model, &oriented, &sites);
        genes +=
            emit_forms(gff, s, strand ? '-' : '+', span, &introns, genes + 1);
        free(introns.items);
        free(sites.points);
      }
      free(spans.items);
    }
    free(reverse);
    cry("%s: %lu transcripts", s->name, genes - previous);
  }
  for (int i = 0; i < 2 * SITES; ++i)
    if (files[i] && fclose(files[i]))
      bye("cannot close %s bedGraph", classifier_names[i % SITES]);
  return genes;
}
