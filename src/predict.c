#include "include/model.h"
#include "include/utils.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>

enum { TILE = 16384, PATHS = 8, SCORE_DROP = 8 };

typedef struct {
  int64_t position;
  float logit;
  unsigned kind;
} Point;

typedef struct {
  Point *points;
  int count;
} Candidates;

typedef struct {
  double score;
  int from, rank;
} Node;

static void *grow(void *items, int *capacity, size_t size, int index) {
  if (index < *capacity)
    return items;
  if (*capacity > INT_MAX / 2)
    bye("point graph allocation exceeds addressable entries");
  *capacity = *capacity ? 2 * *capacity : 256;
  return srealloc(items, (size_t)*capacity, size);
}

static int point_compare(const void *a, const void *b) {
  const Point *x = a, *y = b;
  if (x->position != y->position)
    return x->position < y->position ? -1 : 1;
  // Complete half-open envelopes before opening adjacent transcripts.
  int x_order = x->kind == TTS ? -1 : (int)x->kind;
  int y_order = y->kind == TTS ? -1 : (int)y->kind;
  return x_order < y_order ? -1 : x_order > y_order;
}

static Candidates detect_sites(const Model *model, const Sequence *s) {
  Candidates c = {0};
  int capacity = 0;
  float *power =
      scalloc((TILE + 2 * SITE_RADIUS + 1) * CWT_CHANNELS, sizeof(float));
  float (*logits)[SITES] = scalloc(TILE + 1, sizeof(*logits));
  for (int64_t start = 0; start < s->length; start += TILE) {
    int64_t end = start + TILE < s->length ? start + TILE : s->length;
    int rows = (int)(end - start) + (end == s->length);
    int64_t origin = start > SITE_RADIUS ? start - SITE_RADIUS : 0;
    int64_t limit =
        end + SITE_RADIUS < s->length ? end + SITE_RADIUS : s->length;
    cwt_power(s->bases, s->length, origin, limit, power);
#pragma omp parallel for schedule(static)
    for (int i = 0; i < rows; ++i)
      for (int type = 0; type < SITES; ++type) {
        double z = -INFINITY;
        if (site_candidate(s->bases, s->length, start + i, type)) {
          Sample sample;
          sample_site(s, power, origin, start + i, type, &sample);
          z = site_logit(model, &sample);
          if (!isfinite(z))
            bye("non-finite %s prediction for %s", classifier_names[type],
                s->name);
        }
        logits[i][type] = (float)z;
      }
    for (int i = 0; i < rows; ++i)
      for (int type = 0; type < SITES; ++type)
        if (logits[i][type] >= 0) {
          c.points = grow(c.points, &capacity, sizeof(*c.points), c.count);
          c.points[c.count++] =
              (Point){start + i, logits[i][type], (unsigned)type};
        }
  }
  free(logits);
  free(power);
  if (c.count > 1)
    qsort(c.points, c.count, sizeof(*c.points), point_compare);
  return c;
}

static void reverse_intervals(IntervalList *form) {
  for (int i = 0; i < form->count / 2; ++i) {
    Interval swap = form->items[i];
    form->items[i] = form->items[form->count - 1 - i];
    form->items[form->count - 1 - i] = swap;
  }
}

static double transition_score(const Model *model, const Sequence *s,
                               Point from, Point to, int rescue) {
  int64_t length = to.position - from.position;
  if (length <= 0)
    return -INFINITY;
  if (from.kind == IS && to.kind == IE && length <= model->max_intron) {
    double pair = junction_logit(model, s, from.position, to.position);
    // Weak pair negatives need strong boundary sites and a selected envelope.
    int supported =
        rescue && pair >= -0.5 && from.logit >= log(19) && to.logit >= log(19);
    return pair >= 0 || supported ? to.logit + pair : -INFINITY;
  }
  return (from.kind == TSS || from.kind == IE) &&
                 (to.kind == IS || to.kind == TTS) && length <= model->max_exon
             ? to.logit
             : -INFINITY;
}

static int emit_forms(const Model *model, const Sequence *s,
                      const Sequence *genomic, char strand, const Candidates *c,
                      int first, int last, FILE *output, unsigned long *genes) {
  int count = last - first + 1;
  // A per-point beam bounds exponential splice-path growth in dense graphs.
  Node(*paths)[PATHS] = scalloc(count, sizeof(*paths));
  for (int p = 0; p < count; ++p)
    for (int k = 0; k < PATHS; ++k)
      paths[p][k] = (Node){-INFINITY, -1, 0};
  paths[0][0] = (Node){c->points[first].logit, -1, 0};
  for (int p = 1; p < count; ++p) {
    Point point = c->points[first + p];
    if (point.kind == TSS || (point.kind == TTS && p != count - 1))
      continue;
    int64_t limit = point.kind == IE ? model->max_intron : model->max_exon;
    for (int j = p - 1;
         j >= 0 && point.position - c->points[first + j].position <= limit;
         --j) {
      if (!isfinite(paths[j][0].score))
        continue;
      double edge = transition_score(model, s, c->points[first + j], point, 1);
      if (!isfinite(edge))
        continue;
      for (int k = 0; k < PATHS && isfinite(paths[j][k].score); ++k) {
        double score = paths[j][k].score + edge;
        int rank = PATHS - 1;
        if (score <= paths[p][rank].score)
          continue;
        while (rank && score > paths[p][rank - 1].score) {
          paths[p][rank] = paths[p][rank - 1];
          --rank;
        }
        paths[p][rank] = (Node){score, j, k};
      }
    }
  }
  IntervalList form = {0};
  size_t number = 0;
  for (int k = 0;
       k < PATHS && isfinite(paths[count - 1][k].score) &&
       paths[count - 1][k].score >= paths[count - 1][0].score - SCORE_DROP;
       ++k) {
    form.count = 0;
    int rank = k;
    for (int p = count - 1; p > 0;) {
      Node node = paths[p][rank];
      Point opening = c->points[first + node.from];
      if (opening.kind == TSS || opening.kind == IE)
        push_interval(&form, opening.position, c->points[first + p].position);
      p = node.from;
      rank = node.rank;
    }
    reverse_intervals(&form);
    int status = write_isoform(output, genomic, &form, strand,
                               output ? *genes + 1 : 1, number + 1);
    if (status < 0)
      bye("cannot write transcript for %s", genomic->name);
    number += status;
    if (!output && number)
      break;
  }
  if (output)
    *genes += number != 0;
  free(form.items);
  free(paths);
  return number != 0;
}

static void retain_start(Node *paths, Node node) {
  int rank = 0;
  while (rank < PATHS && paths[rank].from != node.from)
    ++rank;
  if (rank == PATHS)
    rank = PATHS - 1;
  if (node.score <= paths[rank].score)
    return;
  while (rank && node.score > paths[rank - 1].score) {
    paths[rank] = paths[rank - 1];
    --rank;
  }
  paths[rank] = node;
}

static void assemble(const Model *model, const Sequence *s,
                     const Sequence *genomic, char strand, const Candidates *c,
                     FILE *output, unsigned long *genes) {
  // Keep distinct starts, not competing splice paths from the same start.
  // Only PATHS origins survive per point; lower-ranked origins may be pruned.
  // Here from is the TSS index and rank is the preceding completed terminal.
  Node(*nodes)[PATHS] = scalloc(c->count, sizeof(*nodes));
  int best = -1;
  for (int p = 0; p < c->count; ++p) {
    Point point = c->points[p];
    for (int k = 0; k < PATHS; ++k)
      nodes[p][k] = (Node){-INFINITY, -1, -1};
    if (point.kind == TSS) {
      nodes[p][0] =
          (Node){(best < 0 ? 0 : nodes[best][0].score) + point.logit, p, best};
      continue;
    }
    int64_t limit = point.kind == IE ? model->max_intron : model->max_exon;
    for (int j = p - 1;
         j >= 0 && point.position - c->points[j].position <= limit; --j) {
      if (!isfinite(nodes[j][0].score))
        continue;
      double edge = transition_score(model, s, c->points[j], point, 0);
      if (!isfinite(edge))
        continue;
      for (int k = 0; k < PATHS && isfinite(nodes[j][k].score); ++k) {
        Node node = nodes[j][k];
        node.score += edge;
        retain_start(nodes[p], node);
      }
    }
    if (point.kind == TTS)
      for (int k = 0; k < PATHS && isfinite(nodes[p][k].score) &&
                      (best < 0 || nodes[p][k].score > nodes[best][0].score);
           ++k)
        if (emit_forms(model, s, genomic, strand, c, nodes[p][k].from, p, NULL,
                       NULL)) {
          nodes[p][0] = nodes[p][k];
          best = p;
          break;
        }
  }
  IntervalList terminals = {0};
  for (int p = best; p >= 0; p = nodes[p][0].rank)
    push_interval(&terminals, nodes[p][0].from, p);
  for (int i = terminals.count - 1; i >= 0; --i)
    emit_forms(model, s, genomic, strand, c, (int)terminals.items[i].start,
               (int)terminals.items[i].end, output, genes);
  free(terminals.items);
  free(nodes);
}

static void write_points(FILE *const *files, const Sequence *s,
                         const Candidates *c, int reverse) {
  for (int i = 0; i < c->count; ++i) {
    Point point = c->points[reverse ? c->count - 1 - i : i];
    int type = (int)point.kind;
    if (!files[type])
      continue;
    int64_t p = point.position - (type == TTS);
    if (reverse)
      p = s->length - 1 - p;
    if (p >= 0 && p < s->length &&
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
      Candidates c = detect_sites(model, &oriented);
      write_points(files + strand * SITES, s, &c, strand);
      assemble(model, &oriented, s, strand ? '-' : '+', &c, gff, &genes);
      free(c.points);
    }
    free(reverse);
    cry("%s: %lu TSS/TTS transcripts", s->name, genes - previous);
  }
  for (int i = 0; i < 2 * SITES; ++i)
    if (files[i] && fclose(files[i]))
      bye("cannot close %s bedGraph", classifier_names[i % SITES]);
  return genes;
}
