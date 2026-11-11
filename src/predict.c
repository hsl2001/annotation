#include "include/model.h"
#include "include/utils.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

enum { TILE = 16384, STATES = 24 };

extern int write_isoform(FILE *, const Sequence *, const IntervalList *, char,
                         unsigned long, size_t);

typedef struct {
  int64_t position;
  float probability;
  unsigned kind;
} Point;

typedef struct {
  int64_t stop[3], bad;
} OpeningInfo;

typedef struct {
  Point *points;
  OpeningInfo *openings;
  int count;
} Candidates;

typedef struct {
  int64_t length;
  size_t slot;
  int first, best;
} Node;

typedef struct {
  int node, next;
} Link;

typedef struct {
  Node *nodes;
  Link *edges;
  int *lookup;

  int count, capacity, edge_count, edge_capacity, last;
} Graph;

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
  return x->kind < y->kind ? -1 : x->kind > y->kind;
}

static void prepare_openings(Candidates *c, const Sequence *s) {
  qsort(c->points, c->count, sizeof(*c->points), point_compare);
  unsigned count = 0;
  for (int i = 0; i < c->count; ++i)
    if (!(c->points[i].kind & 1)) {
      if (count >= UINT_MAX / 4)
        bye("too many exon openings for %s", s->name);
      c->points[i].kind = (c->points[i].kind & 3) | count++ << 2;
    }
  c->openings = scalloc(count, sizeof(*c->openings));
  int64_t next[3] = {s->length + 1, s->length + 1, s->length + 1};
  int64_t bad = s->length, p = s->length;
  for (int i = c->count - 1; i >= 0; --i) {
    Point *point = &c->points[i];
    if (point->kind & 1)
      continue;
    for (; p >= point->position; --p) {
      if (p < s->length && base_index(s->bases[p]) < 0)
        bad = p;
      if (p + 3 <= s->length && is_stop(s->bases + p))
        next[p % 3] = p + 3;
    }
    OpeningInfo *info = &c->openings[point->kind >> 2];
    info->bad = bad;
    for (int phase = 0; phase < 3; ++phase)
      info->stop[phase] = next[(point->position + 3 - phase) % 3];
  }
}

static Candidates detect_sites(const Model *model, const Sequence *s) {
  Candidates c = {0};
  int capacity = 0;
  float *power =
      scalloc((TILE + 2 * SITE_RADIUS + 1) * CWT_CHANNELS, sizeof(float));
  float (*probabilities)[SITES] = scalloc(TILE + 1, sizeof(*probabilities));
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
            bye("non-finite point prediction for %s", s->name);
        }
        probabilities[i][type] = z >= 0 ? (float)(1 / (1 + exp(-z))) : 0;
      }
    for (int i = 0; i < rows; ++i)
      for (int type = 0; type < SITES; ++type)
        if (probabilities[i][type] >= 0.5f) {
          c.points = grow(c.points, &capacity, sizeof(*c.points), c.count);
          c.points[c.count++] =
              (Point){start + i, probabilities[i][type], (unsigned)type};
        }
  }
  free(probabilities);
  free(power);
  prepare_openings(&c, s);
  return c;
}

static int graph_link(Graph *g, int from, int point, int state, int64_t width) {
  size_t slot = (size_t)point * STATES + state;
  int to = g->lookup[slot];
  if (!to) {
    g->nodes = grow(g->nodes, &g->capacity, sizeof(*g->nodes), g->count + 1);
    g->lookup[slot] = to = ++g->count;
    g->nodes[to] = (Node){.slot = slot};
    if (point > g->last)
      g->last = point;
  }
  if (!from)
    return to;
  g->edges =
      grow(g->edges, &g->edge_capacity, sizeof(*g->edges), g->edge_count + 1);
  int edge = ++g->edge_count;
  g->edges[edge] = (Link){from, g->nodes[to].first};
  g->nodes[to].first = edge;
  int64_t length = g->nodes[from].length + width;
  if (length > g->nodes[to].length) {
    g->nodes[to].length = length;
    g->nodes[to].best = edge;
  }
  return to;
}

static int phase_of(int state) {
  return state >= 21 ? state - 21 : (state > 0) + (state >= 5);
}

static unsigned tail_of(int state) {
  return state >= 21 ? (state == 23 ? 3U : 0U)
         : !state    ? 0U
                     : (unsigned)(state - (state < 5 ? 1 : 5));
}

static int64_t exon_stop(const Sequence *s, const Point *begin,
                         const OpeningInfo *info, int state) {
  int phase = phase_of(state);
  int64_t stop = info->stop[phase];
  if (phase && begin->position + 3 - phase <= s->length) {
    unsigned code = tail_of(state);
    for (int i = 0; i < 3 - phase; ++i) {
      int b = base_index(s->bases[begin->position + i]);
      if (b < 0)
        return stop;
      code = (code << 2) | (unsigned)b;
    }
    if (code == 48 || code == 50 || code == 56)
      stop = begin->position + 3 - phase;
  }
  return stop;
}

static int extend_exon(const Sequence *s, const Point *begin, const Point *end,
                       int state, int64_t stop) {
  int64_t width = end->position - begin->position;
  if (width < 1)
    return -1;
  int phase = phase_of(state), next = (int)((phase + width) % 3);
  if ((end->kind & 3) == STOP ? end->position != stop || next
                              : end->position >= stop)
    return -1;
  if (state >= 21) {
    for (int i = 0; i < width && phase + i < 3; ++i)
      if (s->bases[begin->position + i] != "ATG"[phase + i])
        return -1;
    if (phase + width < 3)
      return state + (int)width;
  }
  unsigned tail = width >= next ? 0 : tail_of(state);
  for (int64_t p = width >= next ? end->position - next : begin->position;
       p < end->position; ++p)
    tail = (tail << 2) | (unsigned)base_index(s->bases[p]);
  return next ? (next == 1 ? 1 : 5) + (int)tail : 0;
}

static void reverse_intervals(IntervalList *form) {
  for (int i = 0; i < form->count / 2; ++i) {
    Interval swap = form->items[i];
    form->items[i] = form->items[form->count - 1 - i];
    form->items[form->count - 1 - i] = swap;
  }
}

static int64_t covered(const IntervalList *reference, const int64_t *prefix,
                       int64_t p) {
  int lo = 0, hi = reference->count;
  while (lo < hi) {
    int mid = lo + (hi - lo) / 2;
    if (reference->items[mid].start < p)
      lo = mid + 1;
    else
      hi = mid;
  }
  if (!lo)
    return 0;
  Interval e = reference->items[lo - 1];
  return prefix[lo - 1] + (p < e.end ? p : e.end) - e.start;
}

static void emit_pair(const Graph *g, const Point *points, int root,
                      int terminal, const Sequence *genomic, char strand,
                      FILE *output, unsigned long *genes) {
  IntervalList reference = {0}, form = {0};
  for (int node = terminal; node != root;) {
    int from = g->edges[g->nodes[node].best].node;
    Point a = points[g->nodes[from].slot / STATES];
    Point b = points[g->nodes[node].slot / STATES];
    if (!(a.kind & 1))
      push_interval(&reference, a.position, b.position);
    node = from;
  }
  reverse_intervals(&reference);
  int64_t *prefix = scalloc((size_t)reference.count + 1, sizeof(*prefix));
  for (int i = 0; i < reference.count; ++i)
    prefix[i + 1] =
        prefix[i] + reference.items[i].end - reference.items[i].start;
  int64_t length = prefix[reference.count];
  double required = CDS_RECIPROCAL_OVERLAP * length;
  if (length != g->nodes[terminal].length || length < 6 || length % 3)
    bye("invalid longest CDS path for %s", genomic->name);
  unsigned long gene = ++*genes;
  size_t number = 1;
  if (write_isoform(output, genomic, &reference, strand, gene, number) != 1)
    bye("cannot write longest coding form for %s", genomic->name);
  typedef struct {
    int node, edge;
    int64_t overlap;
  } Frame;
  int capacity = 0, depth = 1;
  Frame *stack = grow(NULL, &capacity, sizeof(*stack), 0);
  stack[0] = (Frame){.node = terminal, .edge = g->nodes[terminal].first};
  while (depth) {
    Frame *f = &stack[depth - 1];
    if (f->node == root && f->overlap >= required && f->overlap < length) {
      reverse_intervals(&form);
      if (write_isoform(output, genomic, &form, strand, gene, ++number) != 1)
        bye("cannot write alternative coding form for %s", genomic->name);
      reverse_intervals(&form);
    }
    if (!f->edge) {
      if (!(points[g->nodes[f->node].slot / STATES].kind & 1))
        --form.count;
      --depth;
      continue;
    }
    Link edge = g->edges[f->edge];
    f->edge = edge.next;
    Point a = points[g->nodes[edge.node].slot / STATES];
    Point b = points[g->nodes[f->node].slot / STATES];
    int64_t possible = covered(&reference, prefix, a.position);
    int64_t overlap = f->overlap;
    if (!(a.kind & 1))
      overlap += covered(&reference, prefix, b.position) - possible;
    if (possible > g->nodes[edge.node].length)
      possible = g->nodes[edge.node].length;
    if (overlap + possible < required)
      continue;
    if (!(a.kind & 1))
      push_interval(&form, a.position, b.position);
    stack = grow(stack, &capacity, sizeof(*stack), depth);
    stack[depth++] = (Frame){edge.node, g->nodes[edge.node].first, overlap};
  }
  free(stack);
  free(prefix);
  free(reference.items);
  free(form.items);
}

static void assemble(const Model *model, const Sequence *s,
                     const Sequence *genomic, char strand, const Candidates *c,
                     FILE *output, unsigned long *genes) {
  int *heads = scalloc(c->count, sizeof(*heads));
  Link *junctions = NULL;
  int count = 0, capacity = 0;
  for (int donor = c->count - 1, acceptor = -1; donor >= 0; --donor) {
    if ((c->points[donor].kind & 3) == IE) {
      heads[donor] = acceptor;
      acceptor = donor;
    } else if ((c->points[donor].kind & 3) == IS)
      for (int a = acceptor;
           a >= 0 && c->points[a].position - c->points[donor].position <=
                         model->max_intron;
           a = heads[a])
        if (junction_logit(model, s, c->points[donor].position,
                           c->points[a].position) >= 0) {
          junctions = grow(junctions, &capacity, sizeof(*junctions), count + 1);
          junctions[++count] = (Link){a, heads[donor]};
          heads[donor] = count;
        }
  }
  Graph g = {.lookup = scalloc(smul((size_t)c->count, STATES), sizeof(int))};
  for (int start = 0; start < c->count; ++start) {
    if ((c->points[start].kind & 3) != START)
      continue;
    g.last = start;
    int root = graph_link(&g, 0, start, 21, 0);

    for (int p = start; p <= g.last; ++p) {
      const Point *point = &c->points[p];
      int type = point->kind & 3;
      if (type == STOP) {
        int terminal = g.lookup[(size_t)p * STATES];
        if (terminal)
          emit_pair(&g, c->points, root, terminal, genomic, strand, output,
                    genes);
        continue;
      }
      for (int state = 0; state < STATES; ++state) {
        int from = g.lookup[(size_t)p * STATES + state];
        if (!from)
          continue;
        if (type == IS) {
          for (int link = heads[p]; link; link = junctions[link].next)
            graph_link(&g, from, junctions[link].node, state, 0);
        } else {
          const OpeningInfo *info = &c->openings[point->kind >> 2];
          int64_t stop = exon_stop(s, point, info, state);
          for (int j = p + 1; j < c->count && c->points[j].position <= stop &&
                              c->points[j].position <= info->bad;
               ++j) {
            if (!(c->points[j].kind & 1))
              continue;
            int next = extend_exon(s, point, &c->points[j], state, stop);
            if (next >= 0)
              graph_link(&g, from, j, next,
                         c->points[j].position - point->position);
          }
        }
      }
    }
    for (int n = 1; n <= g.count; ++n)
      g.lookup[g.nodes[n].slot] = 0;
    g.count = g.edge_count = 0;
  }
  free(g.lookup);

  free(g.nodes);
  free(g.edges);
  free(heads);
  free(junctions);
}

static void write_points(FILE *const *files, const Sequence *s,
                         const Candidates *c, int reverse) {
  for (int i = 0; i < c->count; ++i) {
    Point point = c->points[reverse ? c->count - 1 - i : i];
    int type = point.kind & 3;
    if (!files[type])
      continue;
    int64_t p = point.position - (type == STOP);
    if (reverse)
      p = s->length - 1 - p;
    if (p >= 0 && p < s->length &&
        fprintf(files[type], "%s\t%lld\t%lld\t%.6g\n", s->name, (long long)p,
                (long long)p + 1, point.probability) < 0)
      bye("cannot write point bedGraph");
  }
}

unsigned long predict(const Model *model, const Sequence *sequences, int count,
                      FILE *gff, const char *prefix) {
  if (!(CDS_RECIPROCAL_OVERLAP > 0 && CDS_RECIPROCAL_OVERLAP <= 1))
    bye("CDS_RECIPROCAL_OVERLAP must be greater than zero and at most one");
  static const char *const names[] = {"Start", "IS", "IE", "Stop"};
  FILE *files[2 * SITES] = {0};
  for (int i = 0; prefix && i < 2 * SITES; ++i) {
    char path[PATH_BUFFER];
    if (snprintf(path, sizeof(path), "%s.%s.%s.bedgraph", prefix,
                 names[i % SITES],
                 i < SITES ? "plus" : "minus") >= (int)sizeof(path))
      bye("output prefix too long");
    if (!(files[i] = fopen(path, "w")))
      bye("cannot create %s", path);
    fprintf(files[i], "track type=bedGraph name=\"%s %c\"\n", names[i % SITES],
            i < SITES ? '+' : '-');
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
      free(c.openings);
      free(c.points);
    }
    free(reverse);
    cry("%s: %lu Start/Stop pairs", s->name, genes - previous);
  }
  for (int i = 0; i < 2 * SITES; ++i)
    if (files[i] && fclose(files[i]))
      bye("cannot close point bedGraph");
  return genes;
}
