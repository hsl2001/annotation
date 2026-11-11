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
  int count, capacity, opening_count;
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

static void *grow(void *items, int *capacity, size_t size) {
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
  c->opening_count = 0;
  for (int i = 0; i < c->count; ++i)
    if ((c->points[i].kind & 3) == START || (c->points[i].kind & 3) == IE) {
      if ((unsigned)c->opening_count >= UINT_MAX / 4)
        bye("too many exon openings for %s", s->name);
      c->points[i].kind = (c->points[i].kind & 3) | (unsigned)c->opening_count++
                                                        << 2;
    }
  c->openings = scalloc(c->opening_count, sizeof(*c->openings));
  int64_t next[3] = {s->length + 1, s->length + 1, s->length + 1};
  int64_t bad = s->length;
  int at = c->count - 1;
  for (int64_t p = s->length; p >= 0 && at >= 0; --p) {
    if (p < s->length) {
      if (base_index(s->bases[p]) < 0)
        bad = p;
      if (p + 3 <= s->length && is_stop(s->bases + p))
        next[p % 3] = p + 3;
    }
    while (at >= 0 && c->points[at].position == p) {
      Point *point = &c->points[at--];
      int type = point->kind & 3;
      if (type != START && type != IE)
        continue;
      OpeningInfo *info = &c->openings[point->kind >> 2];
      info->bad = bad;
      for (int phase = 0; phase < 3; ++phase)
        info->stop[phase] = next[(p + (3 - phase) % 3) % 3];
    }
  }
}

static Candidates detect_sites(const Model *model, const Sequence *s) {
  Candidates c = {0};
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
        double probability = 0;
        if (site_candidate(s->bases, s->length, start + i, type)) {
          Sample sample;
          sample_site(s, power, origin, start + i, type, &sample);
          double z = site_logit(model, &sample);
          if (!isfinite(z))
            bye("non-finite point prediction for %s", s->name);
          if (z >= 0)
            probability = 1 / (1 + exp(-z));
        }
        probabilities[i][type] = (float)probability;
      }
    for (int i = 0; i < rows; ++i)
      for (int type = 0; type < SITES; ++type)
        if (probabilities[i][type] >= 0.5f) {
          if (c.count == c.capacity)
            c.points = grow(c.points, &c.capacity, sizeof(*c.points));
          c.points[c.count++] =
              (Point){start + i, probabilities[i][type], (unsigned)type};
        }
  }
  free(probabilities);
  free(power);
  prepare_openings(&c, s);
  return c;
}

static int graph_node(Graph *g, int point, int state) {
  size_t slot = (size_t)point * STATES + state;
  if (!g->lookup[slot]) {
    if (g->count + 1 >= g->capacity)
      g->nodes = grow(g->nodes, &g->capacity, sizeof(*g->nodes));
    int id = ++g->count;
    g->lookup[slot] = id;
    g->nodes[id] = (Node){.length = -1, .slot = slot};
    if (point > g->last)
      g->last = point;
  }
  return g->lookup[slot];
}

static void graph_link(Graph *g, const Candidates *c, int from, int point,
                       int state) {
  int to = graph_node(g, point, state);
  int previous = (int)(g->nodes[from].slot / STATES);
  int64_t width =
      (c->points[previous].kind & 3) == IS
          ? 0
          : c->points[point].position - c->points[previous].position;
  if (g->edge_count + 1 >= g->edge_capacity)
    g->edges = grow(g->edges, &g->edge_capacity, sizeof(*g->edges));
  int edge = ++g->edge_count;
  g->edges[edge] = (Link){from, g->nodes[to].first};
  g->nodes[to].first = edge;
  int64_t length = g->nodes[from].length + width;
  if (length > g->nodes[to].length) {
    g->nodes[to].length = length;
    g->nodes[to].best = edge;
  }
}

static int phase_of(int state) {
  return state >= 21 ? state - 21 : !state ? 0 : state < 5 ? 1 : 2;
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
      return 21 + phase + (int)width;
  }
  unsigned tail = tail_of(state);
  if (width >= next)
    tail = 0;
  for (int64_t p = width >= next ? end->position - next : begin->position;
       p < end->position; ++p)
    tail = (tail << 2) | (unsigned)base_index(s->bases[p]);
  return next ? (next == 1 ? 1 : 5) + (int)(tail & ((1U << (2 * next)) - 1))
              : 0;
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

static void emit_pair(const Graph *g, const Candidates *c, int root,
                      int terminal, const Sequence *genomic, char strand,
                      FILE *output, unsigned long *genes) {
  IntervalList reference = {0}, form = {0};
  for (int node = terminal; node != root;) {
    int from = g->edges[g->nodes[node].best].node;
    Point a = c->points[g->nodes[from].slot / STATES];
    Point b = c->points[g->nodes[node].slot / STATES];
    if ((a.kind & 3) != IS)
      push_interval(&reference, a.position, b.position);
    node = from;
  }
  reverse_intervals(&reference);
  int64_t *prefix = scalloc((size_t)reference.count + 1, sizeof(*prefix));
  for (int i = 0; i < reference.count; ++i)
    prefix[i + 1] =
        prefix[i] + reference.items[i].end - reference.items[i].start;
  int64_t length = prefix[reference.count];
  if (length != g->nodes[terminal].length || length < 6 || length % 3)
    bye("invalid longest CDS path for %s", genomic->name);
  unsigned long gene = ++*genes;
  size_t number = 1;
  if (write_isoform(output, genomic, &reference, strand, gene, number) != 1)
    bye("cannot write longest coding form for %s", genomic->name);
  typedef struct {
    int node, edge, saved;
    int64_t length, overlap;
  } Frame;
  int capacity = 0, depth = 1;
  Frame *stack = grow(NULL, &capacity, sizeof(*stack));
  stack[0] = (Frame){.node = terminal, .edge = g->nodes[terminal].first};
  while (depth) {
    Frame *f = &stack[depth - 1];
    if (f->node == root) {
      if (f->overlap >= CDS_RECIPROCAL_OVERLAP * length &&
          f->overlap >= CDS_RECIPROCAL_OVERLAP * f->length) {
        reverse_intervals(&form);
        if (form.count != reference.count ||
            memcmp(form.items, reference.items,
                   form.count * sizeof(*form.items)))
          if (write_isoform(output, genomic, &form, strand, gene, ++number) !=
              1)
            bye("cannot write alternative coding form for %s", genomic->name);
        reverse_intervals(&form);
      }
      f->edge = 0;
    }
    if (!f->edge) {
      form.count = f->saved;
      --depth;
      continue;
    }
    Link edge = g->edges[f->edge];
    f->edge = edge.next;
    Point a = c->points[g->nodes[edge.node].slot / STATES];
    Point b = c->points[g->nodes[f->node].slot / STATES];
    int saved = form.count;
    int64_t size = f->length, overlap = f->overlap;
    if ((a.kind & 3) != IS) {
      push_interval(&form, a.position, b.position);
      size += b.position - a.position;
      overlap += covered(&reference, prefix, b.position) -
                 covered(&reference, prefix, a.position);
    }
    int64_t possible = overlap + covered(&reference, prefix, a.position);
    if (possible < CDS_RECIPROCAL_OVERLAP * length ||
        possible < CDS_RECIPROCAL_OVERLAP * size ||
        size + g->nodes[edge.node].length < CDS_RECIPROCAL_OVERLAP * length) {
      form.count = saved;
      continue;
    }
    if (depth == capacity)
      stack = grow(stack, &capacity, sizeof(*stack));
    stack[depth++] =
        (Frame){edge.node, g->nodes[edge.node].first, saved, size, overlap};
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
  for (int donor = 0; donor < c->count; ++donor)
    if ((c->points[donor].kind & 3) == IS)
      for (int acceptor = donor + 1;
           acceptor < c->count &&
           c->points[acceptor].position - c->points[donor].position <=
               model->max_intron;
           ++acceptor)
        if ((c->points[acceptor].kind & 3) == IE &&
            c->points[acceptor].position - c->points[donor].position >= 20 &&
            junction_logit(model, s, c->points[donor].position,
                           c->points[acceptor].position) >= 0) {
          if (count + 1 >= capacity)
            junctions = grow(junctions, &capacity, sizeof(*junctions));
          junctions[++count] = (Link){acceptor, heads[donor]};
          heads[donor] = count;
        }
  Graph g = {.lookup = scalloc(smul((size_t)c->count, STATES), sizeof(int))};
  for (int start = 0; start < c->count; ++start) {
    if ((c->points[start].kind & 3) != START)
      continue;
    g.last = start;
    int root = graph_node(&g, start, 21);
    g.nodes[root].length = 0;
    for (int p = start; p <= g.last; ++p) {
      const Point *point = &c->points[p];
      int type = point->kind & 3;
      if (type == STOP) {
        int terminal = g.lookup[(size_t)p * STATES];
        if (terminal)
          emit_pair(&g, c, root, terminal, genomic, strand, output, genes);
        continue;
      }
      for (int state = 0; state < STATES; ++state) {
        int from = g.lookup[(size_t)p * STATES + state];
        if (!from)
          continue;
        if (type == IS) {
          for (int link = heads[p]; link; link = junctions[link].next)
            graph_link(&g, c, from, junctions[link].node, state);
        } else {
          const OpeningInfo *info = &c->openings[point->kind >> 2];
          int64_t stop = exon_stop(s, point, info, state);
          for (int j = p + 1; j < c->count && c->points[j].position <= stop &&
                              c->points[j].position <= info->bad;
               ++j) {
            int closing = c->points[j].kind & 3;
            if (closing != IS && closing != STOP)
              continue;
            int next = extend_exon(s, point, &c->points[j], state, stop);
            if (next >= 0)
              graph_link(&g, c, from, j, next);
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
  static const char *const names[SITES] = {"Start", "IS", "IE", "Stop"};
  FILE *files[2][SITES] = {{0}};
  for (int strand = 0; prefix && strand < 2; ++strand)
    for (int type = 0; type < SITES; ++type) {
      char path[PATH_BUFFER];
      if (snprintf(path, sizeof(path), "%s.%s.%s.bedgraph", prefix, names[type],
                   strand ? "minus" : "plus") >= (int)sizeof(path))
        bye("output prefix too long");
      if (!(files[strand][type] = fopen(path, "w")))
        bye("cannot create %s", path);
      fprintf(files[strand][type], "track type=bedGraph name=\"%s %c\"\n",
              names[type], strand ? '-' : '+');
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
      write_points(files[strand], s, &c, strand);
      assemble(model, &oriented, s, strand ? '-' : '+', &c, gff, &genes);
      free(c.openings);
      free(c.points);
    }
    free(reverse);
    cry("%s: %lu Start/Stop pairs", s->name, genes - previous);
  }
  for (int strand = 0; strand < 2; ++strand)
    for (int type = 0; type < SITES; ++type)
      if (files[strand][type] && fclose(files[strand][type]))
        bye("cannot close point bedGraph");
  return genes;
}
