/**
 * @file decoder.c
 * @brief Duration grammar construction and streaming Viterbi decoding.
 */
#include "anno_internal.h"

#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------------ */
/* CRF training and decoding                                               */
/* ------------------------------------------------------------------------ */

/**
 * @brief Add a scored predecessor edge to a grammar state.
 * @param[in,out] g Grammar to modify.
 * @param[in] from Predecessor state index.
 * @param[in] to Destination state index.
 * @return None.
 */
static void grammar_edge(Grammar *g, int from, int to) {
  int k = g->predecessors[to]++;
  if (k >= GRAMMAR_MAX_PREDECESSORS)
    anno_fail("grammar has too many predecessors");
  g->from[to][k] = from;
}

/**
 * @brief Add a self-looping chain of grammar states.
 * @param[in,out] g Grammar to modify.
 * @param[in] cls Class assigned to each new state.
 * @param[in] length Number of states in the chain.
 * @return Index of the first new state.
 */
static int grammar_chain(Grammar *g, uint8_t cls, int length) {
  int head = g->count;
  for (int i = 0; i < length; i++) {
    g->cls[g->count] = cls;
    if (i)
      grammar_edge(g, g->count - 1, g->count);
    g->count++;
  }
  grammar_edge(g, g->count - 1, g->count - 1);
  return head;
}

/**
 * @brief Build a grammar with start/single, internal, and end exons.
 * Plus: inter -> start_single [intron (internal intron)* end] -> inter.
 * Minus follows the reverse genomic order; single exons remain start_single.
 * @param[out] g Grammar to initialize.
 * @return None.
 */
void grammar_build(Grammar *g) {
  memset(g, 0, sizeof(*g));
  int inter = grammar_chain(g, ANNO_INTERGENIC, 1);
  for (int strand = 0; strand < ANNO_STRANDS; strand++) {
    uint8_t start_class = strand ? ANNO_MINUS_START_SINGLE_EXON
                                 : ANNO_PLUS_START_SINGLE_EXON;
    uint8_t end_class = strand ? ANNO_MINUS_END_EXON : ANNO_PLUS_END_EXON;
    uint8_t internal_class = strand ? ANNO_MINUS_INTERNAL_EXON
                                    : ANNO_PLUS_INTERNAL_EXON;
    uint8_t intron = strand ? ANNO_MINUS_INTRON : ANNO_PLUS_INTRON;
    int start = grammar_chain(g, start_class, ANNO_MIN_EXON_LENGTH);
    int end = grammar_chain(g, end_class, ANNO_MIN_EXON_LENGTH);
    int internal = grammar_chain(g, internal_class, ANNO_MIN_EXON_LENGTH);
    int intron_chain = grammar_chain(g, intron, ANNO_MIN_INTRON_LENGTH);
    int start_tail = start + ANNO_MIN_EXON_LENGTH - 1;
    int end_tail = end + ANNO_MIN_EXON_LENGTH - 1;
    int internal_tail = internal + ANNO_MIN_EXON_LENGTH - 1;
    int intron_tail = intron_chain + ANNO_MIN_INTRON_LENGTH - 1;
    grammar_edge(g, inter, start);
    grammar_edge(g, start_tail, inter);
    if (strand)
      grammar_edge(g, inter, end);
    grammar_edge(g, strand ? end_tail : start_tail, intron_chain);
    grammar_edge(g, intron_tail, internal);
    grammar_edge(g, internal_tail, intron_chain);
    grammar_edge(g, intron_tail, strand ? start : end);
    if (!strand)
      grammar_edge(g, end_tail, inter);
  }
  int shift = 0;
  for (int s = 0; s < g->count; s++) {
    int bits = 0;
    while ((1 << bits) < g->predecessors[s])
      bits++;
    g->shift[s] = shift;
    g->bits[s] = bits;
    shift += bits;
  }
  if (shift > 32)
    anno_fail("grammar too large");
}

struct Decoder {
  Grammar grammar;
  int64_t length, position;
  double *previous, *current;
  uint32_t *back;
};

/**
 * @brief Create a decoder for a sequence of the given length.
 * @param[in] length Number of sequence positions.
 * @return Newly allocated decoder.
 */
Decoder *decoder_create(int64_t length) {
  Decoder *d = anno_xcalloc(1, sizeof(*d));
  grammar_build(&d->grammar);
  d->length = length;
  d->previous = anno_xcalloc((size_t)d->grammar.count, sizeof(double));
  d->current = anno_xcalloc((size_t)d->grammar.count, sizeof(double));
  d->back = anno_xcalloc(length > 1 ? (size_t)(length - 1) : 1,
                     sizeof(*d->back));
  return d;
}

/**
 * @brief Free a decoder and its buffers.
 * @param[in] d Decoder to free, or NULL.
 * @return None.
 */
void decoder_destroy(Decoder *d) {
  if (!d)
    return;
  free(d->previous);
  free(d->current);
  free(d->back);
  free(d);
}

/**
 * @brief Map a decoded-only internal exon state to its observed class.
 * @param[in] cls Decoder grammar class.
 * @return Observed class used for CRF emissions and transitions.
 */
static uint8_t model_class(uint8_t cls) {
  return cls == ANNO_PLUS_INTERNAL_EXON || cls == ANNO_MINUS_INTERNAL_EXON
             ? ANNO_INTERGENIC : cls;
}

/**
 * @brief Consume one emission row and advance the Viterbi recurrence.
 * @param[in,out] d Streaming decoder state and path backpointers.
 * @param[in] emission Emission scores indexed by observed class.
 * @param[in] transition Learned class-transition scores, or NULL.
 * @return None.
 */
void decoder_step(Decoder *d, const double *emission,
                         const double *transition) {
  const Grammar *g = &d->grammar;
  int64_t i = d->position++;
  if (!i) {
    for (int s = 0; s < g->count; s++)
      d->previous[s] = -INFINITY;
    d->previous[0] = emission[ANNO_INTERGENIC];
    return;
  }
  uint32_t back = 0;
  double top = -INFINITY;
  for (int s = 0; s < g->count; s++) {
    double best = -INFINITY;
    int choice = 0;
    uint8_t to_class = model_class(g->cls[s]);
    for (int k = 0; k < g->predecessors[s]; k++) {
      int previous = g->from[s][k];
      uint8_t from_class = model_class(g->cls[previous]);
      double score = d->previous[previous];
      if (transition && from_class != to_class)
        score += transition[(size_t)from_class * ANNO_CLASSES + to_class];
      if (score > best) {
        best = score;
        choice = k;
      }
    }
    d->current[s] = best + emission[to_class];
    if (d->current[s] > top)
      top = d->current[s];
    back |= (uint32_t)choice << g->shift[s];
  }
  for (int s = 0; s < g->count; s++)
    d->current[s] -= top;
  d->back[i - 1] = back;
  double *swap = d->previous;
  d->previous = d->current;
  d->current = swap;
}



/**
 * @brief Trace back the best path ending in intergenic.
 * @param[in,out] d Completed decoder.
 * @param[out] path Receives one class per sequence position.
 * @return 0 on success, otherwise ERANGE.
 */
int decoder_finish(Decoder *d, uint8_t *path) {
  const Grammar *g = &d->grammar;
  if (d->position != d->length || !isfinite(d->previous[0]))
    return ERANGE;
  int state = 0;
  for (int64_t i = d->length - 1; i > 0; i--) {
    path[i] = g->cls[state];
    int choice =
        (int)((d->back[i - 1] >> g->shift[state]) &
              ((1u << g->bits[state]) - 1u));
    state = g->from[state][choice];
  }
  path[0] = g->cls[state];
  return state == 0 ? 0 : ERANGE;
}


