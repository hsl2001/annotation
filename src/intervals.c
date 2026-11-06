/**
 * @file intervals.c
 * @brief Half-openned genomic intervals
 */
#include "include/surf.h"
#include "include/utils.h"

/**
 * @brief Append an interval to a growable list.
 * @param[in,out] list: Interval list to extend.
 * @param[in] start: Inclusive interval start.
 * @param[in] end: Exclusive interval end.
 * @return None.
 */
void push_interval(IntervalList *list, int64_t start, int64_t end) {
  if (list->count == list->capacity) {
    list->capacity = list->capacity ? list->capacity * 2 : 8;
    list->items =
        srealloc(list->items, (size_t)list->capacity, sizeof(*list->items));
  }
  list->items[list->count++] = (Interval){start, end};
}

/**
 * @brief Compare intervals by start coordinate and then end coordinate (mainly
 * for sorting).
 * @param[in] a: First interval.
 * @param[in] b: Second interval.
 * @return Negative, zero, or positive comparison result.
 */
int compare_interval(const void *a, const void *b) {
  const Interval *x = a, *y = b;
  return x->start < y->start   ? -1
         : x->start > y->start ? 1
         : x->end < y->end     ? -1
                               : x->end > y->end;
}

char *splice_sequence(const Sequence *sequence, const IntervalList *exons,
                      char strand, int64_t *length) {
  *length = 0;
  for (int i = 0; i < exons->count; i++)
    *length += exons->items[i].end - exons->items[i].start;
  char *transcript = scalloc((size_t)*length + 1, 1);
  int64_t at = 0;
  for (int i = 0; i < exons->count; i++) {
    const Interval *e = &exons->items[strand == '+' ? i : exons->count - 1 - i];
    for (int64_t p = 0; p < e->end - e->start; p++)
      transcript[at++] = strand == '+'
                             ? sequence->bases[e->start + p]
                             : complement(sequence->bases[e->end - 1 - p]);
  }
  return transcript;
}
