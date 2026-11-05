/**
 * @file intervals.h
 * @brief Half-openned genomic intervals
 */
#ifndef INTERVALS_H
#define INTERVALS_H

#include "utils.c"

/** Half-open genomic intervals (for both labels and sequences) */
typedef struct {
    int64_t start;
    int64_t end;
} Interval;

/** List of half-open genomic intervals (for both labels and sequences) */
typedef struct {
    Interval *items;
    int count;
    int capacity;
} IntervalList;

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
 * @brief Compare intervals by start coordinate and then end coordinate.
 * @param[in] a: First interval.
 * @param[in] b: Second interval.
 * @return Negative, zero, or positive comparison result.
 */
int compare_inverval(const void *a, const void *b) {
  const Interval *x = a, *y = b;
  return x->start < y->start   ? -1
         : x->start > y->start ? 1
         : x->end < y->end     ? -1
                               : x->end > y->end;
}

#endif /* INTERVALS_H */
