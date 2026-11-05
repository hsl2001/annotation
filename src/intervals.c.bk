/**
 * @file segments.c
 * @brief Shared operations for half-open genomic intervals.
 */
#include "include/anno_internal.h"

#include <stdlib.h>

/**
 * @brief Append a half-open segment to a growable list.
 * @param[in,out] list Segment list to extend.
 * @param[in] start Inclusive segment start.
 * @param[in] end Exclusive segment end.
 * @return None.
 */
void anno_segments_push(SegmentList *list, int64_t start, int64_t end) {
  if (list->count == list->capacity) {
    list->capacity = list->capacity ? list->capacity * 2 : 8;
    list->items =
        anno_xrealloc(list->items, (size_t)list->capacity, sizeof(*list->items));
  }
  list->items[list->count++] = (Segment){start, end};
}

/**
 * @brief Compare segments by start and then end.
 * @param[in] a First segment.
 * @param[in] b Second segment.
 * @return Negative, zero, or positive comparison result.
 */
int anno_segment_compare(const void *a, const void *b) {
  const Segment *x = a, *y = b;
  return x->start < y->start   ? -1
         : x->start > y->start ? 1
         : x->end < y->end     ? -1
                               : x->end > y->end;
}
