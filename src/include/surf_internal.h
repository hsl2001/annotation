/**
 * @file surf_internal.h
 * @brief Private contracts shared within codebase
 */
#ifndef SURF_INTERNAL_H
#define SURF_INTERNAL_H

#include "surf.h"
#include "crf.h"

#include <stddef.h>

#define GRAMMAR_MAX_STATES \
(1 + STRANDS * (3 * MIN_EXON_LENGTH + MIN_INTRON_LENGTH))
#define GRAMMAR_MAX_PREDECESSORS (2 * STRANDS + 1)

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
void push_interval(IntervalList *list, int64_t start, int64_t end);

/**
 * @brief Compare intervals by start coordinate and then end coordinate.
 * @param[in] a: First interval.
 * @param[in] b: Second interval.
 * @return Negative, zero, or positive comparison result.
 */
int compare_interval(const void *a, const void *b);

/** Grammar for genes */
typedef struct {
    int count;
    uint8_t cls[GRAMMAR_MAX_STATES];
    int predecessors[GRAMMAR_MAX_STATES];
    int from[GRAMMAR_MAX_STATES][GRAMMAR_MAX_PREDECESSORS];
    int shift[GRAMMAR_MAX_STATES];
    int bits[GRAMMAR_MAX_STATES];
} Grammar;
// TODO:이 grammar안의 predecessor가 진짜 필요???

/**
* @brief Build the duration-expanded grammar.
* @param[out] grammar: Receives the initialized grammar.
* @return None.
*/
void build_grammar(Grammar *grammar);

// typedef struct Decoder Decoder;

/**
 * @brief Report an error and terminate.
 * @param[in] format: printf-style diagnostic format.
 * @param[in] ... Values referenced by @p format.
 * @return !!This function does not return!!
 */
_Noreturn void bye(const char *format, ...);

/**
 * @brief Write an informational log message to standard error.
 * @param[in] format: printf-style message format.
 * @param[in] ... Values referenced by @p format.
 * @return None.
 */
void cry(const char *format, ...);

/**
 * @brief Safely multiply two sizes, terminating on overflow.
 * @param[in] a: First factor.
 * @param[in] b: Second factor.
 * @return The product of the two factors.
 */
size_t smul(size_t a, size_t b);

/**
 * @brief Safely allocate zero-initialized memory or terminate on failure.
 * @param[in] count: Number of elements.
 * @param[in] size: Size of each element.
 * @return Allocated memory.
 */
void *scalloc(size_t count, size_t size);

/**
 * @brief Safely resize an allocation or terminate on failure.
 * @param[in] old: Existing allocation, or NULL.
 * @param[in] count: Number of elements.
 * @param[in] size: Size of each element.
 * @return Resized allocation.
 */
void *srealloc(void *old, size_t count, size_t size);

/**
 * @brief Safely duplicate a string or terminate on allocation failure.
 * @param[in] text: String to duplicate.
 * @return Caller-owned duplicate.
 */
char *sstrdup(const char *text);

/**
 * @brief Read the monotonic clock.
 * @return Clock value in seconds.
 */
double whattime(void);

/**
 * @brief Convert A/C/G/T to a base-channel index.
 * @param[in] base: Nucleotide character.
 * @return Index 0..3, or -1 for another character.
 */
int base_index(char base);

/**
 * @brief Return the DNA complement.
 * @param[in] base: Nucleotide character.
 * @return Complementary nucleotide, or N for an unknown base.
 */
char complement(char base);

/**
 * @brief Fill channel-major one-hot and normalized CWT features.
 * @param[in] sequence: Sequence bases and length.
 * @param[in] power: Position-major CWT power values.
 * @param[in] power_start: Genomic coordinate of the first power row.
 * @param[in] start: First feature position.
 * @param[in] valid: Number of valid positions.
 * @param[out] features: Receives the channel-major feature matrix.
 * @return None.
 */
void compute_features(const Sequence *sequence, const float *power,
                   int64_t power_start, int64_t start, int64_t valid,
                   float *features);



#endif /* surf_INTERNAL_H */
