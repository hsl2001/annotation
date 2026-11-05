/**
 * @file surf_internal.h
 * @brief Private contracts shared within codebase
 */
#ifndef SURF_INTERNAL_H
#define SURF_INTERNAL_H

#include "surf.h"
#include "crf.h"

#include <stddef.h>

#define SHORTGRAMMER_MAX_STATES \
(1 + STRANDS * (3 * MIN_EXON_LENGTH + MIN_INTRON_LENGTH))
#define SHORTGRAMMAR_MAX_BACKPOINTERS (2 * STRANDS + 1)

/** Grammar for genes */
typedef struct {
    int count;
    uint8_t cls[SHORTGRAMMER_MAX_STATES];
    int backpointers[SHORTGRAMMER_MAX_STATES];
    int from[SHORTGRAMMER_MAX_STATES][SHORTGRAMMAR_MAX_BACKPOINTERS];
    int shift[SHORTGRAMMER_MAX_STATES];
    int bits[SHORTGRAMMER_MAX_STATES];
} ShortGrammar;

/**
* @brief Build the duration-expanded grammar.
* @param[out] grammar: Receives the initialized grammar.
* @return None.
*/
void build_grammar(ShortGrammar *grammar);

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

/**
 * @brief Predict sequences, writing GFF3 and optional class bedGraphs.
 * @param[in] model: Trained CRF model.
 * @param[in] sequences: Genome sequences to predict.
 * @param[in] count: Number of sequences.
 * @param[in,out] gff: Destination GFF3 stream.
 * @param[in] prefix: Optional bedGraph output prefix, or NULL.
 * @return Number of genes written.
 */
unsigned long predict_genome(const Model *model,
                             const Sequence *sequences, int count,
                             FILE *gff, const char *prefix);

#endif /* surf_INTERNAL_H */
