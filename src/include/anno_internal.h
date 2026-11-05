/**
 * @file anno_internal.h
 * @brief Private contracts shared by the anno implementation modules.
 */
#ifndef ANNO_INTERNAL_H
#define ANNO_INTERNAL_H

#include "anno.h"
#include "crf.h"

#include <stddef.h>

#define CRF_TRAIN_WINDOW ANNO_FEATURE_WINDOW
#define GRAMMAR_MAX_STATES \
  (1 + ANNO_STRANDS * (3 * ANNO_MIN_EXON_LENGTH + ANNO_MIN_INTRON_LENGTH))
#define GRAMMAR_MAX_PREDECESSORS (2 * ANNO_STRANDS + 1)

/** Half-open genomic interval shared by label and GFF modules. */
typedef struct {
  int64_t start;
  int64_t end;
} Segment;

/** Growable list of half-open genomic intervals. */
typedef struct {
  Segment *items;
  int count;
  int capacity;
} SegmentList;

/** Duration-expanded decoder grammar; internal exons map to background emissions. */
typedef struct {
  int count;
  uint8_t cls[GRAMMAR_MAX_STATES];
  int predecessors[GRAMMAR_MAX_STATES];
  int from[GRAMMAR_MAX_STATES][GRAMMAR_MAX_PREDECESSORS];
  int shift[GRAMMAR_MAX_STATES];
  int bits[GRAMMAR_MAX_STATES];
} Grammar;

typedef struct Decoder Decoder;

/**
 * @brief Report a fatal application error and terminate.
 * @param[in] format printf-style diagnostic format.
 * @param[in] ... Values referenced by @p format.
 * @return This function does not return.
 */
_Noreturn void anno_fail(const char *format, ...);
/**
 * @brief Write an informational log message to standard error.
 * @param[in] format printf-style message format.
 * @param[in] ... Values referenced by @p format.
 * @return None.
 */
void anno_log(const char *format, ...);
/**
 * @brief Multiply two sizes, terminating on overflow.
 * @param[in] a First factor.
 * @param[in] b Second factor.
 * @return The product of the two factors.
 */
size_t anno_mul_size(size_t a, size_t b);
/**
 * @brief Allocate zero-initialized memory or terminate on failure.
 * @param[in] count Number of elements.
 * @param[in] size Size of each element.
 * @return Allocated memory.
 */
void *anno_xcalloc(size_t count, size_t size);
/**
 * @brief Resize an allocation or terminate on failure.
 * @param[in] old Existing allocation, or NULL.
 * @param[in] count Number of elements.
 * @param[in] size Size of each element.
 * @return Resized allocation.
 */
void *anno_xrealloc(void *old, size_t count, size_t size);
/**
 * @brief Duplicate a string or terminate on allocation failure.
 * @param[in] text String to duplicate.
 * @return Caller-owned duplicate.
 */
char *anno_xstrdup(const char *text);
/**
 * @brief Read the monotonic clock.
 * @return Clock value in seconds.
 */
double anno_now_seconds(void);
/**
 * @brief Convert A/C/G/T to a base-channel index.
 * @param[in] base Nucleotide character.
 * @return Index 0..3, or -1 for another character.
 */
int anno_base_index(char base);
/**
 * @brief Return the DNA complement.
 * @param[in] base Nucleotide character.
 * @return Complementary nucleotide, or N for an unknown base.
 */
char anno_complement(char base);
/**
 * @brief Compute CWT power in position-major layout.
 * @param[in] bases Sequence bases.
 * @param[in] length Total sequence length.
 * @param[in] start Inclusive interval start.
 * @param[in] end Exclusive interval end.
 * @param[out] out Receives the interval's CWT power channels.
 * @return None.
 */
void anno_cwt_power(const char *bases, int64_t length, int64_t start,
                    int64_t end, float *out);
/**
 * @brief Fill channel-major one-hot and normalized CWT features.
 * @param[in] sequence Sequence bases and length.
 * @param[in] power Position-major CWT power values.
 * @param[in] power_start Genomic coordinate of the first power row.
 * @param[in] start First feature position.
 * @param[in] valid Number of valid positions.
 * @param[out] input Receives the channel-major feature matrix.
 * @return None.
 */
void anno_features(const AnnoSequence *sequence, const float *power,
                   int64_t power_start, int64_t start, int64_t valid,
                   float *input);

/**
 * @brief Append an interval to a growable list.
 * @param[in,out] list Interval list to extend.
 * @param[in] start Inclusive interval start.
 * @param[in] end Exclusive interval end.
 * @return None.
 */
void anno_segments_push(SegmentList *list, int64_t start, int64_t end);
/**
 * @brief Compare intervals by start coordinate and then end coordinate.
 * @param[in] a First interval.
 * @param[in] b Second interval.
 * @return Negative, zero, or positive comparison result.
 */
int anno_segment_compare(const void *a, const void *b);

/**
 * @brief Build the duration-expanded grammar.
 * @param[out] grammar Receives the initialized grammar.
 * @return None.
 */
void grammar_build(Grammar *grammar);
/**
 * @brief Convert a decoder grammar to the CRF observed-class representation.
 * @param[in] grammar Duration-expanded decoder grammar.
 * @param[out] classes Receives the observed class for each state.
 * @param[out] counts Receives each state's predecessor count.
 * @param[out] predecessors Receives flattened predecessor indices.
 * @return A CRF grammar view referencing the output arrays.
 */
CrfGrammar crf_grammar_view(const Grammar *grammar, uint8_t *classes,
                            uint8_t *counts, int *predecessors);
/**
 * @brief Create a streaming decoder for a sequence.
 * @param[in] length Number of sequence positions.
 * @return Newly allocated decoder.
 */
Decoder *decoder_create(int64_t length);
/**
 * @brief Release a decoder and its state buffers.
 * @param[in,out] decoder Decoder to release; NULL is accepted.
 * @return None.
 */
void decoder_destroy(Decoder *decoder);
/**
 * @brief Advance the decoder by one emission row.
 * @param[in,out] decoder Decoder state to update.
 * @param[in] emission Scores indexed by observed class.
 * @param[in] transition Class-pair transition weights, or NULL.
 * @return None.
 */
void decoder_step(Decoder *decoder, const double *emission,
                  const double *transition);
/**
 * @brief Trace back the best complete path.
 * @param[in,out] decoder Completed decoder state.
 * @param[out] path Receives one decoded class per sequence position.
 * @return Zero on success, or ERANGE if no complete path exists.
 */
int decoder_finish(Decoder *decoder, uint8_t *path);

/**
 * @brief Check whether a codon is a stop codon.
 * @param[in] codon Three-base codon.
 * @return Nonzero for TAA, TAG, or TGA.
 */
int is_stop(const char *codon);
/**
 * @brief Find the longest complete ORF in spliced transcript coordinates.
 * @param[in] transcript Spliced transcript bases.
 * @param[in] length Transcript length.
 * @param[out] start Receives the ORF start.
 * @param[out] end Receives the exclusive ORF end.
 * @return Nonzero if an ORF was found.
 */
int longest_orf(const char *transcript, int64_t length, int64_t *start,
                int64_t *end);

/**
 * @brief Train a CRF from labeled sequences.
 * @param[in] sequences Sequences with per-base labels.
 * @param[in] count Number of sequences.
 * @param[in] iterations Maximum optimizer iterations.
 * @param[in] threads Number of parallel objective workers.
 * @return Newly allocated trained model.
 */
CrfModel *train_crf(const AnnoSequence *sequences, int count, int iterations,
                    int threads);
/**
 * @brief Serialize a finite CRF model.
 * @param[in] path Destination model path.
 * @param[in] model Model to serialize.
 * @return Zero on success, otherwise an errno-style error code.
 */
int crf_save_model(const char *path, const CrfModel *model);
/**
 * @brief Load and validate a CRF model.
 * @param[in] path Source model path.
 * @param[out] model Receives the loaded model.
 * @return Zero on success, otherwise an errno-style error code.
 */
int crf_load_model(const char *path, CrfModel *model);
/**
 * @brief Predict sequences, writing GFF3 and optional class bedGraphs.
 * @param[in] model Trained CRF model.
 * @param[in] sequences Genome sequences to predict.
 * @param[in] count Number of sequences.
 * @param[in,out] gff Destination GFF3 stream.
 * @param[in] prefix Optional bedGraph output prefix, or NULL.
 * @return Number of genes written.
 */
unsigned long predict_genome(const CrfModel *model,
                             const AnnoSequence *sequences, int count,
                             FILE *gff, const char *prefix);

#endif /* ANNO_INTERNAL_H */
