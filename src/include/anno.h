/**
 * @file anno.h
 * @brief Public API of anno: CWT features, L-BFGS CRF training, grammar
 * decoding, ORF resolution, and GFF3 output.
 *
 * Coordinates are 0-based half-open everywhere inside the program; GFF3
 * input and output use 1-based inclusive coordinates.
 */
#ifndef ANNO_H
#define ANNO_H

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

/** @name Fixed sizes */
/**@{*/
#define ANNO_NUCLEOTIDES 4        /**< A/C/G/T base channels */
#define ANNO_STRANDS 2            /**< plus and minus strands */
#define ANNO_COMPLEX_COMPONENTS 2 /**< real and imaginary components */
#define ANNO_PATH_BUFFER 4096    /**< maximum generated path length */
#define ANNO_CLASSES 7           /**< model output / training-label classes */
#define ANNO_CWT_CHANNELS 3      /**< Morlet CWT scales (power channels) */
#define ANNO_INPUT_CHANNELS 7    /**< A/C/G/T one-hot + CWT power */
#define ANNO_CWT_WINDOW 1024     /**< CWT block processed by one job */
#define ANNO_FEATURE_WINDOW 4096 /**< feature and CRF training window */
#define ANNO_CWT_RADIUS_FACTOR 4 /**< Morlet radius in scale units */
#define ANNO_MORLET_FREQUENCY 6.0 /**< Morlet central frequency */
/**@}*/

/** @name Grammar constants */
/**@{*/
#define ANNO_MIN_EXON_LENGTH 6    /**< minimum coding-exon run (bp)          */
#define ANNO_MIN_INTRON_LENGTH 30 /**< minimum intron length (bp) */

/**@}*/

/** @name Training constants */
/**@{*/
#define ANNO_DEFAULT_ITERATIONS 50 /**< maximum L-BFGS iterations */
/**@}*/

/** Model labels. Start/end refer to transcription order, not genomic order. */
typedef enum {
  ANNO_INTERGENIC = 0,
  ANNO_PLUS_START_SINGLE_EXON = 1,
  ANNO_PLUS_END_EXON = 2,
  ANNO_PLUS_INTRON = 3,
  ANNO_MINUS_START_SINGLE_EXON = 4,
  ANNO_MINUS_END_EXON = 5,
  ANNO_MINUS_INTRON = 6,
  /* Decoded-only states: internal exons use the intergenic emission. */
  ANNO_PLUS_INTERNAL_EXON = 7,
  ANNO_MINUS_INTERNAL_EXON = 8
} AnnoClass;

/** One sequence record. The struct owns every pointer it holds. */
typedef struct {
  char *name;      /**< FASTA ID                                     */
  char *bases;     /**< upper-case bases, NUL terminated             */
  int64_t length;  /**< number of bases                              */
  uint8_t *labels; /**< per-base AnnoClass (training only) or NULL   */
} AnnoSequence;

/**
 * @brief Morlet CWT of the complex base signal for positions [start,end).
 * @param[in] bases Input sequence bases.
 * @param[in] length Total sequence length.
 * @param[in] start First position, inclusive.
 * @param[in] end Last position, exclusive.
 * @param[out] out Receives (end-start) x ANNO_CWT_CHANNELS x 2 floats (real, imaginary).
 * @return None.
 *
 * The reverse-complement strand never needs recomputation: its transform
 * at position i equals the complex conjugate of this one at length-1-i.
 */
void anno_cwt(const char *bases, int64_t length, int64_t start, int64_t end,
              float *out);


/**
 * @brief Convert a per-base class path into GFF3 transcript records.
 * @param[in,out] output Destination GFF3 stream.
 * @param[in] sequence Sequence and optional name metadata.
 * @param[in] path Decoded classes, including inferred INTERNAL_EXON states.
 * @param[in,out] gene_counter Running gene number; updated as records are written.
 * @return number of transcripts written, or -1 on error.
 */
long anno_write_gff(FILE *output, const AnnoSequence *sequence,
                    const uint8_t *path, unsigned long *gene_counter);

/**
 * @brief Read a FASTA(.gz) file.
 * @param[in] path FASTA path.
 * @param[out] count Receives the number of records.
 * @return Array owned by the caller, or NULL only if the process does not fail.
 */
AnnoSequence *anno_read_fasta(const char *path, int *count);
/**
 * @brief Fill sequence labels from a GFF3(.gz) annotation.
 * @param[in] path GFF3 path.
 * @param[in,out] sequences FASTA records whose per-base labels are filled.
 * @param[in] count Number of FASTA records.
 * @return None.
 */
void anno_read_labels(const char *path, AnnoSequence *sequences, int count);
/**
 * @brief Free records returned by anno_read_fasta().
 * @param[in,out] sequences Records to free; NULL is accepted.
 * @param[in] count Number of records.
 * @return None.
 */
void anno_sequences_free(AnnoSequence *sequences, int count);

#endif
