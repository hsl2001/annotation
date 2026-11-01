/**
 * @file anno.h
 * @brief Public API of anno: CWT + 1D UNet emissions, 9-state grammar
 *        Viterbi decoding, longest-ORF resolution and GFF3 output.
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
#define ANNO_CLASSES 9            /**< UNet classes == decoded grammar states */
#define ANNO_CWT_CHANNELS 3       /**< Morlet CWT scales (power channels)    */
#define ANNO_INPUT_CHANNELS 7     /**< A/C/G/T one-hot + CWT power           */
#define ANNO_CWT_WINDOW 1024      /**< CWT block processed by one job        */
#define ANNO_UNET_WINDOW 65536    /**< UNet training / inference window      */
#define ANNO_UNET_MARGIN 4096     /**< inference context discarded per side  */
#define ANNO_MAX_DEPTH 6          /**< maximum number of UNet levels         */
/**@}*/

/** @name Grammar constants */
/**@{*/
#define ANNO_MIN_EXON_LENGTH 6    /**< minimum coding-exon run (bp)          */
#define ANNO_MIN_INTRON_LENGTH 30 /**< minimum intron length (bp)            */
#define ANNO_EMISSION_EPSILON 1e-6f
#define ANNO_SPLICE_BONUS 8.0     /**< log advantage of GT..AG over other splice motifs */
#define ANNO_TRANSCRIPT_PENALTY 8.0 /**< log penalty for opening a transcript */
/**@}*/

/** @name Training constants */
/**@{*/
#define ANNO_DEFAULT_ITERATIONS 1000
#define ANNO_BATCH_WINDOWS 1
#define ANNO_LEARNING_RATE 3e-3f
#define ANNO_FOCAL_GAMMA 2.0f
#define ANNO_DICE_EPSILON 1e-5f
#define ANNO_DICE_WEIGHT 0.25     /**< weight of the soft-F1 term vs focal */
#define ANNO_GRADIENT_CLIP 1.0f
#define ANNO_SEED 2026u
/**@}*/

/** UNet output classes; each class is also one decoded grammar state. */
typedef enum {
  ANNO_INTERGENIC = 0,
  ANNO_PLUS_5UTR = 1,
  ANNO_PLUS_EXON = 2,  /**< coding part of a plus-strand exon  */
  ANNO_PLUS_INTRON = 3,
  ANNO_PLUS_3UTR = 4,
  ANNO_MINUS_5UTR = 5,
  ANNO_MINUS_EXON = 6, /**< coding part of a minus-strand exon */
  ANNO_MINUS_INTRON = 7,
  ANNO_MINUS_3UTR = 8
} AnnoClass;

/** Shape of the 1D UNet. */
typedef struct {
  int depth;                       /**< number of resolution levels       */
  int kernel;                      /**< odd convolution width             */
  int pool;                        /**< average-pool / upsample factor    */
  int channels[ANNO_MAX_DEPTH];    /**< feature channels per level        */
} UnetConfig;

/** 1D UNet with its input normalisation; opaque. */
typedef struct Unet Unet;

/** One sequence record. The struct owns every pointer it holds. */
typedef struct {
  char *name;       /**< FASTA ID                                     */
  char *bases;      /**< upper-case bases, NUL terminated             */
  int64_t length;   /**< number of bases                              */
  uint8_t *labels;  /**< per-base AnnoClass (training only) or NULL   */
} AnnoSequence;

/** Default UNet shape used for training. */
UnetConfig unet_default_config(void);
/** @brief Create a randomly initialised UNet; NULL on invalid config. */
Unet *unet_create(const UnetConfig *config);
/** @brief Release a UNet created by unet_create() or unet_load(). */
void unet_destroy(Unet *unet);
/**
 * @brief Predict class probabilities.
 * @param input       length x ANNO_INPUT_CHANNELS normalised features.
 * @param probability length x ANNO_CLASSES output (rows sum to one).
 * @return 0 on success, an errno value otherwise.
 */
int unet_predict(const Unet *unet, const float *input, int64_t length,
                 float *probability);
/** @brief Read a binary .ann model. @return 0 on success. */
int unet_load(const char *path, Unet **out);
/** @brief Write a binary .ann model. @return 0 on success. */
int unet_save(const char *path, const Unet *unet);

/**
 * @brief Morlet CWT of the complex base signal for positions [start,end).
 * @param out (end-start) x ANNO_CWT_CHANNELS x 2 floats (real, imaginary).
 *
 * The reverse-complement strand never needs recomputation: its transform
 * at position i equals the complex conjugate of this one at length-1-i.
 */
void anno_cwt(const char *bases, int64_t length, int64_t start, int64_t end,
              float *out);

/**
 * @brief Grammar-constrained Viterbi decoding of a probability matrix.
 * @param bases       sequence (for splice motifs).
 * @param probability length x ANNO_CLASSES.
 * @param path        length output classes (AnnoClass).
 * @return 0 on success.
 */
int anno_decode(const char *bases, int64_t length, const float *probability,
                uint8_t *path);

/**
 * @brief Convert a per-base class path into GFF3 transcript records.
 * @param gene_counter running gene number, updated.
 * @return number of transcripts written, or -1 on error.
 */
long anno_write_gff(FILE *output, const AnnoSequence *sequence,
                    const uint8_t *path, unsigned long *gene_counter);

/**
 * @brief Read a FASTA(.gz) file.
 * @return array owned by the caller (anno_sequences_free()).
 */
AnnoSequence *anno_read_fasta(const char *path, int *count);
/** @brief Fill sequences[].labels from a GFF3(.gz) annotation. */
void anno_read_labels(const char *path, AnnoSequence *sequences, int count);
/** @brief Free an array returned by anno_read_fasta(). */
void anno_sequences_free(AnnoSequence *sequences, int count);

#endif
