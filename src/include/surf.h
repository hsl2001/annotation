/**
 * @file surf.h
 * @brief IO of surf:
 * REMEMBER: Coordinates are 0-based, half-opened, in contrast to GFF3
 * inputs and outputs have 1-based inclusive coordinates.
 *
 */
#ifndef SURF_H
#define SURF_H

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

/** Fized constants */
#define NUCLEOTIDES 4        // A, C, G, T
#define STRANDS 2            // Plus & minus strands
#define COMPLEX_COMPONENTS 2 // Real & imaginary components
#define PATH_BUFFER 4096     // Maximum handled path length
#define CLASSES 7            // Number of train labels (classes)
/** Constraints for gene grammar */
#define MIN_EXON_LENGTH 6    // Minimum exon length
#define MIN_INTRON_LENGTH 30 // Minimum intron length

/** Related to trainings */
#define DEFAULT_ITERATIONS 1000 // Maximum L-BFGS iterations
typedef enum {
  INTERGENIC = 0,
  PLUS_START_SINGLE_EXON = 1,
  PLUS_END_EXON = 2,
  PLUS_INTRON = 3,
  MINUS_START_SINGLE_EXON = 4,
  MINUS_END_EXON = 5,
  MINUS_INTRON = 6,

  /* Decode-only states: internal exons use the intergenic emission. */
  PLUS_INTERNAL_EXON = 7,
  MINUS_INTERNAL_EXON = 8
} Class;

/** Sequence records */
typedef struct {
  char *name;      // FASTA ID
  char *bases;     // !!Upper-case bases!!, NUL for termination
  int64_t length;  // number of bases
  uint8_t *labels; // Classes per-base (during training) or NULL
} Sequence;

/**
 * @brief Read a FASTA(.gz) file.
 * @param[in] path: FASTA path.
 * @param[out] count: Receives the number of records.
 * @param[out] output: Receives the allocated sequence array.
 */
void read_fasta(const char *path, int *count, Sequence **output);

/**
 * @brief Free records returned by read_fasta().
 * @param[in,out] sequences Records to free; NULL is accepted.
 * @param[in] count Number of records.
 * @return None.
 */
void free_sequences(Sequence *sequences, int count);

/**
 * @brief Fill sequence labels from a GFF3(.gz) annotation.
 * @param[in] path: GFF3 path.
 * @param[in,out] sequences: FASTA records whose per-base labels are filled.
 * @param[in] count: Number of FASTA records.
 * @return None.
 */
void read_gff3(const char *path, Sequence *sequences, int count);

int longest_orf(const char *transcript, int64_t length, int64_t *start,
                int64_t *end);

/**
 * @brief Convert a per-base class path into GFF3 transcript records.
 * @param[in,out] output: Destination GFF3 stream.
 * @param[in] sequence: Sequence and optional name metadata.
 * @param[in] path: Decoded classes, including inferred INTERNAL_EXON states.
 * @param[in,out] gene_counter: Decoded gene number; updated as records are
 * written.
 * @return number of transcripts written, or -1 on error.
 */
long write_gff(FILE *output, const Sequence *sequence, const uint8_t *path,
               unsigned long *gene_counter);

#endif /* SURF_H */
