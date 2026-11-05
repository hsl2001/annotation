/**
 * @file sequence.c
 * @brief FASTA parsing and sequence-record ownership.
 */
#define _POSIX_C_SOURCE 200809L

#include "include/anno_internal.h"

#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

/* ------------------------------------------------------------------------ */
/* FASTA                                                                    */
/* ------------------------------------------------------------------------ */

/**
 * @brief Read one line from a gzFile, growing the buffer as needed.
 * @param[in] input Open compressed or plain-text input stream.
 * @return Allocated line including its newline, or NULL at end of file.
 */
static char *gz_line(gzFile input) {
  size_t length = 0, capacity = 256;
  char *line = anno_xcalloc(capacity, 1);
  int ch;
  while ((ch = gzgetc(input)) != -1) {
    if (length + 1 >= capacity) {
      capacity = anno_mul_size(capacity, 2);
      line = anno_xrealloc(line, capacity, 1);
    }
    line[length++] = (char)ch;
    if (ch == '\n')
      break;
  }
  if (!length) {
    free(line);
    return NULL;
  }
  line[length] = '\0';
  return line;
}

/**
 * @brief Read FASTA records from a possibly compressed file.
 * @param[in] path FASTA file path.
 * @param[out] count_out Receives the number of records.
 * @return Caller-owned sequence array.
 */
AnnoSequence *anno_read_fasta(const char *path, int *count_out) {
  gzFile input = gzopen(path, "rb");
  if (!input)
    anno_fail("cannot open FASTA %s", path);
  AnnoSequence *sequences = NULL;
  int count = 0;
  char *line = NULL, *name = NULL, *bases = NULL;
  size_t length = 0, capacity = 0;
  while ((line = gz_line(input)) != NULL) {
    size_t line_length = strcspn(line, "\r\n");
    if (line[0] == '>') {
      if (name) {
        if (!length)
          anno_fail("empty FASTA sequence %s", name);
        sequences = anno_xrealloc(sequences, (size_t)count + 1, sizeof(*sequences));
        sequences[count++] = (AnnoSequence){name, bases, (int64_t)length, NULL};
        name = NULL;
        bases = NULL;
        length = capacity = 0;
      }
      char *first = line + 1;
      while (*first && isspace((unsigned char)*first))
        first++;
      char *end = first;
      while (*end && !isspace((unsigned char)*end))
        end++;
      if (end == first)
        anno_fail("empty FASTA sequence ID in %s", path);
      *end = '\0';
      name = anno_xstrdup(first);
      for (int i = 0; i < count; i++)
        if (!strcmp(sequences[i].name, name))
          anno_fail("duplicate FASTA sequence ID %s", name);
    } else if (name) {
      for (size_t i = 0; i < line_length; i++)
        if (!isspace((unsigned char)line[i])) {
          if (length + 1 >= capacity) {
            capacity = capacity ? capacity * 2 : 256;
            bases = anno_xrealloc(bases, capacity, 1);
          }
          bases[length++] = (char)toupper((unsigned char)line[i]);
        }
    } else if (line_length) {
      anno_fail("FASTA sequence before first header in %s", path);
    }
    free(line);
  }
  if (name) {
    if (!length)
      anno_fail("empty FASTA sequence %s", name);
    sequences = anno_xrealloc(sequences, (size_t)count + 1, sizeof(*sequences));
    sequences[count++] = (AnnoSequence){name, bases, (int64_t)length, NULL};
  }
  if (gzclose(input) != Z_OK)
    anno_fail("cannot read FASTA %s", path);
  if (!count)
    anno_fail("no FASTA records in %s", path);
  *count_out = count;
  return sequences;
}

/**
 * @brief Release an array of sequence records and their owned data.
 * @param[in,out] sequences Sequence records, or NULL.
 * @param[in] count Number of records.
 * @return None.
 */
void anno_sequences_free(AnnoSequence *sequences, int count) {
  for (int i = 0; sequences && i < count; i++) {
    free(sequences[i].name);
    free(sequences[i].bases);
    free(sequences[i].labels);
  }
  free(sequences);
}
