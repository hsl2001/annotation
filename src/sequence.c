/**
 * @file sequence.c
 * @brief FASTA parsing and sequence-records
 */
#ifndef SEQUENCE_C
#define SEQUENCE_C

#include "include/surf.h"
#include "include/utils.h"

#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

/**
 * @brief Read one line from a gzFile, growing the buffer as needed.
 * @param[in] input: Open compressed or plain-text input stream.
 * @return Allocated line including its newline, or NULL at end of file.
 */
static char *gz_line(gzFile input) {
  size_t length = 0, capacity = 256;
  char *line = scalloc(capacity, 1);
  int ch;
  while ((ch = gzgetc(input)) != -1) {
    if (length + 1 >= capacity) {
      capacity = smul(capacity, 2);
      line = srealloc(line, capacity, 1);
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

void read_fasta(const char *path, int *count_out, Sequence **out) {
  gzFile input = gzopen(path, "rb");
  if (!input)
    bye("cannot open FASTA %s", path);
  Sequence *sequences = NULL;
  int count = 0;
  char *line = NULL, *name = NULL, *bases = NULL;
  size_t length = 0, capacity = 0;
  while ((line = gz_line(input)) != NULL) {
    size_t line_length = strcspn(line, "\r\n");
    if (line[0] == '>') {
      if (name) {
        if (!length)
          bye("empty FASTA sequence %s", name);
        sequences = srealloc(sequences, (size_t)count + 1, sizeof(*sequences));
        sequences[count++] =
            (Sequence){.name = name, .bases = bases, .length = (int64_t)length};
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
        bye("empty FASTA sequence ID in %s", path);
      *end = '\0';
      name = sstrdup(first);
      for (int i = 0; i < count; i++)
        if (!strcmp(sequences[i].name, name))
          bye("duplicate FASTA sequence ID %s", name);
    } else if (name) {
      for (size_t i = 0; i < line_length; i++)
        if (!isspace((unsigned char)line[i])) {
          if (length + 1 >= capacity) {
            capacity = capacity ? capacity * 2 : 256;
            bases = srealloc(bases, capacity, 1);
          }
          bases[length++] = (char)toupper((unsigned char)line[i]);
        }
    } else if (line_length) {
      bye("FASTA sequence before first header in %s", path);
    }
    free(line);
  }
  if (name) {
    if (!length)
      bye("empty FASTA sequence %s", name);
    sequences = srealloc(sequences, (size_t)count + 1, sizeof(*sequences));
    sequences[count++] =
        (Sequence){.name = name, .bases = bases, .length = (int64_t)length};
  }
  if (gzclose(input) != Z_OK)
    bye("cannot read FASTA %s", path);
  if (!count)
    bye("no FASTA records in %s", path);
  *count_out = count;
  *out = sequences;
}

void free_sequences(Sequence *sequences, int count) {
  for (int i = 0; sequences && i < count; i++) {
    free(sequences[i].name);
    free(sequences[i].bases);
    free(sequences[i].labels);
    for (int strand = 0; strand < 2; ++strand) {
      free(sequences[i].introns[strand].items);
      free(sequences[i].spans[strand].items);
    }
  }
  free(sequences);
}

#endif /* SEQUENCE_C */
