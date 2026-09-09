#ifndef ANNO_DATA_H
#define ANNO_DATA_H

#include <stdio.h>

typedef struct {
  char *name, *bases;
  int length;
  unsigned char *coding[2], *states[2];
} Sequence;

typedef struct {
  Sequence *sequences;
  int count, selected_transcripts;
} Dataset;

Dataset *dataset_read(const char *path);
void dataset_free(Dataset *dataset);
void annotation_read(Dataset *dataset, FILE *stream);

#endif