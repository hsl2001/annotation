#define _POSIX_C_SOURCE 200809L
#include "anno_signal.h"
#include "kseq.h"
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <zlib.h>

KSEQ_INIT(gzFile, gzread)

#define POWER_THRESHOLD 1.0
#define EXON_WINDOW 200
#define MIN_EXON_PEAKS 50

static FILE *output_file(const char *directory, const char *name) {
  char *path = anno_alloc(strlen(directory) + strlen(name) + 2, 1);
  sprintf(path, "%s/%s", directory, name);
  FILE *stream = fopen(path, "wb");
  if (!stream) anno_fail("Cannot create %s", path);
  free(path);
  return stream;
}

static void write_matrix(FILE *stream, float *values, size_t count) {
  _Static_assert(sizeof(float) == 4 && FLT_RADIX == 2 && FLT_MANT_DIG == 24 && FLT_MAX_EXP == 128,
                 "CWT export requires IEEE binary32 floats");
  uint16_t endian = 1;
  if (*(unsigned char *)&endian) {
    if (fwrite(values, sizeof(float), count, stream) != count) anno_fail("Cannot write CWT matrix");
    return;
  }
  for (size_t index = 0; index < count; index++) {
    unsigned char bytes[4];
    memcpy(bytes, values + index, 4);
    unsigned char swap = bytes[0]; bytes[0] = bytes[3]; bytes[3] = swap;
    swap = bytes[1]; bytes[1] = bytes[2]; bytes[2] = swap;
    if (fwrite(bytes, 1, 4, stream) != 4) anno_fail("Cannot write CWT matrix");
  }
}

static void annotate(const char *name, const float *powers, int length, unsigned long *genes) {
  float *scores = anno_alloc(length, sizeof(*scores));
  for (int position = 0; position < length; position++) {
    double score = 0.0;
    for (int scale = 0; scale < WAVE_COUNT; scale++)
      score += powers[(size_t)position * WAVE_COUNT + scale] / WAVE_COUNT;
    scores[position] = (float)score;
  }
  int *peaks = anno_alloc(length, sizeof(*peaks)), peak_count = 0;
  for (int position = 0; position < length; position++) {
    if (position < 1 || position + 1 >= length || !isfinite(scores[position]) ||
        scores[position] < POWER_THRESHOLD || scores[position] < scores[position - 1] ||
        scores[position] < scores[position + 1]) continue;
    if (scores[position] == scores[position - 1] && position > 1 && scores[position - 1] >= scores[position - 2]) continue;
    peaks[peak_count++] = position;
  }
  int exon_count = 0;
  for (int peak_index = 0; peak_index < peak_count;) {
    int cluster_start = peak_index;
    int cluster_end = cluster_start + 1;
    while (cluster_end < peak_count &&
           peaks[cluster_end] - peaks[cluster_start] <= EXON_WINDOW)
      cluster_end++;
    if (cluster_end - cluster_start < MIN_EXON_PEAKS) {
      peak_index++;
      continue;
    }
    int start = peaks[cluster_start], end = peaks[cluster_end - 1];
    unsigned long gene = ++*genes;
    double confidence = fmin(scores[start], scores[end]);
    int score = (int)lrint(fmin(1000.0, fmax(0.0, confidence * 100.0)));
    printf("%s\t%d\t%d\texon%lu\t%d\t.\n", name, start, end, gene, score);
    exon_count++;
    peak_index = cluster_end;
  }
    fprintf(stderr, "Fixed threshold: %s, threshold %.6g, %d peaks, %d exon clusters\n",
      name, POWER_THRESHOLD, peak_count, exon_count);
  free(peaks);
  free(scores);
}

int main(int argc, char **argv) {
  int help = argc == 2 && (!strcmp(argv[1], "--help") || !strcmp(argv[1], "-h"));
  int export = argc == 4 && !strcmp(argv[1], "cwt");
  if (help || (!export && argc != 2)) {
    fprintf(stderr, "Usage: %s <genome.fasta[.gz]|-> > candidates.bed\n"
                    "       %s cwt <genome.fasta[.gz]|-> <new_directory>\n", argv[0], argv[0]);
    return help ? 0 : 1;
  }
  const char *fasta = argv[export ? 2 : 1];
  gzFile input = gzopen(!strcmp(fasta, "-") ? "/dev/stdin" : fasta, "rb");
  if (!input) anno_fail("Cannot open FASTA: %s", fasta);
  FILE *matrix = NULL, *index = NULL;
  if (export) {
    if (mkdir(argv[3], 0777)) anno_fail("Cannot create new directory %s: %s", argv[3], strerror(errno));
    matrix = output_file(argv[3], "matrix.bin");
    index = output_file(argv[3], "contigs.tsv");
    fprintf(index, "# anno-cwt-v1\n# kernel_widths\t%d,%d,%d", wave_sizes[0], wave_sizes[1], wave_sizes[2]);
    fputs("\nseqid\tlength\tplus_offset\tminus_offset\n", index);
  }
  Wavelets wavelets;
  wavelets_init(&wavelets);
  kseq_t *record = kseq_init(input);
  char **names = NULL;
  size_t records = 0, rows = 0;
  unsigned long genes = 0;
  int result;
  while ((result = kseq_read(record)) >= 0) {
    if (!record->seq.l) continue;
    if (record->seq.l > INT_MAX - MAX_WAVE_SIZE) anno_fail("FASTA contig too long");
    int length = (int)record->seq.l;
    for (size_t previous = 0; previous < records; previous++)
      if (!strcmp(names[previous], record->name.s)) anno_fail("Duplicate FASTA identifier: %s", record->name.s);
    if (records == SIZE_MAX / sizeof(*names)) anno_fail("Too many FASTA records");
    char **grown = realloc(names, (records + 1) * sizeof(*names));
    if (!grown) anno_fail("Out of memory");
    names = grown;
    names[records] = strdup(record->name.s);
    if (!names[records++]) anno_fail("Out of memory");
    if ((size_t)length > (SIZE_MAX / (WAVE_COUNT * 8) - rows) / 2) anno_fail("CWT matrix size overflow");
    if (export) fprintf(index, "%s\t%d\t%zu\t%zu\n", record->name.s, length, rows, rows + length);
    float *powers = export ? NULL : anno_alloc((size_t)length * WAVE_COUNT, sizeof(float));
    int strand_count = export ? 2 : 1;
    for (int strand = 0; strand < strand_count; strand++) {
      char *reverse = strand ? reverse_complement(record->seq.s, length) : NULL;
      const char *bases = strand ? reverse : record->seq.s;
      double features[1024 * CWT_CHANNELS];
      float values[1024 * CWT_CHANNELS];
      for (int start = 0; start < length;) {
        int count = length - start < 1024 ? length - start : 1024;
        cwt_extract(&wavelets, bases, length, start, count, features);
        if (export) {
          size_t count_values = (size_t)count * CWT_CHANNELS;
          for (size_t value = 0; value < count_values; value++) values[value] = (float)features[value];
          write_matrix(matrix, values, count_values);
        } else {
          for (int offset = 0; offset < count; offset++) {
            int position = start + offset;
            int known = position >= MAX_WAVE_SIZE / 2 &&
                        position + (MAX_WAVE_SIZE - 1) / 2 < length;
            if (known) for (int tap = -MAX_WAVE_SIZE / 2; tap <= (MAX_WAVE_SIZE - 1) / 2; tap++)
              if (base_index(bases[position + tap]) < 0) known = 0;
            for (int scale = 0; scale < WAVE_COUNT; scale++) {
              double real = features[offset * CWT_CHANNELS + 2 * scale];
              double imaginary = features[offset * CWT_CHANNELS + 2 * scale + 1];
              powers[(size_t)position * WAVE_COUNT + scale] =
                  known ? (float)(real * real + imaginary * imaginary) : NAN;
            }
          }
        }
        start += count;
      }
      if (!export) annotate(record->name.s, powers, length, &genes);
      rows += length;
      free(reverse);
    }
    fprintf(stderr, "%s: %s, %d bp\n", export ? "CWT" : "Annotated", record->name.s, length);
    free(powers);
  }
  int error;
  gzerror(input, &error);
  if (result < -1 || (error != Z_OK && error != Z_STREAM_END)) anno_fail("Malformed or unreadable FASTA");
  kseq_destroy(record);
  if (gzclose(input) != Z_OK || !records) anno_fail("No readable nonempty FASTA records");
  for (size_t previous = 0; previous < records; previous++) free(names[previous]);
  free(names);
  if (export) {
    int failed = ferror(matrix) || ferror(index);
    if (fclose(matrix)) failed = 1;
    if (fclose(index)) failed = 1;
    if (failed) anno_fail("Cannot finish CWT output");
    fprintf(stderr, "Saved %zu rows x %d complex64 scales\n", rows, WAVE_COUNT);
  } else {
    if (fflush(stdout) || ferror(stdout)) anno_fail("Cannot write GFF output");
    fprintf(stderr, "Predicted %lu exon candidates from threshold peak clusters\n", genes);
  }
  return 0;
}