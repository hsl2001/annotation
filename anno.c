#define _POSIX_C_SOURCE 200809L
#include "anno.h"
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

const int wave_sizes[WAVE_COUNT] = {4, 5,8, 9};

void anno_fail(const char *format, ...) {
  va_list arguments;
  va_start(arguments, format);
  fputs("Error: ", stderr);
  vfprintf(stderr, format, arguments);
  fputc('\n', stderr);
  va_end(arguments);
  exit(EXIT_FAILURE);
}

void *anno_alloc(size_t count, size_t size) {
  if (size && count > SIZE_MAX / size)
    anno_fail("Allocation size overflow");
  void *memory = calloc(count ? count : 1, size);
  if (!memory)
    anno_fail("Out of memory");
  return memory;
}

int base_index(char base) {
  switch (base) {
  case 'A': case 'a': return 0;
  case 'C': case 'c': return 1;
  case 'G': case 'g': return 2;
  case 'T': case 't': return 3;
  default: return -1;
  }
}

double complex base_signal(char base) {
  const double complex mapping[4] = {1.0, I, -I, -1.0};
  int index = base_index(base);
  return index < 0 ? 0.0 : mapping[index];
}

char *reverse_complement(const char *sequence, int length) {
  char *result = anno_alloc((size_t)length + 1, sizeof(char));
  for (int position = 0; position < length; position++) {
    int base = base_index(sequence[length - position - 1]);
    result[position] = base < 0 ? 'N' : "TGCA"[base];
  }
  return result;
}

void wavelets_init(Wavelets *wavelets) {
  memset(wavelets, 0, sizeof(*wavelets));
  for (int scale = 0; scale < WAVE_COUNT; scale++) {
    int width = wave_sizes[scale];
    double complex mean = 0.0;
    for (int tap = 0; tap < width; tap++) {
      double coordinate = (tap - width / 2.0 + 0.5) / (width / 8.0);
      double complex value = exp(-0.5 * coordinate * coordinate) *
                             (cexp(I * 6.0 * coordinate) - exp(-18.0));
      wavelets->kernel[scale][tap] = value;
      mean += value / width;
    }
    double energy = 0.0;
    for (int tap = 0; tap < width; tap++) {
      wavelets->kernel[scale][tap] -= mean;
      double magnitude = cabs(wavelets->kernel[scale][tap]);
      energy += magnitude * magnitude;
    }
    for (int tap = 0; tap < width; tap++)
      wavelets->kernel[scale][tap] =
          conj(wavelets->kernel[scale][tap]) / sqrt(energy);
  }
}

void cwt_extract(const Wavelets *wavelets, const char *sequence, int length,
                 int start, int count, double *features) {
  memset(features, 0, (size_t)count * CWT_CHANNELS * sizeof(double));
  for (int offset = 0; offset < count; offset++) {
    long position = (long)start + offset;
    if (position < 0 || position >= length)
      continue;
    for (int scale = 0; scale < WAVE_COUNT; scale++) {
      int width = wave_sizes[scale];
      double complex coefficient = 0.0;
      for (int tap = 0; tap < width; tap++) {
        long context = position + tap - width / 2;
        if (context >= 0 && context < length)
          coefficient += base_signal(sequence[context]) *
                         wavelets->kernel[scale][tap];
      }
      features[(size_t)offset * CWT_CHANNELS + 2 * scale] = creal(coefficient);
      features[(size_t)offset * CWT_CHANNELS + 2 * scale + 1] = cimag(coefficient);
    }
  }
}

int anno_call_exons_with_params(const float *powers, int length, int wave_count,
                                const AnnoParameters *parameters, Exon *exons) {
  float *scores = anno_alloc(length, sizeof(*scores));
  for (int position = 0; position < length; position++) {
    double sum = 0.0;
    for (int scale = 0; scale < wave_count; scale++)
      sum += powers[(size_t)position * wave_count + scale];
    scores[position] = (float)(sum / wave_count);
  }
  int *hits = anno_alloc(length, sizeof(*hits));
  for (int scale = 0; scale < wave_count; scale++) {
    for (int position = 1; position + 1 < length; position++) {
      float value = powers[(size_t)position * wave_count + scale];
      float previous = powers[(size_t)(position - 1) * wave_count + scale];
      float next = powers[(size_t)(position + 1) * wave_count + scale];
        if (!isfinite(value) || value < parameters->power_threshold ||
          value < previous || value < next)
        continue;
      if (value == previous && position > 1 &&
          previous >= powers[(size_t)(position - 2) * wave_count + scale])
        continue;
      hits[position]++;
    }
  }
  int *peaks = anno_alloc(length, sizeof(*peaks)), peak_count = 0;
  for (int position = 0; position < length; position++)
    if (hits[position] == wave_count) peaks[peak_count++] = position;
  int exon_count = 0;
  for (int peak_index = 0; peak_index < peak_count;) {
    int cluster_start = peak_index;
    int cluster_end = cluster_start + 1;
    while (cluster_end < peak_count &&
              peaks[cluster_end] - peaks[cluster_start] <= parameters->exon_window)
      cluster_end++;
            if (cluster_end - cluster_start < parameters->min_exon_peaks) {
      peak_index++;
      continue;
    }
    int start = peaks[cluster_start], end = peaks[cluster_end - 1];
    exons[exon_count].start = start;
    exons[exon_count].end = end;
    exons[exon_count].score = (float)fmin(scores[start], scores[end]);
    exon_count++;
    peak_index = cluster_end;
  }
  free(peaks);
  free(hits);
  free(scores);
  return exon_count;
}

int anno_call_exons(const float *powers, int length, int wave_count, Exon *exons) {
  const AnnoParameters defaults = {
      POWER_THRESHOLD, EXON_WINDOW, MIN_EXON_PEAKS};
  return anno_call_exons_with_params(powers, length, wave_count, &defaults, exons);
}

#ifndef ANNO_NO_MAIN
#include "kseq.h"
#include "ketopt.h"
#include <zlib.h>

KSEQ_INIT(gzFile, gzread)

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

static void annotate(const char *name, const float *powers, int length,
                     const AnnoParameters *parameters, unsigned long *genes) {
  Exon *exons = anno_alloc(length, sizeof(*exons));
  int exon_count = anno_call_exons_with_params(powers, length, WAVE_COUNT, parameters, exons);
  for (int index = 0; index < exon_count; index++) {
    unsigned long gene = ++*genes;
    int score = (int)lrint(fmin(1000.0, fmax(0.0, exons[index].score * 100.0)));
    printf("%s\t%d\t%d\texon%lu\t%d\t.\n", name, exons[index].start, exons[index].end, gene, score);
  }
    fprintf(stderr, "Boundary peaks: %s, threshold %.6g, window %d, min peaks %d, %d exon clusters\n",
      name, parameters->power_threshold, parameters->exon_window,
      parameters->min_exon_peaks, exon_count);
  free(exons);
}

static int parse_int_option(const char *name, const char *text, int minimum) {
  char *end = NULL;
  long value = strtol(text, &end, 10);
  if (*text == '\0' || *end != '\0' || value < minimum || value > INT_MAX)
    anno_fail("Invalid %s: %s", name, text);
  return (int)value;
}

static float parse_float_option(const char *name, const char *text, float minimum) {
  char *end = NULL;
  float value = strtof(text, &end);
  if (*text == '\0' || *end != '\0' || !isfinite(value) || value < minimum)
    anno_fail("Invalid %s: %s", name, text);
  return value;
}

int main(int argc, char **argv) {
  AnnoParameters parameters = {POWER_THRESHOLD, EXON_WINDOW, MIN_EXON_PEAKS};
  ketopt_t options = KETOPT_INIT;
  int option;
  while ((option = ketopt(&options, argc, argv, 0, "t:w:m:h", NULL)) >= 0) {
    if (option == 'h') {
      fprintf(stderr, "Usage: %s [-t FLOAT] [-w INT] [-m INT] <genome.fasta[.gz]|->\n"
                      "       %s cwt <genome.fasta[.gz]|-> <new_directory>\n", argv[0], argv[0]);
      return 0;
    }
    if (option == '?') anno_fail("Unknown option");
    if (option == ':') anno_fail("Missing option value");
    if (option == 't') parameters.power_threshold = parse_float_option("-t", options.arg, 0.0f);
    else if (option == 'w') parameters.exon_window = parse_int_option("-w", options.arg, 1);
    else if (option == 'm') parameters.min_exon_peaks = parse_int_option("-m", options.arg, 1);
  }
  int argument = options.ind;
  int export = argument < argc && !strcmp(argv[argument], "cwt");
  if (argc - argument != (export ? 3 : 1)) {
    fprintf(stderr, "Usage: %s [--threshold FLOAT] [--window INT] [--min-peaks INT] <genome.fasta[.gz]|->\n"
                    "       %s cwt <genome.fasta[.gz]|-> <new_directory>\n", argv[0], argv[0]);
    return 1;
  }
  const char *fasta = argv[argument + (export ? 1 : 0)];
  gzFile input = gzopen(!strcmp(fasta, "-") ? "/dev/stdin" : fasta, "rb");
  if (!input) anno_fail("Cannot open FASTA: %s", fasta);
  FILE *matrix = NULL, *index = NULL;
  if (export) {
    if (mkdir(argv[3], 0777)) anno_fail("Cannot create new directory %s: %s", argv[3], strerror(errno));
    matrix = output_file(argv[3], "matrix.bin");
    index = output_file(argv[3], "contigs.tsv");
    fputs("# anno-cwt-v1\n# kernel_widths\t", index);
    for (int scale = 0; scale < WAVE_COUNT; scale++)
      fprintf(index, "%s%d", scale ? "," : "", wave_sizes[scale]);
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
      if (!export) annotate(record->name.s, powers, length, &parameters, &genes);
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
#endif