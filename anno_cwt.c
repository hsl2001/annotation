#define _POSIX_C_SOURCE 200809L
#include "kseq.h"
#include <complex.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <zlib.h>

#define RFFT_IMPLEMENTATION
#include "rfft.h"

KSEQ_INIT(gzFile, gzread)

#ifndef CWT_SCALES
#define CWT_SCALES 4, 5, 6, 7, 8, 9, 16, 32, 64, 128
#endif
#ifndef MAX_WAVE_SIZE
#define MAX_WAVE_SIZE 9
#endif
static const double wave_scales[] = {CWT_SCALES};
#define WAVE_COUNT (sizeof(wave_scales) / sizeof(wave_scales[0]))
#define CWT_CHANNELS (2 * WAVE_COUNT)

typedef struct {
  double complex *kernel[WAVE_COUNT];
  int widths[WAVE_COUNT], max_width;
} Wavelets;
typedef struct { char *name, *seq; int length; } Contig;

#define fail(...) do { \
  fputs("Error: ", stderr); \
  fprintf(stderr, __VA_ARGS__); \
  fputc('\n', stderr); \
  exit(EXIT_FAILURE); \
} while (0)

static void *alloc(size_t count, size_t size) {
  if (size && count > SIZE_MAX / size) fail("Allocation size overflow");
  void *memory = calloc(count ? count : 1, size);
  if (!memory) fail("Out of memory");
  return memory;
}

static int base_index(char base) {
  switch (base) {
  case 'A': case 'a': return 0;
  case 'C': case 'c': return 1;
  case 'G': case 'g': return 2;
  case 'T': case 't': return 3;
  default: return -1;
  }
}

static double complex base_signal(char base) {
  const double complex mapping[4] = {1.0, I, -I, -1.0}; // 이거 나중에 꼭 매핑 바꿔봐야 함.
  int index = base_index(base);
  return index < 0 ? 0.0 : mapping[index];
}

static double complex morlet_integral(double lower, double upper) {
  static const double nodes[8] = {
    0.09501250983763744, 0.2816035507792589, 0.4580167776572274, 0.6178762444026438,
    0.7554044083550030, 0.8656312023878318, 0.9445750230732326, 0.9894009349916499
  };
  static const double weights[8] = {
    0.1894506104550685, 0.1826034150449236, 0.1691565193950025, 0.1495959888165767,
    0.1246289712555339, 0.0951585116824928, 0.0622535239386479, 0.0271524594117541
  };
  lower = fmax(lower, -8.0);
  upper = fmin(upper, 8.0);
  if (lower >= upper) return 0.0;
  double normalization = 1.0 / sqrt(sqrt(acos(-1.0)) *
                                   (1.0 + exp(-36.0) - 2.0 * exp(-27.0)));
  int segments = (int)ceil((upper - lower) * 2.0);
  double step = (upper - lower) / segments;
  double complex integral = 0.0;
  for (int segment = 0; segment < segments; segment++) {
    double middle = lower + (segment + 0.5) * step;
    for (int node = 0; node < 8; node++) {
      double left = middle - nodes[node] * step / 2.0;
      double right = middle + nodes[node] * step / 2.0;
      integral += weights[node] * step / 2.0 * (
        exp(-0.5 * left * left) * (cexp(-I * 6.0 * left) - exp(-18.0)) +
        exp(-0.5 * right * right) * (cexp(-I * 6.0 * right) - exp(-18.0)));
    }
  }
  return normalization * integral;
}

static void wavelets_init(Wavelets *wavelets) {
  memset(wavelets, 0, sizeof(*wavelets));
  for (size_t scale = 0; scale < WAVE_COUNT; scale++) {
    double dilation = wave_scales[scale];
    if (!isfinite(dilation) || dilation <= 0.0 || dilation != floor(dilation))
      fail("CWT scales must be positive natural numbers");
    double radius_value = ceil(8.0 * dilation + 0.5);
    if (radius_value > (INT_MAX - 1) / 2) fail("CWT scale exceeds supported kernel size");
    int radius = (int)radius_value;
    int width = 2 * radius + 1;
    size_t capacity = MAX_WAVE_SIZE > 0 ? (size_t)MAX_WAVE_SIZE : 1;
    if (capacity > (size_t)width) capacity = (size_t)width;
    while (capacity < (size_t)width) capacity *= 2;
    wavelets->kernel[scale] = alloc(capacity, sizeof(double complex));
    wavelets->widths[scale] = width;
    if (width > wavelets->max_width) wavelets->max_width = width;
    for (int tap = 0; tap < width; tap++) {
      double displacement = tap - radius;
      wavelets->kernel[scale][tap] = sqrt(dilation) * morlet_integral(
        (displacement - 0.5) / dilation, (displacement + 0.5) / dilation);
    }
  }
}

static size_t convolution_fft_size(size_t signal_length, size_t kernel_length) {
  if (signal_length > SIZE_MAX - (kernel_length - 1)) fail("CWT FFT size overflow");
  size_t required = signal_length + kernel_length - 1;
  size_t size = 1;
  while (size < required) {
    if (size > SIZE_MAX / 2) fail("CWT FFT size overflow");
    size <<= 1;
  }
  return size;
}

static void cwt_extract(const Wavelets *wavelets, const char *sequence, int length,
                        int start, int count, double *features) {
  if (count <= 0) return;
  memset(features, 0, (size_t)count * CWT_CHANNELS * sizeof(double));
  if (start >= length || (long)start + count <= 0) return;

  long context_start = (long)start - wavelets->max_width / 2;
  long context_end = (long)start + count + wavelets->max_width / 2;
  if (context_start < 0) context_start = 0;
  if (context_end > length) context_end = length;
  size_t signal_length = (size_t)(context_end - context_start);
  size_t fft_size = convolution_fft_size(signal_length, (size_t)wavelets->max_width);
  double complex *signal_fft = alloc(fft_size, sizeof(*signal_fft));
  double complex *kernel_fft = alloc(fft_size, sizeof(*kernel_fft));

  for (size_t index = 0; index < signal_length; index++)
    signal_fft[index] = base_signal(sequence[context_start + (long)index]);
  fft_transform(signal_fft, fft_size, false);

  for (size_t scale = 0; scale < WAVE_COUNT; scale++) {
    int width = wavelets->widths[scale];
    memset(kernel_fft, 0, fft_size * sizeof(*kernel_fft));
    for (int tap = 0; tap < width; tap++)
      kernel_fft[width - 1 - tap] = wavelets->kernel[scale][tap];
    fft_transform(kernel_fft, fft_size, false);
    for (size_t index = 0; index < fft_size; index++) kernel_fft[index] *= signal_fft[index];
    fft_transform(kernel_fft, fft_size, true);

    for (int offset = 0; offset < count; offset++) {
      long position = (long)start + offset;
      if (position < 0 || position >= length) continue;
      size_t convolution_index = (size_t)(position - context_start) + width - 1 - width / 2;
      double complex coefficient = kernel_fft[convolution_index] / (double)fft_size;
      features[(size_t)offset * CWT_CHANNELS + 2 * scale] = creal(coefficient);
      features[(size_t)offset * CWT_CHANNELS + 2 * scale + 1] = cimag(coefficient);
    }
  }

  free(signal_fft);
  free(kernel_fft);
}

static Contig *read_fasta(const char *path, int *count) {
  gzFile input = gzopen(!strcmp(path, "-") ? "/dev/stdin" : path, "rb");
  if (!input) fail("Cannot open FASTA: %s", path);
  kseq_t *record = kseq_init(input);
  Contig *contigs = NULL;
  *count = 0;
  while (kseq_read(record) >= 0) {
    if (!record->seq.l) continue;
    if (record->seq.l > INT_MAX - 1024) fail("Contig too long: %s", record->name.s);
    Contig *grown = realloc(contigs, (size_t)(*count + 1) * sizeof(*contigs));
    if (!grown) fail("Out of memory");
    contigs = grown;
    Contig *contig = contigs + (*count)++;
    contig->name = strdup(record->name.s);
    contig->seq = strdup(record->seq.s);
    contig->length = (int)record->seq.l;
    if (!contig->name || !contig->seq) fail("Out of memory");
  }
  kseq_destroy(record);
  if (gzclose(input) != Z_OK || !*count) fail("No readable FASTA records: %s", path);
  return contigs;
}

static FILE *create_in(const char *directory, const char *name) {
  char *path = alloc(strlen(directory) + strlen(name) + 2, 1);
  sprintf(path, "%s/%s", directory, name);
  FILE *file = fopen(path, "wb");
  if (!file) fail("Cannot create %s", path);
  free(path);
  return file;
}

static void write_floats(FILE *stream, const float *values, size_t count) {
  _Static_assert(sizeof(float) == 4 && FLT_RADIX == 2 && FLT_MANT_DIG == 24,
                 "CWT export requires binary32 floats");
  uint16_t endian = 1;
  if (*(unsigned char *)&endian) {
    if (fwrite(values, sizeof(float), count, stream) != count) fail("Cannot write CWT matrix");
    return;
  }
  for (size_t index = 0; index < count; index++) {
    unsigned char bytes[4], swapped[4];
    memcpy(bytes, values + index, 4);
    for (int byte = 0; byte < 4; byte++) swapped[byte] = bytes[3 - byte];
    if (fwrite(swapped, 1, 4, stream) != 4) fail("Cannot write CWT matrix");
  }
}

static void export_cwt(const Wavelets *wavelets, const char *fasta, const char *directory) {
  int count;
  Contig *contigs = read_fasta(fasta, &count);
  if (mkdir(directory, 0777)) fail("Cannot create directory %s", directory);
  FILE *matrix = create_in(directory, "matrix.bin");
  FILE *index = create_in(directory, "contigs.tsv");
  fputs("# anno-cwt-v2\n# scales\t", index);
  for (size_t scale = 0; scale < WAVE_COUNT; scale++) fprintf(index, "%s%.17g", scale ? "," : "", wave_scales[scale]);
  fputs("\nseqid\tlength\tplus_offset\tminus_offset\n", index);
  double *features = alloc(1024 * CWT_CHANNELS, sizeof(*features));
  float *values = alloc(1024 * CWT_CHANNELS, sizeof(*values));
  size_t rows = 0;
  for (int contig_index = 0; contig_index < count; contig_index++) {
    Contig *contig = contigs + contig_index;
    fprintf(index, "%s\t%d\t%zu\t%zu\n", contig->name, contig->length, rows, rows + contig->length);
    char *reverse = alloc((size_t)contig->length + 1, 1);
    for (int position = 0; position < contig->length; position++) {
      int base = base_index(contig->seq[contig->length - 1 - position]);
      reverse[position] = base < 0 ? 'N' : "TGCA"[base];
    }
    for (int strand = 0; strand < 2; strand++) {
      const char *bases = strand ? reverse : contig->seq;
      for (int start = 0; start < contig->length; start += 1024) {
        int length = contig->length - start < 1024 ? contig->length - start : 1024;
        cwt_extract(wavelets, bases, contig->length, start, length, features);
        for (size_t value = 0; value < (size_t)length * CWT_CHANNELS; value++) values[value] = (float)features[value];
        write_floats(matrix, values, (size_t)length * CWT_CHANNELS);
      }
    }
    free(reverse);
    rows += 2 * (size_t)contig->length;
  }
  if (fclose(matrix)) fail("Cannot finish CWT matrix");
  if (fclose(index)) fail("Cannot finish CWT index");
  for (int contig_index = 0; contig_index < count; contig_index++) {
    free(contigs[contig_index].name);
    free(contigs[contig_index].seq);
  }
  free(contigs);
  free(features);
  free(values);
  fprintf(stderr, "Saved %zu rows x %zu complex64 scales\n", rows, WAVE_COUNT);
}

int main(int argc, char **argv) {
  if (argc != 3) {
    fprintf(stderr, "Usage: %s <genome.fa[.gz]> <new_directory>\n", argv[0]);
    return EXIT_FAILURE;
  }
  Wavelets wavelets;
  wavelets_init(&wavelets);
  export_cwt(&wavelets, argv[1], argv[2]);
  for (size_t scale = 0; scale < WAVE_COUNT; scale++) free(wavelets.kernel[scale]);
  return EXIT_SUCCESS;
}