#define _POSIX_C_SOURCE 200809L
#include "anno.h"
#include "kseq.h"
#include "ketopt.h"
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <zlib.h>

#define RFFT_IMPLEMENTATION
#include "rfft.h"
#define LBFGS_IMPLEMENTATION
#include "lbfgs.h"
#define CRF_IMPLEMENTATION
#include "crf.h"

KSEQ_INIT(gzFile, gzread)

const int wave_sizes[WAVE_COUNT] = {4, 5, 6, 7, 8, 9, 16, 32, 64, 128};

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
  if (size && count > SIZE_MAX / size) anno_fail("Allocation size overflow");
  void *memory = calloc(count ? count : 1, size);
  if (!memory) anno_fail("Out of memory");
  return memory;
}

static void *anno_resize(void *memory, size_t count, size_t size) {
  if (size && count > SIZE_MAX / size) anno_fail("Allocation size overflow");
  void *resized = realloc(memory, count * size);
  if (!resized) anno_fail("Out of memory");
  return resized;
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

static double complex base_signal(char base) {
  const double complex mapping[4] = {1.0, I, -I, -1.0};
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

void wavelets_init(Wavelets *wavelets) {
  memset(wavelets, 0, sizeof(*wavelets));
  for (int scale = 0; scale < WAVE_COUNT; scale++) {
    double dilation = wave_sizes[scale];
    if (!isfinite(dilation) || dilation <= 0.0 || dilation != floor(dilation))
      anno_fail("CWT scales must be positive natural numbers");
    double radius_value = ceil(8.0 * dilation + 0.5);
    if (radius_value > (INT_MAX - 1) / 2.0) anno_fail("CWT scale exceeds supported kernel size");
    int radius = (int)radius_value;
    int width = 2 * radius + 1;
    size_t capacity = MAX_WAVE_SIZE > 0 ? (size_t)MAX_WAVE_SIZE : 1;
    if (capacity > (size_t)width) capacity = (size_t)width;
    while (capacity < (size_t)width) capacity *= 2;
    wavelets->kernel[scale] = anno_alloc(capacity, sizeof(double complex));
    wavelets->widths[scale] = width;
    if (width > wavelets->max_width) wavelets->max_width = width;
    for (int tap = 0; tap < width; tap++) {
      double displacement = tap - radius;
      wavelets->kernel[scale][tap] = sqrt(dilation) * morlet_integral(
        (displacement - 0.5) / dilation, (displacement + 0.5) / dilation);
    }
  }
}

static void wavelets_free(Wavelets *wavelets) {
  for (int scale = 0; scale < WAVE_COUNT; scale++) free(wavelets->kernel[scale]);
}

static size_t convolution_fft_size(size_t signal_length, size_t kernel_length) {
  if (signal_length > SIZE_MAX - (kernel_length - 1)) anno_fail("CWT FFT size overflow");
  size_t required = signal_length + kernel_length - 1;
  size_t size = 1;
  while (size < required) {
    if (size > SIZE_MAX / 2) anno_fail("CWT FFT size overflow");
    size <<= 1;
  }
  return size;
}

void cwt_extract(const Wavelets *wavelets, const char *sequence, int length, int start, int count, double *features) {
  if (count <= 0) return;
  memset(features, 0, (size_t)count * CWT_CHANNELS * sizeof(double));
  if (start >= length || (long)start + count <= 0) return;

  long context_start = (long)start - wavelets->max_width / 2;
  long context_end = (long)start + count + wavelets->max_width / 2;
  if (context_start < 0) context_start = 0;
  if (context_end > length) context_end = length;
  size_t signal_length = (size_t)(context_end - context_start);
  size_t fft_size = convolution_fft_size(signal_length, (size_t)wavelets->max_width);
  double complex *signal_fft = anno_alloc(fft_size, sizeof(*signal_fft));
  double complex *kernel_fft = anno_alloc(fft_size, sizeof(*kernel_fft));

  for (size_t index = 0; index < signal_length; index++)
    signal_fft[index] = base_signal(sequence[context_start + (long)index]);
  fft_transform(signal_fft, fft_size, false);

  for (int scale = 0; scale < WAVE_COUNT; scale++) {
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

static double cwt_power(const double *coefficient) {
  return coefficient[0] * coefficient[0] + coefficient[1] * coefficient[1];
}

static float *cwt_powers(const Wavelets *wavelets, const char *sequence, int length, int block) {
  float *powers = anno_alloc((size_t)length, WAVE_COUNT * sizeof(*powers));
  double *features = anno_alloc((size_t)block * CWT_CHANNELS, sizeof(*features));
  for (int start = 0; start < length; start += block) {
    int count = length - start < block ? length - start : block;
    cwt_extract(wavelets, sequence, length, start, count, features);
    for (int offset = 0; offset < count; offset++) {
      int position = start + offset;
      const double *coefficient = features + (size_t)offset * CWT_CHANNELS;
      for (int scale = 0; scale < WAVE_COUNT; scale++) {
        powers[(size_t)position * WAVE_COUNT + scale] = (float)cwt_power(coefficient + 2 * scale);
      }
    }
  }
  free(features);
  return powers;
}

Contig *read_fasta(const char *path, int *count) {
  gzFile input = gzopen(!strcmp(path, "-") ? "/dev/stdin" : path, "rb");
  if (!input) anno_fail("Cannot open FASTA: %s", path);
  kseq_t *record = kseq_init(input);
  Contig *contigs = NULL;
  *count = 0;
  while (kseq_read(record) >= 0) {
    if (!record->seq.l) continue;
    if (record->seq.l > (size_t)INT_MAX) anno_fail("Contig too long: %s", record->name.s);
    contigs = anno_resize(contigs, (size_t)*count + 1, sizeof(*contigs));
    Contig *c = contigs + (*count)++;
    memset(c, 0, sizeof(*c));
    c->name = strdup(record->name.s);
    c->seq = strdup(record->seq.s);
    c->length = (int)record->seq.l;
    if (!c->name || !c->seq) anno_fail("Out of memory");
  }
  kseq_destroy(record);
  if (gzclose(input) != Z_OK || !*count) anno_fail("No readable FASTA records: %s", path);
  return contigs;
}

static void free_contigs(Contig *contigs, int count) {
  for (int contig = 0; contig < count; contig++) {
    free(contigs[contig].name);
    free(contigs[contig].seq);
    free(contigs[contig].cwt_power);
    free(contigs[contig].label);
  }
  free(contigs);
}

typedef struct { int contig, start, end; char strand, *parent; } Cds;
typedef struct { int first, count, total; } Group;

static int by_parent(const void *x, const void *y) {
  const Cds *a = x, *b = y;
  int order = strcmp(a->parent, b->parent);
  return order ? order : a->start - b->start;
}

static int by_total(const void *x, const void *y) { return ((const Group *)y)->total - ((const Group *)x)->total; }

/* Label CDS and introns from GFF3, one transcript per locus (longest CDS wins). */
int label_cds(Contig *contigs, int count, const char *path, int *skipped) {
  gzFile input = gzopen(path, "rb");
  if (!input) anno_fail("Cannot open GFF: %s", path);
  kstream_t *stream = ks_init(input);
  kstring_t line = {0, 0, NULL};
  Cds *cds = NULL;
  int n = 0, separator;
  while (ks_getuntil(stream, KS_SEP_LINE, &line, &separator) >= 0) {
    char *field[9], *cursor = line.s;
    int fields = 0;
    if (!line.l || line.s[0] == '#') continue;
    while (fields < 9 && cursor) {
      field[fields++] = cursor;
      if ((cursor = strchr(cursor, '\t'))) *cursor++ = 0;
    }
    if (fields < 9 || strcmp(field[2], "CDS")) continue;
    char *parent = strstr(field[8], "Parent=");
    if (!parent) continue;
    parent += 7;
    parent[strcspn(parent, ";,")] = 0;
    int contig = 0;
    while (contig < count && strcmp(contigs[contig].name, field[0])) contig++;
    int start = atoi(field[3]), end = atoi(field[4]);
    if (contig == count || start < 1 || end < start || end > contigs[contig].length) continue;
    cds = anno_resize(cds, (size_t)n + 1, sizeof(*cds));
    if (!(cds[n].parent = strdup(parent))) anno_fail("Out of memory");
    cds[n].contig = contig, cds[n].start = start, cds[n].end = end, cds[n].strand = field[6][0];
    n++;
  }
  free(line.s);
  ks_destroy(stream);
  gzclose(input);
  qsort(cds, n, sizeof(*cds), by_parent);
  Group *groups = anno_alloc(n, sizeof(*groups));
  int g = 0;
  for (int i = 0; i < n; i++) {
    if (!i || strcmp(cds[i].parent, cds[i - 1].parent)) groups[g++].first = i;
    groups[g - 1].count++, groups[g - 1].total += cds[i].end - cds[i].start + 1;
  }
  qsort(groups, g, sizeof(*groups), by_total);
  int used = 0;
  for (int k = 0; k < g; k++) {
    Cds *s = cds + groups[k].first;
    Contig *c = contigs + s->contig;
    int m = groups[k].count, minus = s->strand == '-', ok = s->strand == '+' || minus;
    for (int j = 1; j < m; j++) ok &= s[j].contig == s->contig && s[j].strand == s->strand && s[j].start > s[j - 1].end;
    for (int p = s->start - 1; ok && p < s[m - 1].end; p++) ok &= !c->label[p];
    if (!ok) continue;
    used++;
    /* Reading frame is recomputed from CDS lengths in translation order; the GFF phase column is ignored. */
    int *startpos = anno_alloc(m, sizeof(int));
    for (int jj = 0, frame = 0; jj < m; jj++) {
      int j = minus ? m - 1 - jj : jj;
      startpos[j] = frame;
      frame = (frame + s[j].end - s[j].start + 1) % 3;
    }
    for (int j = 0; j < m; j++) {
      int pos = startpos[j], len = s[j].end - s[j].start + 1;
      for (int t = 0; t < len; t++, pos = (pos + 1) % 3)
        c->label[minus ? s[j].end - 1 - t : s[j].start - 1 + t] = (minus ? 7 : 1) + pos;
      int q = minus ? (startpos[j] + len - 1) % 3 : startpos[j];
      for (int p = j ? s[j - 1].end : s[j].start; p < s[j].start - 1; p++) c->label[p] = (minus ? 10 : 4) + q;
    }
    free(startpos);
  }
  for (int i = 0; i < n; i++) free(cds[i].parent);
  free(cds);
  free(groups);
  *skipped = g - used;
  return used;
}

static int is_cds(int state) { return (state >= 1 && state <= 3) || (state >= 7 && state <= 9); }

void write_gff(const Contig *c, const uint8_t *path, unsigned long *genes) {
  for (int i = 0; i < c->length;) {
    if (!path[i]) { i++; continue; }
    int minus = path[i] >= 7, start = i, coding = 0;
    while (i < c->length && path[i] && (path[i] >= 7) == minus) coding += is_cds(path[i++]);
    if (coding < MIN_CDS) continue;
    unsigned long g = ++*genes;
    char strand = minus ? '-' : '+';
    printf("%s\tanno\tgene\t%d\t%d\t.\t%c\t.\tID=gene%lu\n", c->name, start + 1, i, strand, g);
    printf("%s\tanno\tmRNA\t%d\t%d\t.\t%c\t.\tID=mRNA%lu;Parent=gene%lu\n", c->name, start + 1, i, strand, g, g);
    for (int j = start; j < i;) {
      if (!is_cds(path[j])) { j++; continue; }
      int k = j;
      while (k < i && is_cds(path[k])) k++;
      int phase = (3 - (minus ? path[k - 1] - 7 : path[j] - 1)) % 3;
      printf("%s\tanno\texon\t%d\t%d\t.\t%c\t.\tParent=mRNA%lu\n", c->name, j + 1, k, strand, g);
      printf("%s\tanno\tCDS\t%d\t%d\t.\t%c\t%d\tID=cds%lu;Parent=mRNA%lu\n", c->name, j + 1, k, strand, phase, g, g);
      j = k;
    }
  }
}

static void save_model(const char *path, const double *w) {
  float *values = anno_alloc(WEIGHTS, sizeof(*values));
  for (int i = 0; i < WEIGHTS; i++) values[i] = (float)w[i];
  FILE *file = fopen(path, "wb");
  if (!file || fwrite("ANNOCRFE", 1, 8, file) != 8 ||
      fwrite(values, sizeof(*values), WEIGHTS, file) != WEIGHTS || fclose(file))
    anno_fail("Cannot write model: %s", path);
  free(values);
}

static void load_model(const char *path, double *w) {
  float *values = anno_alloc(WEIGHTS, sizeof(*values));
  FILE *file = fopen(path, "rb");
  char magic[8];
  if (!file || fread(magic, 1, 8, file) != 8) anno_fail("Cannot read model: %s", path);
  if (memcmp(magic, "ANNOCRFE", 8)) anno_fail("Incompatible model; retrain with this executable: %s", path);
  if (fread(values, sizeof(*values), WEIGHTS, file) != WEIGHTS || fgetc(file) != EOF)
    anno_fail("Cannot read model: %s", path);
  fclose(file);
  for (int i = 0; i < WEIGHTS; i++) {
    w[i] = values[i];
    if (!isfinite(w[i])) anno_fail("Non-finite model weight: %s", path);
  }
  free(values);
}

static FILE *create_in(const char *directory, const char *name) {
  char *path = anno_alloc(strlen(directory) + strlen(name) + 2, 1);
  sprintf(path, "%s/%s", directory, name);
  FILE *file = fopen(path, "wb");
  if (!file) anno_fail("Cannot create %s", path);
  free(path);
  return file;
}

static void write_floats(FILE *stream, const float *values, size_t count) {
  _Static_assert(sizeof(float) == 4 && FLT_RADIX == 2 && FLT_MANT_DIG == 24, "CWT export requires binary32 floats");
  uint16_t endian = 1;
  if (*(unsigned char *)&endian) {
    if (fwrite(values, sizeof(float), count, stream) != count) anno_fail("Cannot write CWT matrix");
    return;
  }
  for (size_t index = 0; index < count; index++) {
    unsigned char bytes[4], swapped[4];
    memcpy(bytes, values + index, 4);
    for (int j = 0; j < 4; j++) swapped[j] = bytes[3 - j];
    if (fwrite(swapped, 1, 4, stream) != 4) anno_fail("Cannot write CWT matrix");
  }
}

/* Raw complex64 CWT export for plot_cwt.py: plus strand rows then reverse-complement rows. */
static void export_cwt(const Wavelets *wavelets, const char *fasta, const char *directory, int block) {
  int count;
  Contig *contigs = read_fasta(fasta, &count);
  if (mkdir(directory, 0777)) anno_fail("Cannot create directory %s: %s", directory, strerror(errno));
  FILE *matrix = create_in(directory, "matrix.bin"), *index = create_in(directory, "contigs.tsv");
  fputs("# anno-cwt-v2\n# scales\t", index);
  for (int scale = 0; scale < WAVE_COUNT; scale++) fprintf(index, "%s%d", scale ? "," : "", wave_sizes[scale]);
  fputs("\nseqid\tlength\tplus_offset\tminus_offset\n", index);
  double *features = anno_alloc((size_t)block * CWT_CHANNELS, sizeof(*features));
  float *values = anno_alloc((size_t)block * CWT_CHANNELS, sizeof(*values));
  size_t rows = 0;
  for (int c = 0; c < count; c++) {
    Contig *g = contigs + c;
    fprintf(index, "%s\t%d\t%zu\t%zu\n", g->name, g->length, rows, rows + g->length);
    char *reverse = anno_alloc((size_t)g->length + 1, 1);
    for (int i = 0; i < g->length; i++) {
      int base = base_index(g->seq[g->length - 1 - i]);
      reverse[i] = base < 0 ? 'N' : "TGCA"[base];
    }
    for (int strand = 0; strand < 2; strand++) {
      const char *bases = strand ? reverse : g->seq;
      for (int start = 0; start < g->length; start += block) {
        int n = g->length - start < block ? g->length - start : block;
        cwt_extract(wavelets, bases, g->length, start, n, features);
        for (size_t v = 0; v < (size_t)n * CWT_CHANNELS; v++) values[v] = (float)features[v];
        write_floats(matrix, values, (size_t)n * CWT_CHANNELS);
      }
    }
    free(reverse);
    rows += 2 * (size_t)g->length;
  }
  if (fclose(matrix) || fclose(index)) anno_fail("Cannot finish CWT export");
  free_contigs(contigs, count);
  free(features);
  free(values);
  fprintf(stderr, "Saved %zu rows x %d complex64 scales\n", rows, WAVE_COUNT);
}

static int usage(const char *program, int status) {
  fprintf(stderr, "Usage: %s [-m model] [-i initial_model] [-e iterations] [-p threads] [-b block] <genome.fa[.gz]> [<train.gff3[.gz]>] > genes.gff3\n"
                  "       %s cwt <genome.fa[.gz]> <new_directory>\n", program, program);
  return status;
}

int main(int argc, char **argv) {
  const char *model = "anno.model";
  const char *initial_model = NULL;
  int iterations = DEFAULT_ITERATIONS, threads = DEFAULT_THREADS, block = DEFAULT_BLOCK, option;
  ketopt_t options = KETOPT_INIT;
  while ((option = ketopt(&options, argc, argv, 0, "m:i:e:p:b:h", NULL)) >= 0) {
    if (option == 'm') model = options.arg;
    else if (option == 'i') initial_model = options.arg;
    else if (option == 'e' && atoi(options.arg) > 0) iterations = atoi(options.arg);
    else if (option == 'p' && atoi(options.arg) > 0) threads = atoi(options.arg);
    else if (option == 'b' && atoi(options.arg) > 0) block = atoi(options.arg);
    else return usage(argv[0], option != 'h');
  }
  int rest = argc - options.ind;
  char **args = argv + options.ind;
  int export = rest == 3 && !strcmp(args[0], "cwt");
  if (!export && (rest < 1 || rest > 2)) return usage(argv[0], 1);
  if (initial_model && rest != 2) anno_fail("Initial model is only supported for training");
  crf_init();
  Wavelets wavelets;
  wavelets_init(&wavelets);
  if (export) {
    export_cwt(&wavelets, args[1], args[2], block);
    wavelets_free(&wavelets);
    return 0;
  }
  double *weights = anno_alloc(WEIGHTS, sizeof(*weights));
  if (rest == 1) load_model(model, weights);
  if (initial_model) {
    load_model(initial_model, weights);
    fprintf(stderr, "Warm start from %s\n", initial_model);
  }
  int count;
  Contig *contigs = read_fasta(args[0], &count);
  for (int c = 0; c < count; c++)
    contigs[c].cwt_power = cwt_powers(&wavelets, contigs[c].seq, contigs[c].length, block);
  if (rest == 2) {
    for (int c = 0; c < count; c++) contigs[c].label = anno_alloc(contigs[c].length, 1);
    int skipped, used = label_cds(contigs, count, args[1], &skipped);
    fprintf(stderr, "Training on %d transcripts (%d overlapping or inconsistent skipped)\n", used, skipped);
    crf_train(contigs, count, weights, iterations, block, threads);
    save_model(model, weights);
  }
  puts("##gff-version 3");
  unsigned long genes = 0;
  for (int c = 0; c < count; c++) {
    uint8_t *path = anno_alloc(contigs[c].length, 1);
    crf_decode(weights, contigs + c, path, block);
    write_gff(contigs + c, path, &genes);
    free(path);
    fprintf(stderr, "Annotated %s: %d bp\n", contigs[c].name, contigs[c].length);
  }
  if (fflush(stdout) || ferror(stdout)) anno_fail("Cannot write GFF");
  fprintf(stderr, "Predicted %lu genes\n", genes);
  wavelets_free(&wavelets);
  free_contigs(contigs, count);
  free(weights);
  return 0;
}
