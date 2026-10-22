#ifdef __APPLE__
#define _DARWIN_C_SOURCE
#endif
#define _POSIX_C_SOURCE 200809L
#include "anno.h"
#include "gff_cwt_model.h"
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
#include <zlib.h>

#define RFFT_IMPLEMENTATION
#include "rfft.h"

const int wave_sizes[WAVE_COUNT] = {CWT_SCALES};
#define CWT_IMPLEMENTATION
#include "cwt.h"
#define CRF_IMPLEMENTATION
#include "crf.h"

_Static_assert(sizeof(GffCwtModel) == GFF_CWT_BASIS * CWT_CHANNELS * sizeof(float),
               "GffCwtModel must be a packed float array");

KSEQ_INIT(gzFile, gzread)

void anno_fail(const char *format, ...) {
  va_list arguments;
  va_start(arguments, format);
  fputs("Error: ", stderr);
  vfprintf(stderr, format, arguments);
  fputc('\n', stderr);
  va_end(arguments);
  exit(EXIT_FAILURE);
}

static void check_allocation_size(size_t count, size_t size) {
  if (size && count > SIZE_MAX / size) anno_fail("Allocation size overflow");
}

static char *anno_strdup(const char *text) {
  char *copy = strdup(text);
  if (!copy) anno_fail("Out of memory");
  return copy;
}

void *anno_alloc(size_t count, size_t size) {
  check_allocation_size(count ? count : 1, size);
  void *memory = calloc(count ? count : 1, size);
  if (!memory) anno_fail("Out of memory");
  return memory;
}

static void *anno_resize(void *memory, size_t count, size_t size) {
  check_allocation_size(count, size);
  void *resized = realloc(memory, count * size);
  if (!resized) anno_fail("Out of memory");
  return resized;
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
    *c = (Contig){
      .name = anno_strdup(record->name.s),
      .seq = anno_strdup(record->seq.s),
      .length = (int)record->seq.l
    };
  }
  kseq_destroy(record);
  if (gzclose(input) != Z_OK || !*count) anno_fail("No readable FASTA records: %s", path);
  return contigs;
}

static void free_contigs(Contig *contigs, int count) {
  for (int contig = 0; contig < count; contig++) {
    free(contigs[contig].name);
    free(contigs[contig].seq);
    free(contigs[contig].cwt);
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

static int by_total(const void *x, const void *y) {
  return ((const Group *)y)->total - ((const Group *)x)->total;
}

static int valid_cds(const char *seq, const Cds *s, int m, int total, int minus) {
  if (total < MIN_CDS || total % 3) return 0;
  for (int jj = 0, index = 0, codon = 0; jj < m; jj++) {
    int j = minus ? m - 1 - jj : jj, len = s[j].end - s[j].start + 1;
    for (int t = 0; t < len; t++, index++) {
      int base = base_index(seq[minus ? s[j].end - 1 - t : s[j].start - 1 + t]);
      if (base < 0) return 0;
      codon = (codon * 4 + (minus ? 3 - base : base)) & 63;
      if (index % 3 == 2) {
        int stop = codon == 48 || codon == 50 || codon == 56;
        if (index == 2 && codon != 14) return 0;
        if (index + 1 < total ? stop : !stop) return 0;
      }
    }
  }
  return 1;
}

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
    cds[n].parent = anno_strdup(parent);
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
    for (int j = 1; j < m; j++) ok &= s[j].contig == s->contig && s[j].strand == s->strand && s[j].start > s[j - 1].end + 2;
    if (ok) ok = valid_cds(c->seq, s, m, groups[k].total, minus);
    for (int p = s->start - 1; ok && p < s[m - 1].end; p++) ok &= !c->label[p];
    if (!ok) continue;
    used++;
    /* Reading frame is recomputed from CDS lengths in translation order; the GFF phase column is ignored. */
    for (int jj = 0, frame = 0; jj < m; jj++) {
      int j = minus ? m - 1 - jj : jj, len = s[j].end - s[j].start + 1, start_frame = frame;
      for (int t = 0; t < len; t++, frame = (frame + 1) % 3)
        c->label[minus ? s[j].end - 1 - t : s[j].start - 1 + t] = (minus ? 7 : 1) + frame;
      int q = minus ? (start_frame + len - 1) % 3 : start_frame;
      for (int p = j ? s[j - 1].end : s[j].start; p < s[j].start - 1; p++)
        c->label[p] = (minus ? 10 : 4) + q;
    }
  }
  for (int i = 0; i < n; i++) free(cds[i].parent);
  free(cds);
  free(groups);
  *skipped = g - used;
  return used;
}

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
  uint32_t count = WEIGHTS;
  FILE *file = fopen(path, "wb");
  if (!file || fwrite(&count, sizeof(count), 1, file) != 1 ||
      fwrite(w, sizeof(*w), WEIGHTS, file) != WEIGHTS || fclose(file))
    anno_fail("Cannot write model: %s", path);
}

static void load_model(const char *path, double *w) {
  uint32_t count = 0;
  FILE *file = fopen(path, "rb");
  if (!file || fread(&count, sizeof(count), 1, file) != 1 || count != WEIGHTS ||
      fread(w, sizeof(*w), WEIGHTS, file) != WEIGHTS || fgetc(file) != EOF || fclose(file))
    anno_fail("Cannot read model: %s (expected %d weights)", path, WEIGHTS);
  for (int i = 0; i < WEIGHTS; i++)
    if (!isfinite(w[i])) anno_fail("Non-finite model weight: %s", path);
}

static void write_gff_cwt_header(const char *path, const GffCwtModel *model) {
  const float *values = (const float *)model;
  FILE *file = fopen(path, "w");
  if (!file) anno_fail("Cannot write embedded GFF->CWT header: %s", path);
  fprintf(file, "#ifndef GFF_CWT_MODEL_H\n#define GFF_CWT_MODEL_H\n\n");
  fprintf(file, "static const float embedded_gff_cwt_values[%d] = {\n",
          GFF_CWT_BASIS * CWT_CHANNELS);
  for (int i = 0; i < GFF_CWT_BASIS * CWT_CHANNELS; i++) {
    if (!isfinite(values[i])) {
      fclose(file);
      anno_fail("Non-finite GFF->CWT coefficient at index %d", i);
    }
    fprintf(file, "  %#.9gf,%s", values[i], i % 6 == 5 ? "\n" : " ");
  }
  fprintf(file, "\n};\n\n#endif\n");
  if (fclose(file)) anno_fail("Cannot finish embedded GFF->CWT header: %s", path);
}

static int pretrain_cwt(int argc, char **argv) {
  if (argc < 5 || (argc - 3) % 2) {
    fprintf(stderr, "Usage: %s pretrain-cwt <output-header.h> <fasta> <gff3> [<fasta> <gff3> ...]\n",
            argv[0]);
    return 1;
  }
  Wavelets wavelets;
  wavelets_init(&wavelets);
  GffCwtAccumulator accumulator;
  gff_cwt_accumulator_init(&accumulator);
  for (int arg = 3; arg < argc; arg += 2) {
    int count;
    Contig *contigs = read_fasta(argv[arg], &count);
    for (int c = 0; c < count; c++)
      contigs[c].label = anno_alloc((size_t)contigs[c].length, sizeof(*contigs[c].label));
    int skipped = 0;
    int used = label_cds(contigs, count, argv[arg + 1], &skipped);
    if (!used) anno_fail("No valid CDS annotations found in %s", argv[arg + 1]);
    enum { PRETRAIN_CHUNK_SIZE = 1 << 20 };
    int halo = wavelets.max_width / 2;
    for (int c = 0; c < count; c++) {
      const Contig *contig = contigs + c;
      for (int start = 0; start < contig->length;) {
        int chunk_length = contig->length - start > PRETRAIN_CHUNK_SIZE
          ? PRETRAIN_CHUNK_SIZE : contig->length - start;
        int source_start = start > halo ? start - halo : 0;
        int core_end = start + chunk_length;
        int source_end = contig->length - core_end > halo ? core_end + halo : contig->length;
        int source_length = source_end - source_start;
        float *chunk_cwt = anno_alloc((size_t)source_length * CWT_CHANNELS, sizeof(*chunk_cwt));
        cwt_transform(&wavelets, contig->seq + source_start, source_length, 1, chunk_cwt);
        gff_cwt_accumulate_range(
          &accumulator, contig, start, chunk_length,
          chunk_cwt + (size_t)(start - source_start) * CWT_CHANNELS);
        free(chunk_cwt);
        start = core_end;
      }
    }
    fprintf(stderr, "Accumulated %s / %s: %d contigs, %d transcripts, %d skipped\n",
            argv[arg], argv[arg + 1], count, used, skipped);
    free_contigs(contigs, count);
  }
  GffCwtModel model;
  uint64_t valid_bases = accumulator.valid_bases;
  gff_cwt_accumulator_finish(&accumulator, &model);
  write_gff_cwt_header(argv[2], &model);
  fprintf(stderr, "Embedded GFF->CWT model trained from %llu valid bases: %s\n",
          (unsigned long long)valid_bases, argv[2]);
  wavelets_free(&wavelets);
  return 0;
}

static int usage(const char *program, int status) {
  fprintf(stderr,
          "Usage:\n"
          "  %s pretrain-cwt <output-header.h> <fasta> <gff3> [<fasta> <gff3> ...]\n"
          "  %s [-m model] [-i initial_model] [-e iterations] [-p threads] -t <genome.fa[.gz]> > genes.gff3\n"
          "  %s [-m model] <genome.fa[.gz]> > genes.gff3\n",
          program, program, program);
  return status;
}

static int parse_threads(const char *value) {
  char *end;
  errno = 0;
  long count = strtol(value, &end, 10);
  if (errno || end == value || *end || count < 1 || count > INT_MAX)
    anno_fail("Invalid thread count: %s", value);
  return (int)count;
}

int main(int argc, char **argv) {
  if (argc > 1 && !strcmp(argv[1], "pretrain-cwt"))
    return pretrain_cwt(argc, argv);
  const char *model = "anno.model";
  const char *initial_model = NULL;
  int iterations = DEFAULT_ITERATIONS, threads = 1, train = 0, option;
  ketopt_t options = KETOPT_INIT;
  while ((option = ketopt(&options, argc, argv, 0, "m:i:e:p:th", NULL)) >= 0) {
    if (option == 'h') return usage(argv[0], 0);
    if (option == 'p') threads = parse_threads(options.arg);
    else if (option == 'm') model = options.arg;
    else if (option == 'i') initial_model = options.arg;
    else if (option == 'e' && atoi(options.arg) > 0) iterations = atoi(options.arg);
    else if (option == 't') train = 1;
    else return usage(argv[0], 1);
  }
  int rest = argc - options.ind;
  char **args = argv + options.ind;
  if (rest != 1)
    return usage(argv[0], 1);
  if (initial_model && !train)
    anno_fail("Initial model is only supported for training");

  crf_init();
  double *weights = anno_alloc(WEIGHTS, sizeof(*weights));
  if (train && initial_model) {
    load_model(initial_model, weights);
    fprintf(stderr, "Warm start from %s\n", initial_model);
  } else if (!train) {
    load_model(model, weights);
  }

  int count;
  Contig *contigs = read_fasta(args[0], &count);
  if (train) {
    Wavelets wavelets;
    wavelets_init(&wavelets);
    GffCwtModel embedded_gff_cwt;
    memcpy(&embedded_gff_cwt, embedded_gff_cwt_values, sizeof(embedded_gff_cwt));
    crf_set_gff_cwt(&embedded_gff_cwt);
    for (int c = 0; c < count; c++) {
      contigs[c].cwt = anno_alloc((size_t)contigs[c].length * CWT_CHANNELS, sizeof(float));
      cwt_transform(&wavelets, contigs[c].seq, contigs[c].length, threads, contigs[c].cwt);
    }
    fprintf(stderr, "Self-training on %d FASTA contigs without GFF labels\n", count);
    crf_train(contigs, count, weights, iterations, threads);
    save_model(model, weights);
    wavelets_free(&wavelets);
  }

  puts("##gff-version 3");
  unsigned long genes = 0;
  for (int c = 0; c < count; c++) {
    uint8_t *path = anno_alloc((size_t)contigs[c].length, 1);
    crf_decode(weights, contigs + c, path);
    write_gff(contigs + c, path, &genes);
    free(path);
    fprintf(stderr, "Annotated %s: %d bp\n", contigs[c].name, contigs[c].length);
  }
  if (fflush(stdout) || ferror(stdout)) anno_fail("Cannot write GFF");
  fprintf(stderr, "Predicted %lu genes\n", genes);
  free_contigs(contigs, count);
  free(weights);
  return 0;
}
