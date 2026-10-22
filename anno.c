#ifdef __APPLE__
#define _DARWIN_C_SOURCE
#endif
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
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>

#define RFFT_IMPLEMENTATION
#include "rfft.h"
#define LBFGS_IMPLEMENTATION
#include "lbfgs.h"

const int wave_sizes[WAVE_COUNT] = {CWT_SCALES};
#define CWT_IMPLEMENTATION
#include "cwt.h"
#define CRF_IMPLEMENTATION
#include "crf.h"

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

typedef struct { float *data; size_t size; } CwtStore;

static CwtStore cwt_store(Contig *contigs, int count) {
  size_t bases = 0;
  for (int i = 0; i < count; i++) {
    if ((size_t)contigs[i].genome_length > SIZE_MAX - bases) anno_fail("CWT size overflow");
    bases += (size_t)contigs[i].genome_length;
  }
  if (bases > SIZE_MAX / (CWT_CHANNELS * sizeof(float))) anno_fail("CWT size overflow");
  size_t bytes = bases * CWT_CHANNELS * sizeof(float);
  FILE *backing = tmpfile();
  if (!backing) anno_fail("Cannot create CWT scratch file: %s", strerror(errno));
  void *mapping;
  size_t mapping_size;
  float *data = cwt_map_file(backing, bytes, 0, &mapping, &mapping_size, "CWT scratch");
  if (fclose(backing)) anno_fail("Cannot close CWT scratch file: %s", strerror(errno));
  size_t offset = 0;
  for (int i = 0; i < count; i++) {
    contigs[i].cwt = data + offset * CWT_CHANNELS;
    offset += (size_t)contigs[i].genome_length;
  }
  return (CwtStore){data, mapping_size};
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
    c->name = anno_strdup(record->name.s);
    c->seq = anno_strdup(record->seq.s);
    c->length = (int)record->seq.l;
    c->genome_length = c->length;
    if (record->comment.l) {
      const char *crop = strstr(record->comment.s, "crop=");
      int offset = 0, length = 0;
      if (crop && sscanf(crop + 5, "%d:%d", &offset, &length) == 2 &&
          offset >= 0 && length > 0 && offset <= c->genome_length - length) {
        c->offset = offset;
        c->length = length;
      }
    }
  }
  kseq_destroy(record);
  if (gzclose(input) != Z_OK || !*count) anno_fail("No readable FASTA records: %s", path);
  return contigs;
}

static void free_contigs(Contig *contigs, int count) {
  for (int contig = 0; contig < count; contig++) {
    free(contigs[contig].name);
    free(contigs[contig].seq);
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
    if (contig == count || start < 1 || end < start || end > contigs[contig].genome_length) continue;
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
    int off = c->offset;
    printf("%s\tanno\tgene\t%d\t%d\t.\t%c\t.\tID=gene%lu\n", c->name, off + start + 1, off + i, strand, g);
    printf("%s\tanno\tmRNA\t%d\t%d\t.\t%c\t.\tID=mRNA%lu;Parent=gene%lu\n", c->name, off + start + 1, off + i, strand, g, g);
    for (int j = start; j < i;) {
      if (!is_cds(path[j])) { j++; continue; }
      int k = j;
      while (k < i && is_cds(path[k])) k++;
      int phase = (3 - (minus ? path[k - 1] - 7 : path[j] - 1)) % 3;
      printf("%s\tanno\texon\t%d\t%d\t.\t%c\t.\tParent=mRNA%lu\n", c->name, off + j + 1, off + k, strand, g);
      printf("%s\tanno\tCDS\t%d\t%d\t.\t%c\t%d\tID=cds%lu;Parent=mRNA%lu\n", c->name, off + j + 1, off + k, strand, phase, g, g);
      j = k;
    }
  }
}

typedef struct {
  uint32_t index;
  float weight;
} ModelEntry;

static void save_model(const char *path, const GffCwtModel *gff_cwt, const double *w) {
  ModelEntry *entries = anno_alloc(WEIGHTS, sizeof(*entries));
  uint32_t header[3] = {WEIGHTS, GFF_CWT_BASIS * CWT_CHANNELS, 0};
  for (int i = 0; i < WEIGHTS; i++)
    if ((float)w[i]) entries[header[2]++] = (ModelEntry){(uint32_t)i, (float)w[i]};
  FILE *file = fopen(path, "wb");
  if (!file || fwrite(header, sizeof(*header), 3, file) != 3 ||
      fwrite(gff_cwt, sizeof(*gff_cwt), 1, file) != 1 ||
      fwrite(entries, sizeof(*entries), header[2], file) != header[2] || fclose(file))
    anno_fail("Cannot write model: %s", path);
  free(entries);
  crf_save_weights_text("anno.weights.txt", w);
}

static void load_model(const char *path, GffCwtModel *gff_cwt, double *w) {
  uint32_t header[3];
  FILE *file = fopen(path, "rb");
  if (!file || fread(header, sizeof(*header), 3, file) != 3 ||
      header[0] != WEIGHTS || header[1] != GFF_CWT_BASIS * CWT_CHANNELS ||
      header[2] > WEIGHTS || fread(gff_cwt, sizeof(*gff_cwt), 1, file) != 1)
    anno_fail("Cannot read model: %s (expected %d-weight feature layout)", path, WEIGHTS);
  const float *gc_floats = (const float *)gff_cwt;
  for (size_t i = 0; i < sizeof(*gff_cwt) / sizeof(float); i++)
    if (!isfinite(gc_floats[i])) anno_fail("Non-finite GFF->CWT weight: %s", path);
  ModelEntry *entries = anno_alloc(header[2], sizeof(*entries));
  if (fread(entries, sizeof(*entries), header[2], file) != header[2] || fgetc(file) != EOF)
    anno_fail("Truncated or oversized model: %s", path);
  fclose(file);
  memset(w, 0, WEIGHTS * sizeof(*w));
  for (uint32_t i = 0; i < header[2]; i++) {
    if (entries[i].index >= WEIGHTS) anno_fail("Model weight index out of range: %s", path);
    if (!isfinite(entries[i].weight)) anno_fail("Non-finite model weight: %s", path);
    w[entries[i].index] = entries[i].weight;
  }
  free(entries);
}

static int usage(const char *program, int status) {
  fprintf(stderr, "Usage: %s [-m model] [-i initial_model] [-e iterations] [-p threads] <genome.fa[.gz]> [<train.gff3[.gz]>] > genes.gff3\n",
          program);
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
  const char *model = "anno.model";
  const char *initial_model = NULL;
  int iterations = DEFAULT_ITERATIONS, threads = 1, option;
  ketopt_t options = KETOPT_INIT;
  while ((option = ketopt(&options, argc, argv, 0, "m:i:e:p:h", NULL)) >= 0) {
    if (option == 'h') return usage(argv[0], 0);
    if (option == 'p') threads = parse_threads(options.arg);
    else if (option == 'm') model = options.arg;
    else if (option == 'i') initial_model = options.arg;
    else if (option == 'e' && atoi(options.arg) > 0) iterations = atoi(options.arg);
    else return usage(argv[0], 1);
  }
  int rest = argc - options.ind;
  char **args = argv + options.ind;
  if (rest < 1 || rest > 2)
    return usage(argv[0], 1);
  if (initial_model && rest != 2)
    anno_fail("Initial model is only supported for training");
  Wavelets wavelets;
  wavelets_init(&wavelets);
  crf_init();
  GffCwtModel gff_cwt = {0};
  double *weights = anno_alloc(WEIGHTS, sizeof(*weights));
  if (rest == 1) {
    load_model(model, &gff_cwt, weights);
    crf_set_gff_cwt(&gff_cwt);
  }
  if (initial_model) {
    load_model(initial_model, &gff_cwt, weights);
    crf_set_gff_cwt(&gff_cwt);
    fprintf(stderr, "Warm start from %s\n", initial_model);
  }
  int count;
  Contig *contigs = read_fasta(args[0], &count);
  CwtStore cwt = cwt_store(contigs, count);
  for (int c = 0; c < count; c++)
    cwt_transform(&wavelets, contigs[c].seq, contigs[c].genome_length, threads, contigs[c].cwt, NULL);
  if (rest == 2) {
    for (int c = 0; c < count; c++) contigs[c].label = anno_alloc(contigs[c].genome_length, 1);
    int skipped, used = label_cds(contigs, count, args[1], &skipped);
    fprintf(stderr, "Training on %d transcripts (%d overlapping or inconsistent skipped)\n", used, skipped);
    gff_cwt_train(contigs, count, &gff_cwt);
    crf_set_gff_cwt(&gff_cwt);
    double raw_energy = 0.0;
    double res_energy = gff_cwt_eval_residual(contigs, count, &gff_cwt, &raw_energy);
    fprintf(stderr, "GFF->CWT model trained: raw energy %.6g, residual energy %.6g\n",
            raw_energy, res_energy);
    crf_train(contigs, count, weights, iterations, threads);
    save_model(model, &gff_cwt, weights);
  }
  puts("##gff-version 3");
  unsigned long genes = 0;
  for (int c = 0; c < count; c++) {
    uint8_t *path = anno_alloc(contigs[c].length, 1);
    crf_decode_threads(weights, contigs + c, path, threads);
    write_gff(contigs + c, path, &genes);
    free(path);
    fprintf(stderr, "Annotated %s: %d bp\n", contigs[c].name, contigs[c].length);
  }
  if (fflush(stdout) || ferror(stdout)) anno_fail("Cannot write GFF");
  fprintf(stderr, "Predicted %lu genes\n", genes);
  wavelets_free(&wavelets);
  free_contigs(contigs, count);
  if (munmap(cwt.data, cwt.size)) anno_fail("Cannot release CWT scratch mapping: %s", strerror(errno));
  free(weights);
  return 0;
}
