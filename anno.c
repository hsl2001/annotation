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

KSEQ_INIT(gzFile, gzread)

const int wave_sizes[WAVE_COUNT] = {4, 5, 8, 9};
static const float power_edges[7] = {0.25f, 0.5f, 0.75f, 1.0f, 1.25f, 1.5f, 2.0f};

/* Weight layout: transition x boundary 4-mer, transition x (peaks at i-1, peaks at i),
   state x nucleotide, state x oriented hexamer, state x per-scale CWT code. */
enum { K4 = 257, KPK = (WAVE_COUNT + 1) * (WAVE_COUNT + 1), KNUC = 5, KHEX = 4097, KCWT = 16 * WAVE_COUNT };
enum { OFF_T4 = 0, OFF_TPK = OFF_T4 + STATES * STATES * K4, OFF_NUC = OFF_TPK + STATES * STATES * KPK,
       OFF_HEX = OFF_NUC + STATES * KNUC, OFF_CWT = OFF_HEX + STATES * KHEX, WEIGHTS = OFF_CWT + STATES * KCWT };

typedef struct { uint16_t k4, hexf, hexr, cwt; uint8_t nuc, pk; } Site;

static int trans_a[STATES * STATES], trans_b[STATES * STATES], ntrans, preds[STATES][4], npreds[STATES];

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

void wavelets_init(Wavelets *wavelets) {
  memset(wavelets, 0, sizeof(*wavelets));
  for (int scale = 0; scale < WAVE_COUNT; scale++) {
    int width = wave_sizes[scale];
    double complex mean = 0.0;
    for (int tap = 0; tap < width; tap++) {
      double coordinate = (tap - width / 2.0 + 0.5) / (width / 8.0);
      double complex value = exp(-0.5 * coordinate * coordinate) * (cexp(I * 6.0 * coordinate) - exp(-18.0));
      wavelets->kernel[scale][tap] = value;
      mean += value / width;
    }
    double energy = 0.0;
    for (int tap = 0; tap < width; tap++) {
      wavelets->kernel[scale][tap] -= mean;
      energy += cabs(wavelets->kernel[scale][tap]) * cabs(wavelets->kernel[scale][tap]);
    }
    for (int tap = 0; tap < width; tap++) wavelets->kernel[scale][tap] = conj(wavelets->kernel[scale][tap]) / sqrt(energy);
  }
}

void cwt_extract(const Wavelets *wavelets, const char *sequence, int length, int start, int count, double *features) {
  memset(features, 0, (size_t)count * CWT_CHANNELS * sizeof(double));
  for (int offset = 0; offset < count; offset++) {
    long position = (long)start + offset;
    if (position < 0 || position >= length) continue;
    for (int scale = 0; scale < WAVE_COUNT; scale++) {
      int width = wave_sizes[scale];
      double complex coefficient = 0.0;
      for (int tap = 0; tap < width; tap++) {
        long context = position + tap - width / 2;
        if (context >= 0 && context < length) coefficient += base_signal(sequence[context]) * wavelets->kernel[scale][tap];
      }
      features[(size_t)offset * CWT_CHANNELS + 2 * scale] = creal(coefficient);
      features[(size_t)offset * CWT_CHANNELS + 2 * scale + 1] = cimag(coefficient);
    }
  }
}

/* Per position and scale: bit 3 = local power maximum, bits 0-2 = power bin. */
uint16_t *cwt_features(const Wavelets *wavelets, const char *sequence, int length) {
  uint16_t *codes = anno_alloc(length, sizeof(*codes));
  double *features = anno_alloc((size_t)1026 * CWT_CHANNELS, sizeof(*features));
  for (int start = 0; start < length; start += 1024) {
    int count = length - start < 1024 ? length - start : 1024;
    cwt_extract(wavelets, sequence, length, start - 1, count + 2, features);
    for (int offset = 0; offset < count; offset++) {
      if (base_index(sequence[start + offset]) < 0) continue;
      uint16_t code = 0;
      for (int scale = 0; scale < WAVE_COUNT; scale++) {
        double power[3];
        for (int tap = 0; tap < 3; tap++) {
          const double *value = features + (size_t)(offset + tap) * CWT_CHANNELS + 2 * scale;
          power[tap] = value[0] * value[0] + value[1] * value[1];
        }
        int bin = 0;
        for (int edge = 0; edge < 7; edge++) bin += power[1] > power_edges[edge];
        int peak = power[1] > 0 && power[1] >= power[0] && power[1] >= power[2];
        code |= (uint16_t)((peak << 3 | bin) << (4 * scale));
      }
      codes[start + offset] = code;
    }
  }
  free(features);
  return codes;
}

static int kmer(const char *seq, int length, int from, int k, int reverse) {
  int code = 0;
  for (int j = 0; j < k; j++) {
    int at = reverse ? from + k - 1 - j : from + j;
    int base = at < 0 || at >= length ? -1 : base_index(seq[at]);
    if (base < 0) return 1 << (2 * k);
    code = code * 4 + (reverse ? 3 - base : base);
  }
  return code;
}

static int peak_count(uint16_t code) {
  int count = 0;
  for (int scale = 0; scale < WAVE_COUNT; scale++) count += code >> (4 * scale + 3) & 1;
  return count;
}

static Site site_at(const Contig *c, int i) {
  int base = base_index(c->seq[i]);
  Site s = {kmer(c->seq, c->length, i - 2, 4, 0), kmer(c->seq, c->length, i - 5, 6, 0),
            kmer(c->seq, c->length, i, 6, 1), c->cwt[i], base < 0 ? 4 : base,
            (i ? peak_count(c->cwt[i - 1]) : 0) * (WAVE_COUNT + 1) + peak_count(c->cwt[i])};
  return s;
}

static int allowed(int a, int b) {
  if (a == 0) return b == 0 || b == 1 || b == 9;
  if (a <= 3) return b == 1 + a % 3 || b == 4 + a % 3 || (a == 3 && b == 0);
  if (a <= 6) return b == a || b == a - 3;
  if (a <= 9) return b == 7 + (a + 1) % 3 || b == 10 + (a + 1) % 3 || (a == 7 && b == 0);
  return b == a || b == a - 3;
}

static int is_stop(const char *seq, int length, int at, int reverse) {
  int code = kmer(seq, length, at, 3, reverse);
  return code == 48 || code == 50 || code == 56; /* TAA TAG TGA */
}

/* In-frame stop codons may only terminate a CDS. */
static int blocked(const char *seq, int length, int i, int a, int b) {
  if (a == 3 && b) return is_stop(seq, length, i - 3, 0);
  if ((a == 7 || a == 12) && b == 9) return is_stop(seq, length, i, 1);
  return 0;
}

static void init_transitions(void) {
  for (int a = 0; a < STATES; a++)
    for (int b = 0; b < STATES; b++)
      if (allowed(a, b)) trans_a[ntrans] = a, trans_b[ntrans++] = b, preds[b][npreds[b]++] = a;
}

static float emission(const float *w, const Site *s, int b) {
  float e = w[OFF_NUC + b * KNUC + s->nuc] + w[OFF_HEX + b * KHEX + (b >= 7 ? s->hexr : s->hexf)];
  for (int k = 0; k < WAVE_COUNT; k++) e += w[OFF_CWT + b * KCWT + k * 16 + (s->cwt >> 4 * k & 15)];
  return e;
}

static float transition(const float *w, const Site *s, int a, int b) {
  return w[OFF_T4 + (a * STATES + b) * K4 + s->k4] + w[OFF_TPK + (a * STATES + b) * KPK + s->pk];
}

static void add_emission(float *g, const Site *s, int b, float v) {
  g[OFF_NUC + b * KNUC + s->nuc] += v;
  g[OFF_HEX + b * KHEX + (b >= 7 ? s->hexr : s->hexf)] += v;
  for (int k = 0; k < WAVE_COUNT; k++) g[OFF_CWT + b * KCWT + k * 16 + (s->cwt >> 4 * k & 15)] += v;
}

static void add_transition(float *g, const Site *s, int a, int b, float v) {
  g[OFF_T4 + (a * STATES + b) * K4 + s->k4] += v;
  g[OFF_TPK + (a * STATES + b) * KPK + s->pk] += v;
}

size_t crf_weights(void) { return WEIGHTS; }

Contig *read_fasta(const char *path, int *count) {
  gzFile input = gzopen(!strcmp(path, "-") ? "/dev/stdin" : path, "rb");
  if (!input) anno_fail("Cannot open FASTA: %s", path);
  kseq_t *record = kseq_init(input);
  Contig *contigs = NULL;
  *count = 0;
  while (kseq_read(record) >= 0) {
    if (!record->seq.l) continue;
    if (record->seq.l > INT_MAX - MAX_WAVE_SIZE) anno_fail("Contig too long: %s", record->name.s);
    contigs = realloc(contigs, (*count + 1) * sizeof(*contigs));
    if (!contigs) anno_fail("Out of memory");
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

typedef struct { int contig, start, end, phase; char strand, *parent; } Cds;
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
    cds = realloc(cds, (n + 1) * sizeof(*cds));
    if (!cds || !(cds[n].parent = strdup(parent))) anno_fail("Out of memory");
    cds[n].contig = contig, cds[n].start = start, cds[n].end = end, cds[n].strand = field[6][0];
    cds[n++].phase = field[7][0] >= '0' && field[7][0] <= '2' ? field[7][0] - '0' : 0;
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
    for (int j = 0; j < m; j++) {
      int pos = (3 - s[j].phase) % 3, len = s[j].end - s[j].start + 1;
      for (int t = 0; t < len; t++, pos = (pos + 1) % 3)
        c->label[minus ? s[j].end - 1 - t : s[j].start - 1 + t] = (minus ? 7 : 1) + pos;
      int q = minus ? (3 - s[j].phase + len - 1) % 3 : (3 - s[j].phase) % 3;
      for (int p = j ? s[j - 1].end : s[j].start; p < s[j].start - 1; p++) c->label[p] = (minus ? 10 : 4) + q;
    }
  }
  for (int i = 0; i < n; i++) free(cds[i].parent);
  free(cds);
  free(groups);
  *skipped = g - used;
  return used;
}

typedef struct { Site *sites; float *emis, *grad; double *psi, *alpha, *beta, *scale; } Work;

/* Scaled forward-backward on window [s, e); leaves alpha, beta, psi, scale in k. Returns log Z. */
static double window_posteriors(const float *w, const Contig *c, int s, int e, Work *k) {
  int len = e - s;
  double logz = 0;
  for (int i = 0; i < len; i++) {
    k->sites[i] = site_at(c, s + i);
    for (int st = 0; st < STATES; st++) k->emis[i * STATES + st] = emission(w, k->sites + i, st);
  }
  for (int i = 0; i < len; i++) {
    const Site *site = k->sites + i;
    const float *em = k->emis + i * STATES;
    double *al = k->alpha + (size_t)i * STATES, *psi = k->psi + (size_t)i * ntrans, shift = -INFINITY;
    if (!i) {
      for (int st = 0; st < STATES; st++) shift = fmax(shift, em[st]);
      for (int st = 0; st < STATES; st++) al[st] = exp(em[st] - shift);
    } else {
      for (int t = 0; t < ntrans; t++) {
        psi[t] = blocked(c->seq, c->length, s + i, trans_a[t], trans_b[t]) ? -INFINITY
               : transition(w, site, trans_a[t], trans_b[t]) + em[trans_b[t]];
        shift = fmax(shift, psi[t]);
      }
      memset(al, 0, STATES * sizeof(*al));
      for (int t = 0; t < ntrans; t++) {
        psi[t] = psi[t] == -INFINITY ? 0 : exp(psi[t] - shift);
        al[trans_b[t]] += al[trans_a[t] - STATES] * psi[t];
      }
    }
    double total = 0;
    for (int st = 0; st < STATES; st++) total += al[st];
    for (int st = 0; st < STATES; st++) al[st] /= total;
    k->scale[i] = total;
    logz += log(total) + shift;
  }
  for (int st = 0; st < STATES; st++) k->beta[(size_t)(len - 1) * STATES + st] = 1;
  for (int i = len - 2; i >= 0; i--) {
    double *be = k->beta + (size_t)i * STATES, *next = be + STATES, *psi = k->psi + (size_t)(i + 1) * ntrans;
    memset(be, 0, STATES * sizeof(*be));
    for (int t = 0; t < ntrans; t++) be[trans_a[t]] += psi[t] * next[trans_b[t]];
    for (int st = 0; st < STATES; st++) be[st] /= k->scale[i + 1];
  }
  return logz;
}

/* Accumulates observed minus expected features for one window.
   Returns the window log-likelihood or NAN when the gold path violates the grammar. */
static double window_gradient(const float *w, const Contig *c, int s, int e, Work *k) {
  int len = e - s;
  for (int i = 1; i < len; i++) {
    int a = c->label[s + i - 1], b = c->label[s + i];
    if (!allowed(a, b) || blocked(c->seq, c->length, s + i, a, b)) return NAN;
  }
  double gold = -window_posteriors(w, c, s, e, k);
  for (int i = 0; i < len; i++) {
    const Site *site = k->sites + i;
    const double *al = k->alpha + (size_t)i * STATES, *be = k->beta + (size_t)i * STATES, *psi = k->psi + (size_t)i * ntrans;
    int b = c->label[s + i];
    gold += k->emis[i * STATES + b];
    add_emission(k->grad, site, b, 1);
    for (int st = 0; st < STATES; st++) add_emission(k->grad, site, st, -(float)(al[st] * be[st]));
    if (!i) continue;
    gold += transition(w, site, c->label[s + i - 1], b);
    add_transition(k->grad, site, c->label[s + i - 1], b, 1);
    for (int t = 0; t < ntrans; t++) {
      double mu = al[trans_a[t] - STATES] * psi[t] * be[trans_b[t]] / k->scale[i];
      if (mu > 0) add_transition(k->grad, site, trans_a[t], trans_b[t], -(float)mu);
    }
  }
  return gold;
}

static int window_count(const Contig *c) { return (c->length + WINDOW - 1) / WINDOW; }

static Work work_alloc(void) {
  Work k = {anno_alloc(WINDOW, sizeof(Site)), anno_alloc((size_t)WINDOW * STATES, sizeof(float)),
            anno_alloc(WEIGHTS, sizeof(float)), anno_alloc((size_t)WINDOW * ntrans, sizeof(double)),
            anno_alloc((size_t)WINDOW * STATES, sizeof(double)), anno_alloc((size_t)WINDOW * STATES, sizeof(double)),
            anno_alloc(WINDOW, sizeof(double))};
  return k;
}

static void work_free(Work *k) {
  free(k->sites), free(k->emis), free(k->grad), free(k->psi), free(k->alpha), free(k->beta), free(k->scale);
}

void crf_train(Contig *contigs, int count, float *w, int epochs) {
  int windows = 0;
  for (int c = 0; c < count; c++) windows += window_count(contigs + c);
  int *order = anno_alloc(windows, sizeof(*order));
  for (int i = 0; i < windows; i++) order[i] = i;
  Work k = work_alloc();
  float *history = anno_alloc(WEIGHTS, sizeof(*history));
  srand(1);
  for (int epoch = 1; epoch <= epochs; epoch++) {
    double likelihood = 0;
    long positions = 0;
    int skipped = 0;
    for (int i = windows - 1; i > 0; i--) {
      int j = rand() % (i + 1), t = order[i];
      order[i] = order[j], order[j] = t;
    }
    for (int n = 0; n < windows; n++) {
      int c = 0, index = order[n];
      while (index >= window_count(contigs + c)) index -= window_count(contigs + c++);
      int s = index * WINDOW, e = s + WINDOW < contigs[c].length ? s + WINDOW : contigs[c].length;
      double value = window_gradient(w, contigs + c, s, e, &k);
      if (isnan(value)) { skipped++; continue; }
      likelihood += value, positions += e - s;
      for (int f = 0; f < WEIGHTS; f++) {
        if (k.grad[f] == 0) continue;
        float g = k.grad[f] / (e - s) - L2_PENALTY * w[f];
        history[f] += g * g;
        w[f] += LEARNING_RATE * g / sqrtf(history[f] + 1e-8f);
        k.grad[f] = 0;
      }
    }
    fprintf(stderr, "Epoch %d: %d windows trained, %d skipped, log-likelihood/bp %.4f\n",
            epoch, windows - skipped, skipped, positions ? likelihood / positions : 0.0);
  }
  free(order), free(history), work_free(&k);
}

/* Viterbi with 2-bit back pointers per state (every state has at most three predecessors). */
void crf_decode(const float *w, const Contig *c, uint8_t *path) {
  int n = c->length;
  uint32_t *back = anno_alloc(n, sizeof(*back));
  double prev[STATES], cur[STATES]; /* chromosome-scale sums exceed float resolution */
  for (int i = 0; i < n; i++) {
    Site site = site_at(c, i);
    uint32_t bits = 0;
    for (int b = 0; b < STATES; b++) {
      double best = i ? -INFINITY : 0;
      int arg = 0;
      for (int j = 0; i && j < npreds[b]; j++) {
        int a = preds[b][j];
        if (blocked(c->seq, n, i, a, b)) continue;
        double value = prev[a] + transition(w, &site, a, b);
        if (value > best) best = value, arg = j;
      }
      cur[b] = best + emission(w, &site, b);
      bits |= (uint32_t)arg << 2 * b;
    }
    back[i] = bits;
    memcpy(prev, cur, sizeof(cur));
  }
  int b = 0;
  for (int st = 1; st < STATES; st++) if (prev[st] > prev[b]) b = st;
  for (int i = n - 1; i >= 0; i--) path[i] = b, b = preds[b][back[i] >> 2 * b & 3];
  free(back);
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

static void save_model(const char *path, const float *w) {
  FILE *file = fopen(path, "wb");
  if (!file || fwrite("ANNOCRF1", 1, 8, file) != 8 || fwrite(w, sizeof(*w), WEIGHTS, file) != WEIGHTS || fclose(file))
    anno_fail("Cannot write model: %s", path);
}

static void load_model(const char *path, float *w) {
  FILE *file = fopen(path, "rb");
  char magic[8];
  if (!file || fread(magic, 1, 8, file) != 8 || memcmp(magic, "ANNOCRF1", 8) ||
      fread(w, sizeof(*w), WEIGHTS, file) != WEIGHTS || fgetc(file) != EOF)
    anno_fail("Cannot read model: %s", path);
  fclose(file);
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
static void export_cwt(const Wavelets *wavelets, const char *fasta, const char *directory) {
  int count;
  Contig *contigs = read_fasta(fasta, &count);
  if (mkdir(directory, 0777)) anno_fail("Cannot create directory %s: %s", directory, strerror(errno));
  FILE *matrix = create_in(directory, "matrix.bin"), *index = create_in(directory, "contigs.tsv");
  fputs("# anno-cwt-v1\n# kernel_widths\t", index);
  for (int scale = 0; scale < WAVE_COUNT; scale++) fprintf(index, "%s%d", scale ? "," : "", wave_sizes[scale]);
  fputs("\nseqid\tlength\tplus_offset\tminus_offset\n", index);
  double *features = anno_alloc(1024 * CWT_CHANNELS, sizeof(*features));
  float *values = anno_alloc(1024 * CWT_CHANNELS, sizeof(*values));
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
      for (int start = 0; start < g->length; start += 1024) {
        int n = g->length - start < 1024 ? g->length - start : 1024;
        cwt_extract(wavelets, bases, g->length, start, n, features);
        for (size_t v = 0; v < (size_t)n * CWT_CHANNELS; v++) values[v] = (float)features[v];
        write_floats(matrix, values, (size_t)n * CWT_CHANNELS);
      }
    }
    free(reverse);
    rows += 2 * (size_t)g->length;
  }
  if (fclose(matrix) || fclose(index)) anno_fail("Cannot finish CWT export");
  fprintf(stderr, "Saved %zu rows x %d complex64 scales\n", rows, WAVE_COUNT);
}

static int usage(const char *program, int status) {
  fprintf(stderr, "Usage: %s [-m model] [-e epochs] <genome.fa[.gz]> [<train.gff3[.gz]>] > genes.gff3\n"
                  "       %s cwt <genome.fa[.gz]> <new_directory>\n", program, program);
  return status;
}

int main(int argc, char **argv) {
  const char *model = "anno.model";
  int epochs = 3, option;
  ketopt_t options = KETOPT_INIT;
  while ((option = ketopt(&options, argc, argv, 0, "m:e:h", NULL)) >= 0) {
    if (option == 'm') model = options.arg;
    else if (option == 'e' && atoi(options.arg) > 0) epochs = atoi(options.arg);
    else return usage(argv[0], option != 'h');
  }
  int rest = argc - options.ind;
  char **args = argv + options.ind;
  int export = rest == 3 && !strcmp(args[0], "cwt");
  if (!export && (rest < 1 || rest > 2)) return usage(argv[0], 1);
  init_transitions();
  Wavelets wavelets;
  wavelets_init(&wavelets);
  if (export) {
    export_cwt(&wavelets, args[1], args[2]);
    return 0;
  }
  int count;
  Contig *contigs = read_fasta(args[0], &count);
  for (int c = 0; c < count; c++) contigs[c].cwt = cwt_features(&wavelets, contigs[c].seq, contigs[c].length);
  float *weights = anno_alloc(WEIGHTS, sizeof(*weights));
  if (rest == 2) {
    for (int c = 0; c < count; c++) contigs[c].label = anno_alloc(contigs[c].length, 1);
    int skipped, used = label_cds(contigs, count, args[1], &skipped);
    fprintf(stderr, "Training on %d transcripts (%d overlapping or inconsistent skipped)\n", used, skipped);
    crf_train(contigs, count, weights, epochs);
    save_model(model, weights);
  } else {
    load_model(model, weights);
  }
  puts("##gff-version 3");
  unsigned long genes = 0;
  for (int c = 0; c < count; c++) {
    uint8_t *path = anno_alloc(contigs[c].length, 1);
    crf_decode(weights, contigs + c, path);
    write_gff(contigs + c, path, &genes);
    free(path);
    fprintf(stderr, "Annotated %s: %d bp\n", contigs[c].name, contigs[c].length);
  }
  if (fflush(stdout) || ferror(stdout)) anno_fail("Cannot write GFF");
  fprintf(stderr, "Predicted %lu genes\n", genes);
  return 0;
}
