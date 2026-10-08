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

KSEQ_INIT(gzFile, gzread)

const int wave_sizes[WAVE_COUNT] = {4, 5, 6, 7, 8, 9};
enum { CWT_BLOCK = 1024, K2 = 17, K4 = 257, KNUC = 5, SPLICE_RADIUS = 6,
       FLANKS = SPLICE_RADIUS, EMISSION_FEATURES = 2, TRANSITION_FEATURES = 2 + FLANKS };
enum { OFF_ENUC = 0, OFF_E4 = OFF_ENUC + STATES * KNUC,
       OFF_EP = OFF_E4 + STATES * K4, OFF_ED = OFF_EP + STATES * WAVE_COUNT,
       OFF_TNUC = OFF_ED + STATES * WAVE_COUNT,
       OFF_T4 = OFF_TNUC + STATES * STATES * KNUC,
       OFF_TP = OFF_T4 + STATES * STATES * K4, OFF_TD = OFF_TP + STATES * STATES * WAVE_COUNT,
       OFF_SPLICE = OFF_TD + STATES * STATES * WAVE_COUNT,
       WEIGHTS = OFF_SPLICE + 4 * FLANKS * K2 };

typedef struct {
  uint16_t k4[2];
  float cwt_power[WAVE_COUNT];
  int16_t derivative[WAVE_COUNT];
  uint8_t nuc, flank[2][FLANKS];
} Site;

enum { HIDDEN_STATES = 47 };
static const int hidden_offsets[STATES] = {0, 1, 4, 8, 10, 11, 14, 18, 20, 28, 33, 41, 46};
static const int hidden_counts[STATES] = {1, 3, 4, 2, 1, 3, 4, 2, 8, 5, 8, 5, 1};
static int trans_a[STATES * STATES], trans_b[STATES * STATES], ntrans;
static uint8_t hidden_state[HIDDEN_STATES], hidden_context[HIDDEN_STATES];
typedef struct { uint8_t previous, current; uint16_t transition; } HiddenEdge;
static HiddenEdge hidden_edges[KNUC][HIDDEN_STATES * 4];
static int hidden_edge_count[KNUC];
static uint8_t initial_states[KNUC][HIDDEN_STATES];

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
    if (radius_value > (INT_MAX - 1) / 2) anno_fail("CWT scale exceeds supported kernel size");
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

static float *cwt_powers(const Wavelets *wavelets, const char *sequence, int length) {
  float *powers = anno_alloc((size_t)length, WAVE_COUNT * sizeof(*powers));
  double *features = anno_alloc((size_t)CWT_BLOCK * CWT_CHANNELS, sizeof(*features));
  for (int start = 0; start < length; start += CWT_BLOCK) {
    int count = length - start < CWT_BLOCK ? length - start : CWT_BLOCK;
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

static int kmer(const char *seq, int length, int from, int size, int reverse) {
  int code = 0;
  for (int offset = 0; offset < size; offset++) {
    int position = from + (reverse ? size - 1 - offset : offset);
    int base = position < 0 || position >= length ? -1 : base_index(seq[position]);
    if (base < 0) return 1 << (2 * size);
    code = code * 4 + (reverse ? 3 - base : base);
  }
  return code;
}

static Site site_at(const Contig *c, int i) {
  int base = base_index(c->seq[i]);
  Site s = {.k4 = {kmer(c->seq, c->length, i - 2, 4, 0),
                   kmer(c->seq, c->length, i - 2, 4, 1)}, .nuc = base < 0 ? 4 : base};
  for (int offset = 0; offset < FLANKS; offset++) {
    int displacement = -SPLICE_RADIUS + 2 * offset;
    s.flank[0][offset] = kmer(c->seq, c->length, i + displacement, 2, 0);
    s.flank[1][offset] = kmer(c->seq, c->length, i - displacement - 2, 2, 1);
  }
  memcpy(s.cwt_power, c->cwt_power + (size_t)i * WAVE_COUNT, sizeof(s.cwt_power));
  if (i > 0 && i + 1 < c->length && base_index(c->seq[i - 1]) >= 0 &&
      base >= 0 && base_index(c->seq[i + 1]) >= 0) {
    for (int scale = 0; scale < WAVE_COUNT; scale++) {
      double left = c->cwt_power[(size_t)(i - 1) * WAVE_COUNT + scale];
      double right = c->cwt_power[(size_t)(i + 1) * WAVE_COUNT + scale];
      double slope = (right - left) / (right + left + 1e-12);
      s.derivative[scale] = (int16_t)lrint(slope * INT16_MAX);
    }
  }
  return s;
}

static int boundary_kind(int previous, int current) {
  if ((previous >= 1 && previous <= 3 && current >= 4 && current <= 6) ||
      (previous >= 10 && current >= 7 && current <= 9)) return 0;
  if ((previous >= 4 && previous <= 6 && current >= 1 && current <= 3) ||
      (previous >= 7 && previous <= 9 && current >= 10)) return 1;
  if ((previous == 0 && current == 1) || (previous == 7 && current == 0)) return 2;
  if ((previous == 3 && current == 0) || (previous == 0 && current == 9)) return 3;
  return -1;
}

static int allowed(int a, int b) {
  if (a == 0) return b == 0 || b == 1 || b == 9;
  if (a <= 3) return b == 1 + a % 3 || b == 4 + a % 3 || (a == 3 && b == 0);
  if (a <= 6) return b == a || b == a - 3;
  if (a <= 9) return b == 7 + (a + 1) % 3 || b == 10 + (a + 1) % 3 || (a == 7 && b == 0);
  return b == a || b == a - 3;
}

/* Prefix classes carry start/stop recognition through introns without adding learned features. */
static int grammar_next(int previous, int current, int base) {
  int state = hidden_state[previous], context = hidden_context[previous];
  if (!allowed(state, current)) return -1;
  if (!state) {
    if (!current) return 0;
    if (current == 1) return base == 0 ? hidden_offsets[1] + 2 : -1;
    return base == 3 ? hidden_offsets[9] + 3 : base == 1 ? hidden_offsets[9] + 4 : -1;
  }
  if (!current)
    return (state == 3 || state == 7) && context == 1 ? 0 : -1;
  if (state <= 6) {
    if (current >= 4) {
      if (state == 3 && context == 1) return -1;
      return hidden_offsets[current] + (current == 4 ? 0 : context);
    }
    if (current == 1) {
      if (state == 3 && context == 1) return -1;
      return hidden_offsets[1] + (base == 3);
    }
    if (current == 2) {
      if (context == 2) return base == 3 ? hidden_offsets[2] + 3 : -1;
      int next_context = context == 1 && base == 0 ? 1 : context == 1 && base == 2 ? 2 : 0;
      return hidden_offsets[2] + next_context;
    }
    if (context == 3) return base == 2 ? hidden_offsets[3] : -1;
    int stop = (context == 1 && (base == 0 || base == 2)) || (context == 2 && base == 0);
    return hidden_offsets[3] + stop;
  }
  if (current >= 10) return hidden_offsets[current] + (current == 12 ? 0 : context);
  if (current == 9)
    return hidden_offsets[9] + (base == 3 ? 1 : base == 1 ? 2 : 0);
  if (current == 8) {
    if (context == 3)
      return base == 3 ? hidden_offsets[8] + 5 : base == 1 ? hidden_offsets[8] + 7 : -1;
    if (context == 4) return base == 3 ? hidden_offsets[8] + 6 : -1;
    int next_context = context == 1 && base == 3 ? 1 : context == 2 && base == 3 ? 2 :
                       context == 1 && base == 1 ? 3 : context == 2 && base == 0 ? 4 : 0;
    return hidden_offsets[8] + next_context;
  }
  if (context >= 5) return base == 0 ? hidden_offsets[7] : -1;
  if (context >= 1 && context <= 3 && base == 0) return -1;
  return hidden_offsets[7] + (context == 4 && base == 3);
}

static int initial_allowed(int base, int hidden, int full_start) {
  return full_start ? grammar_next(0, hidden_state[hidden], base) == hidden : initial_states[base][hidden];
}

static int terminal_allowed(int hidden) {
  int state = hidden_state[hidden];
  return !state || ((state == 3 || state == 7) && hidden_context[hidden] == 1);
}

static void init_transitions(void) {
  ntrans = 0;
  for (int a = 0; a < STATES; a++)
    for (int b = 0; b < STATES; b++)
      if (allowed(a, b)) trans_a[ntrans] = a, trans_b[ntrans++] = b;
  for (int state = 0; state < STATES; state++)
    for (int context = 0; context < hidden_counts[state]; context++) {
      int hidden = hidden_offsets[state] + context;
      hidden_state[hidden] = (uint8_t)state;
      hidden_context[hidden] = (uint8_t)context;
    }
  memset(hidden_edge_count, 0, sizeof(hidden_edge_count));
  memset(initial_states, 0, sizeof(initial_states));
  for (int base = 0; base < KNUC; base++)
    for (int previous = 0; previous < HIDDEN_STATES; previous++)
      for (int edge = 0; edge < ntrans; edge++) {
        if (trans_a[edge] != hidden_state[previous]) continue;
        int current = grammar_next(previous, trans_b[edge], base);
        if (current < 0) continue;
        if (hidden_edge_count[base] >= HIDDEN_STATES * 4) anno_fail("Too many grammar transitions");
        hidden_edges[base][hidden_edge_count[base]++] = (HiddenEdge){(uint8_t)previous, (uint8_t)current, (uint16_t)edge};
        initial_states[base][current] = 1;
      }
}

static void emission_indices(const Site *site, int state, int *indices) {
  indices[0] = OFF_ENUC + state * KNUC + site->nuc;
  indices[1] = OFF_E4 + state * K4 + site->k4[state >= 7];
}

static int transition_indices(const Site *site, int previous, int current, int *indices) {
  int pair = previous * STATES + current, count = 2;
  int reverse = previous >= 7 || current >= 7, kind = boundary_kind(previous, current);
  indices[0] = OFF_TNUC + pair * KNUC + site->nuc;
  indices[1] = OFF_T4 + pair * K4 + site->k4[reverse];
  if (kind >= 0)
    for (int offset = 0; offset < FLANKS; offset++)
      indices[count++] = OFF_SPLICE + (kind * FLANKS + offset) * K2 + site->flank[reverse][offset];
  return count;
}

static double feature_score(const float *weights, const int *indices, int count) {
  double score = 0;
  for (int index = 0; index < count; index++) score += weights[indices[index]];
  return score;
}

static float derivative_value(const Site *site, int scale) {
  return site->derivative[scale] / (float)INT16_MAX;
}

static void add_features(float *gradient, const int *indices, int count, float value) {
  for (int index = 0; index < count; index++) gradient[indices[index]] += value;
}

static double emission(const float *weights, const Site *site, int state) {
  int indices[EMISSION_FEATURES];
  emission_indices(site, state, indices);
  double score = feature_score(weights, indices, EMISSION_FEATURES);
  for (int scale = 0; scale < WAVE_COUNT; scale++) {
    score += (double)weights[OFF_EP + state * WAVE_COUNT + scale] * site->cwt_power[scale];
    score += (double)weights[OFF_ED + state * WAVE_COUNT + scale] * derivative_value(site, scale);
  }
  return score;
}

static double transition(const float *weights, const Site *site, int previous, int current) {
  int indices[TRANSITION_FEATURES];
  int count = transition_indices(site, previous, current, indices);
  int pair = previous * STATES + current;
  double score = feature_score(weights, indices, count);
  for (int scale = 0; scale < WAVE_COUNT; scale++) {
    score += (double)weights[OFF_TP + pair * WAVE_COUNT + scale] * site->cwt_power[scale];
    score += (double)weights[OFF_TD + pair * WAVE_COUNT + scale] * derivative_value(site, scale);
  }
  return score;
}

static void add_emission(float *gradient, const Site *site, int state, float value) {
  int indices[EMISSION_FEATURES];
  emission_indices(site, state, indices);
  add_features(gradient, indices, EMISSION_FEATURES, value);
  for (int scale = 0; scale < WAVE_COUNT; scale++) {
    gradient[OFF_EP + state * WAVE_COUNT + scale] += value * site->cwt_power[scale];
    gradient[OFF_ED + state * WAVE_COUNT + scale] += value * derivative_value(site, scale);
  }
}

static void add_transition(float *gradient, const Site *site, int previous, int current, float value) {
  int indices[TRANSITION_FEATURES];
  int count = transition_indices(site, previous, current, indices);
  add_features(gradient, indices, count, value);
  int pair = previous * STATES + current;
  for (int scale = 0; scale < WAVE_COUNT; scale++) {
    gradient[OFF_TP + pair * WAVE_COUNT + scale] += value * site->cwt_power[scale];
    gradient[OFF_TD + pair * WAVE_COUNT + scale] += value * derivative_value(site, scale);
  }
}

Contig *read_fasta(const char *path, int *count) {
  gzFile input = gzopen(!strcmp(path, "-") ? "/dev/stdin" : path, "rb");
  if (!input) anno_fail("Cannot open FASTA: %s", path);
  kseq_t *record = kseq_init(input);
  Contig *contigs = NULL;
  *count = 0;
  while (kseq_read(record) >= 0) {
    if (!record->seq.l) continue;
    if (record->seq.l > INT_MAX - CWT_BLOCK) anno_fail("Contig too long: %s", record->name.s);
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

typedef struct { Site *sites; double *emis, *psi, *alpha, *beta, *scale; } Work;

static double log_add(double left, double right) {
  if (left == -INFINITY) return right;
  if (right == -INFINITY) return left;
  double maximum = fmax(left, right);
  return maximum + log1p(exp(-fabs(left - right)));
}

/* Normalized log-space forward-backward on [s, e). */
static double window_posteriors(const float *w, const Contig *c, int s, int e, Work *k) {
  int len = e - s;
  double logz = 0;
  for (int i = 0; i < len; i++) {
    k->sites[i] = site_at(c, s + i);
    for (int st = 0; st < STATES; st++) {
      k->emis[i * STATES + st] = emission(w, k->sites + i, st);
      if (!isfinite(k->emis[i * STATES + st])) anno_fail("Non-finite CRF emission at position %d", s + i + 1);
    }
  }
  for (int i = 0; i < len; i++) {
    const Site *site = k->sites + i;
    const double *em = k->emis + i * STATES;
    double *al = k->alpha + (size_t)i * HIDDEN_STATES, *psi = k->psi + (size_t)i * ntrans;
    if (!i) {
      for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
        al[hidden] = initial_allowed(site->nuc, hidden, s == 0) ? em[hidden_state[hidden]] : -INFINITY;
    } else {
      for (int t = 0; t < ntrans; t++) {
        psi[t] = transition(w, site, trans_a[t], trans_b[t]) + em[trans_b[t]];
        if (!isfinite(psi[t])) anno_fail("Non-finite CRF transition at position %d", s + i + 1);
      }
      for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) al[hidden] = -INFINITY;
      const HiddenEdge *edges = hidden_edges[site->nuc];
      for (int edge = 0; edge < hidden_edge_count[site->nuc]; edge++) {
        const HiddenEdge *step = edges + edge;
        al[step->current] = log_add(al[step->current], al[(int)step->previous - HIDDEN_STATES] + psi[step->transition]);
      }
    }
    if (e == c->length && i == len - 1)
      for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
        if (!terminal_allowed(hidden)) al[hidden] = -INFINITY;
    double total = -INFINITY;
    for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) total = log_add(total, al[hidden]);
    if (!isfinite(total)) anno_fail("Invalid CRF partition at position %d", s + i + 1);
    for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) al[hidden] -= total;
    k->scale[i] = total;
    logz += total;
  }
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
    k->beta[(size_t)(len - 1) * HIDDEN_STATES + hidden] = e < c->length || terminal_allowed(hidden) ? 0 : -INFINITY;
  for (int i = len - 2; i >= 0; i--) {
    double *be = k->beta + (size_t)i * HIDDEN_STATES, *next = be + HIDDEN_STATES, *psi = k->psi + (size_t)(i + 1) * ntrans;
    int base = k->sites[i + 1].nuc;
    const HiddenEdge *edges = hidden_edges[base];
    for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) be[hidden] = -INFINITY;
    for (int edge = 0; edge < hidden_edge_count[base]; edge++) {
      const HiddenEdge *step = edges + edge;
      be[step->previous] = log_add(be[step->previous], psi[step->transition] + next[step->current]);
    }
    for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) be[hidden] -= k->scale[i + 1];
  }
  return logz;
}

static double gold_grammar_mass(const Contig *contig, int start, int end) {
  double previous[HIDDEN_STATES], current[HIDDEN_STATES];
  int base = base_index(contig->seq[start]);
  if (base < 0) base = 4;
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
    previous[hidden] = hidden_state[hidden] == contig->label[start] && initial_allowed(base, hidden, start == 0) ? 0 : -INFINITY;
  for (int position = start + 1; position < end; position++) {
    base = base_index(contig->seq[position]);
    if (base < 0) base = 4;
    for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) current[hidden] = -INFINITY;
    const HiddenEdge *edges = hidden_edges[base];
    for (int edge = 0; edge < hidden_edge_count[base]; edge++) {
      const HiddenEdge *step = edges + edge;
      if (hidden_state[step->current] == contig->label[position])
        current[step->current] = log_add(current[step->current], previous[step->previous]);
    }
    memcpy(previous, current, sizeof(previous));
  }
  double mass = -INFINITY;
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
    if (end < contig->length || terminal_allowed(hidden)) mass = log_add(mass, previous[hidden]);
  return mass;
}

/* Accumulates observed minus expected features for one window.
   Returns the window log-likelihood or NAN when the gold path violates the grammar. */
static double window_gradient(const float *w, const Contig *c, int s, int e, Work *k, float *gradient) {
  int len = e - s;
  double mass = gold_grammar_mass(c, s, e);
  if (mass == -INFINITY) return NAN;
  double gold = mass - window_posteriors(w, c, s, e, k);
  for (int i = 0; i < len; i++) {
    const Site *site = k->sites + i;
    const double *al = k->alpha + (size_t)i * HIDDEN_STATES, *be = k->beta + (size_t)i * HIDDEN_STATES, *psi = k->psi + (size_t)i * ntrans;
    int b = c->label[s + i];
    gold += k->emis[i * STATES + b];
    add_emission(gradient, site, b, 1);
    double state_marginals[STATES] = {0};
    for (int hidden = 0; hidden < HIDDEN_STATES; hidden++)
      state_marginals[hidden_state[hidden]] += exp(al[hidden] + be[hidden]);
    for (int st = 0; st < STATES; st++) add_emission(gradient, site, st, -(float)state_marginals[st]);
    if (!i) continue;
    gold += transition(w, site, c->label[s + i - 1], b);
    add_transition(gradient, site, c->label[s + i - 1], b, 1);
    double edge_marginals[STATES * STATES] = {0};
    const HiddenEdge *edges = hidden_edges[site->nuc];
    for (int edge = 0; edge < hidden_edge_count[site->nuc]; edge++) {
      const HiddenEdge *step = edges + edge;
      edge_marginals[step->transition] += exp(al[(int)step->previous - HIDDEN_STATES] + psi[step->transition] + be[step->current] - k->scale[i]);
    }
    for (int t = 0; t < ntrans; t++)
      if (edge_marginals[t] > 0) add_transition(gradient, site, trans_a[t], trans_b[t], -(float)edge_marginals[t]);
  }
  return gold;
}

static int window_count(const Contig *c) { return (c->length + WINDOW - 1) / WINDOW; }

static Work work_alloc(void) {
  return (Work){
    .sites = anno_alloc(WINDOW, sizeof(Site)),
    .emis = anno_alloc((size_t)WINDOW * STATES, sizeof(double)),
    .psi = anno_alloc((size_t)WINDOW * ntrans, sizeof(double)),
    .alpha = anno_alloc((size_t)WINDOW * HIDDEN_STATES, sizeof(double)),
    .beta = anno_alloc((size_t)WINDOW * HIDDEN_STATES, sizeof(double)),
    .scale = anno_alloc(WINDOW, sizeof(double))
  };
}

static void work_free(Work *k) {
  free(k->sites), free(k->emis), free(k->psi), free(k->alpha), free(k->beta), free(k->scale);
}

void crf_train(Contig *contigs, int count, float *w, int epochs) {
  int windows = 0;
  for (int c = 0; c < count; c++) windows += window_count(contigs + c);
  int *order = anno_alloc(windows, sizeof(*order));
  for (int i = 0; i < windows; i++) order[i] = i;
  Work work = work_alloc();
  float *gradient = anno_alloc(WEIGHTS, sizeof(*gradient));
  double *history = anno_alloc(WEIGHTS, sizeof(*history));
  double *average = anno_alloc(WEIGHTS, sizeof(*average));
  long averaged = 0;
  srand(1);
  for (int epoch = 1; epoch <= epochs; epoch++) {
    double likelihood = 0;
    long positions = 0;
    int skipped = 0;
    for (int i = windows - 1; i > 0; i--) {
      int j = rand() % (i + 1), t = order[i];
      order[i] = order[j], order[j] = t;
    }
    for (int window = 0; window < windows; window++) {
      int contig_index = 0, window_index = order[window];
      while (window_index >= window_count(contigs + contig_index))
        window_index -= window_count(contigs + contig_index++);
      Contig *sequence = contigs + contig_index;
      int start = window_index * WINDOW;
      int end = start + WINDOW < sequence->length ? start + WINDOW : sequence->length;
      int length = end - start;
      memset(gradient, 0, WEIGHTS * sizeof(*gradient));
      double score = window_gradient(w, sequence, start, end, &work, gradient);
      if (isnan(score)) {
        skipped++;
        continue;
      }
      if (!isfinite(score)) anno_fail("Non-finite training likelihood");
      likelihood += score;
      positions += length;
      for (int feature = 0; feature < WEIGHTS; feature++) {
        if (!isfinite(gradient[feature])) anno_fail("Non-finite training gradient");
        float normalized = gradient[feature] / length;
        if (normalized == 0) continue;
        history[feature] += (double)normalized * normalized;
        w[feature] += LEARNING_RATE * normalized / sqrt(history[feature] + 1e-8);
        if (!isfinite(w[feature])) anno_fail("Non-finite training weight");
      }
      if (epoch == epochs) {
        for (int feature = 0; feature < WEIGHTS; feature++) average[feature] += w[feature];
        averaged++;
      }
    }
    fprintf(stderr, "Epoch %d: %d windows trained, %d skipped, log-likelihood/bp %.4f\n",
            epoch, windows - skipped, skipped, positions ? likelihood / positions : 0.0);
  }
  for (int f = 0; averaged && f < WEIGHTS; f++) w[f] = (float)(average[f] / averaged);
  work_free(&work);
  free(order), free(gradient), free(history), free(average);
}

static void viterbi_step(const float *weights, const Contig *contig, int position,
                         const double *previous, double *current, uint8_t *back) {
  Site site = site_at(contig, position);
  double emissions[STATES], transitions[STATES * STATES];
  for (int state = 0; state < STATES; state++) {
    emissions[state] = emission(weights, &site, state);
    if (!isfinite(emissions[state])) anno_fail("Non-finite prediction emission");
  }
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) {
    current[hidden] = position == 0 && initial_allowed(site.nuc, hidden, 1) ? emissions[hidden_state[hidden]] : -INFINITY;
    if (back) back[hidden] = 0;
  }
  if (position) {
    for (int edge = 0; edge < ntrans; edge++) {
      transitions[edge] = transition(weights, &site, trans_a[edge], trans_b[edge]) + emissions[trans_b[edge]];
      if (!isfinite(transitions[edge])) anno_fail("Non-finite prediction transition");
    }
    const HiddenEdge *edges = hidden_edges[site.nuc];
    for (int edge = 0; edge < hidden_edge_count[site.nuc]; edge++) {
      const HiddenEdge *step = edges + edge;
      double score = previous[step->previous] + transitions[step->transition];
      if (score > current[step->current]) {
        current[step->current] = score;
        if (back) back[step->current] = step->previous;
      }
    }
  }
  double maximum = -INFINITY;
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) maximum = fmax(maximum, current[hidden]);
  if (!isfinite(maximum)) anno_fail("No finite prediction path");
  for (int hidden = 0; hidden < HIDDEN_STATES; hidden++) current[hidden] -= maximum;
}

/* Exact whole-contig MAP; recomputed checkpoint blocks bound traceback memory. */
void crf_decode(const float *w, const Contig *c, uint8_t *path) {
  int length = c->length, blocks = window_count(c);
  double *checkpoints = anno_alloc((size_t)blocks * HIDDEN_STATES, sizeof(*checkpoints));
  int capacity = length < WINDOW ? length : WINDOW;
  uint8_t *back = anno_alloc((size_t)capacity * HIDDEN_STATES, sizeof(*back));
  double previous[HIDDEN_STATES] = {0}, current[HIDDEN_STATES];
  for (int position = 0; position < length; position++) {
    if (position % WINDOW == 0)
      memcpy(checkpoints + (size_t)(position / WINDOW) * HIDDEN_STATES, previous, sizeof(previous));
    viterbi_step(w, c, position, previous, current, NULL);
    memcpy(previous, current, sizeof(previous));
  }
  int best = 0;
  for (int hidden = 1; hidden < HIDDEN_STATES; hidden++)
    if (terminal_allowed(hidden) && previous[hidden] > previous[best]) best = hidden;
  if (!isfinite(previous[best])) anno_fail("No complete prediction path");
  for (int block = blocks - 1; block >= 0; block--) {
    int start = block * WINDOW, end = length - start < WINDOW ? length : start + WINDOW;
    memcpy(previous, checkpoints + (size_t)block * HIDDEN_STATES, sizeof(previous));
    for (int position = start; position < end; position++) {
      viterbi_step(w, c, position, previous, current, back + (size_t)(position - start) * HIDDEN_STATES);
      memcpy(previous, current, sizeof(previous));
    }
    for (int position = end - 1; position >= start; position--) {
      path[position] = hidden_state[best];
      best = back[(size_t)(position - start) * HIDDEN_STATES + best];
    }
  }
  free(checkpoints);
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
  if (!file || fwrite("ANNOCRFD", 1, 8, file) != 8 ||
      fwrite(w, sizeof(*w), WEIGHTS, file) != WEIGHTS || fclose(file))
    anno_fail("Cannot write model: %s", path);
}

static void load_model(const char *path, float *w) {
  FILE *file = fopen(path, "rb");
  char magic[8];
  if (!file || fread(magic, 1, 8, file) != 8) anno_fail("Cannot read model: %s", path);
  if (memcmp(magic, "ANNOCRFD", 8)) anno_fail("Incompatible model; retrain with this executable: %s", path);
  if (fread(w, sizeof(*w), WEIGHTS, file) != WEIGHTS || fgetc(file) != EOF)
    anno_fail("Cannot read model: %s", path);
  for (int feature = 0; feature < WEIGHTS; feature++)
    if (!isfinite(w[feature])) anno_fail("Non-finite model weight: %s", path);
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
  fputs("# anno-cwt-v2\n# scales\t", index);
  for (int scale = 0; scale < WAVE_COUNT; scale++) fprintf(index, "%s%d", scale ? "," : "", wave_sizes[scale]);
  fputs("\nseqid\tlength\tplus_offset\tminus_offset\n", index);
  double *features = anno_alloc(CWT_BLOCK * CWT_CHANNELS, sizeof(*features));
  float *values = anno_alloc(CWT_BLOCK * CWT_CHANNELS, sizeof(*values));
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
      for (int start = 0; start < g->length; start += CWT_BLOCK) {
        int n = g->length - start < CWT_BLOCK ? g->length - start : CWT_BLOCK;
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
  fprintf(stderr, "Usage: %s [-m model] [-i initial_model] [-e epochs] <genome.fa[.gz]> [<train.gff3[.gz]>] > genes.gff3\n"
                  "       %s cwt <genome.fa[.gz]> <new_directory>\n", program, program);
  return status;
}

int main(int argc, char **argv) {
  const char *model = "anno.model";
  const char *initial_model = NULL;
  int epochs = DEFAULT_EPOCHS, option;
  ketopt_t options = KETOPT_INIT;
  while ((option = ketopt(&options, argc, argv, 0, "m:i:e:h", NULL)) >= 0) {
    if (option == 'm') model = options.arg;
    else if (option == 'i') initial_model = options.arg;
    else if (option == 'e' && atoi(options.arg) > 0) epochs = atoi(options.arg);
    else return usage(argv[0], option != 'h');
  }
  int rest = argc - options.ind;
  char **args = argv + options.ind;
  int export = rest == 3 && !strcmp(args[0], "cwt");
  if (!export && (rest < 1 || rest > 2)) return usage(argv[0], 1);
  if (initial_model && rest != 2) anno_fail("Initial model is only supported for training");
  init_transitions();
  Wavelets wavelets;
  wavelets_init(&wavelets);
  if (export) {
    export_cwt(&wavelets, args[1], args[2]);
    wavelets_free(&wavelets);
    return 0;
  }
  float *weights = anno_alloc(WEIGHTS, sizeof(*weights));
  if (rest == 1) load_model(model, weights);
  if (initial_model) {
    load_model(initial_model, weights);
    fprintf(stderr, "Warm start from %s; optimizer accumulators are reset\n", initial_model);
  }
  int count;
  Contig *contigs = read_fasta(args[0], &count);
  for (int c = 0; c < count; c++)
    contigs[c].cwt_power = cwt_powers(&wavelets, contigs[c].seq, contigs[c].length);
  if (rest == 2) {
    for (int c = 0; c < count; c++) contigs[c].label = anno_alloc(contigs[c].length, 1);
    int skipped, used = label_cds(contigs, count, args[1], &skipped);
    fprintf(stderr, "Training on %d transcripts (%d overlapping or inconsistent skipped)\n", used, skipped);
    crf_train(contigs, count, weights, epochs);
    save_model(model, weights);
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
  wavelets_free(&wavelets);
  free_contigs(contigs, count);
  free(weights);
  return 0;
}
