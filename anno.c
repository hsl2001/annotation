#define _POSIX_C_SOURCE 200809L

#include "anno.h"

#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/stat.h>
#include <zlib.h>

static const char *const class_names[ANNO_CLASSES] = {
  "intergenic", "plus_intron", "plus_five_prime_UTR",
  "plus_three_prime_UTR", "minus_intron", "minus_five_prime_UTR",
  "minus_three_prime_UTR"
};
static const uint32_t cwt_scales[ANNO_CWT_CHANNELS] = {4, 5, 6};
static const char model_magic[8] = {'A', 'N', 'N', 'O', '7', 'S', 'T', '1'};
static const uint32_t model_version = 4;
static const uint32_t model_endian = 0x01020304u;

struct AnnoModel {
  AnnoUnetConfig config;
  AnnoModelMetadata metadata;
  float *weights;
};

typedef struct {
  int64_t length;
  int64_t padded;
  int64_t *sizes;
  float *padded_input;
  float **enc;
  float **down;
  float **up;
  float **dec;
} AnnoForwardCache;

typedef struct {
  char *name;
  char *sequence;
  int64_t length;
  float *cwt;
  uint8_t *labels;
} Sequence;

typedef struct {
  char *parent;
  int sequence;
  int kind;
  char strand;
  int64_t start;
  int64_t end;
} Feature;

typedef struct {
  Feature *items;
  size_t count;
  size_t capacity;
} FeatureList;

typedef struct {
  int64_t start;
  int64_t end;
} Segment;

typedef struct {
  int64_t start;
  int64_t end;
  int64_t transcript_start;
  int64_t transcript_end;
  int64_t length;
  char strand;
} Orf;

typedef struct {
  FILE *files[ANNO_CLASSES];
} BedgraphWriter;

typedef struct {
  const char *model;
  const char *bedgraph;
  int iterations;
  int threads;
} Options;

enum {
  FEATURE_EXON,
  FEATURE_CDS,
  FEATURE_5UTR,
  FEATURE_3UTR,
  FEATURE_INTRON
};

/** Input: const char *format, .... Output: does not return. */
static _Noreturn void fail(const char *format, ...) {
  va_list arguments;
  va_start(arguments, format);
  fputs("Error: ", stderr);
  vfprintf(stderr, format, arguments);
  fputc('\n', stderr);
  va_end(arguments);
  exit(EXIT_FAILURE);
}

/** Input: size_t count, size_t size. Output: void * result. */
static void *checked_calloc(size_t count, size_t size) {
  if (size && count > SIZE_MAX / size) fail("allocation size overflow");
  void *memory = calloc(count ? count : 1, size);
  if (!memory) fail("out of memory");
  return memory;
}

/** Input: void *old, size_t count, size_t size. Output: void * result. */
static void *checked_realloc(void *old, size_t count, size_t size) {
  if (size && count > SIZE_MAX / size) fail("allocation size overflow");
  void *memory = realloc(old, count * size);
  if (!memory) fail("out of memory");
  return memory;
}

/** Input: const char *text. Output: char * result. */
static char *copy_string(const char *text) {
  if (!text) fail("invalid null string");
  char *copy = strdup(text);
  if (!copy) fail("out of memory");
  return copy;
}

/** Input: char base. Output: int result. */
static int base_index(char base) {
  switch (base) {
  case 'A': case 'a': return 0;
  case 'C': case 'c': return 1;
  case 'G': case 'g': return 2;
  case 'T': case 't': return 3;
  default: return -1;
  }
}

/** Input: char base. Output: char result. */
static char complement(char base) {
  static const char complements[] = "TGCA";
  int index = base_index(base);
  return index < 0 ? 'N' : complements[index];
}

/** Input: const char *codon. Output: int result. */
static int is_stop(const char *codon) {
  return (((codon[0] == 'T' || codon[0] == 't') &&
           (codon[1] == 'A' || codon[1] == 'a') &&
           (codon[2] == 'A' || codon[2] == 'a' ||
            codon[2] == 'G' || codon[2] == 'g')) ||
         ((codon[0] == 'T' || codon[0] == 't') &&
          (codon[1] == 'G' || codon[1] == 'g') &&
          (codon[2] == 'A' || codon[2] == 'a')));
}

/** Input: const char *codon. Output: int result. */
static int is_atg(const char *codon) {
  return (codon[0] == 'A' || codon[0] == 'a') &&
         (codon[1] == 'T' || codon[1] == 't') &&
         (codon[2] == 'G' || codon[2] == 'g');
}

/** Input: const Sequence *sequences, int count,                           const char *name. Output: int result. */
static int sequence_index(const Sequence *sequences, int count,
                          const char *name) {
  for (int i = 0; i < count; i++)
    if (!strcmp(sequences[i].name, name)) return i;
  return -1;
}

/** Input: Sequence *sequences, int count. Output: none; output parameters are updated or released. */
static void sequence_free(Sequence *sequences, int count) {
  if (!sequences) return;
  for (int i = 0; i < count; i++) {
    free(sequences[i].name);
    free(sequences[i].sequence);
    free(sequences[i].cwt);
    free(sequences[i].labels);
  }
  free(sequences);
}

/** Input: const char *path, int *count_out. Output: Sequence * result. */
static Sequence *read_fasta(const char *path, int *count_out) {
  gzFile input = gzopen(!strcmp(path, "-") ? "/dev/stdin" : path, "rb");
  if (!input) fail("cannot open FASTA %s: %s", path, strerror(errno));
  Sequence *sequences = NULL;
  int count = 0;
  size_t capacity = 0;
  Sequence *current = NULL;
  for (;;) {
    char buffer[8192];
    char *result = gzgets(input, buffer, sizeof(buffer));
    if (!result) {
      if (gzeof(input)) break;
      fail("cannot read FASTA %s", path);
    }
    size_t length = strlen(buffer);
    while (length && (buffer[length - 1] == '\n' || buffer[length - 1] == '\r'))
      buffer[--length] = '\0';
    if (buffer[0] == '>') {
      char *name = buffer + 1;
      while (*name && isspace((unsigned char)*name)) name++;
      char *end = name;
      while (*end && !isspace((unsigned char)*end)) end++;
      if (!*name || *end) *end = '\0';
      if (!*name) fail("empty FASTA sequence ID");
      if ((size_t)count == capacity) {
        capacity = capacity ? capacity * 2 : 8;
        sequences = checked_realloc(sequences, capacity, sizeof(*sequences));
      }
      current = &sequences[count++];
      *current = (Sequence){.name = copy_string(name)};
      continue;
    }
    if (!current) {
      if (length) fail("FASTA sequence appears before a header");
      continue;
    }
    for (size_t i = 0; i < length; i++)
      if (!isspace((unsigned char)buffer[i])) {
        current->sequence = checked_realloc(current->sequence,
                                            (size_t)current->length + 2, 1);
        current->sequence[current->length++] =
          (char)toupper((unsigned char)buffer[i]);
      }
  }
  if (gzclose(input) != Z_OK || !count) fail("no readable FASTA records in %s", path);
  for (int i = 0; i < count; i++) {
    if (!sequences[i].length) fail("empty FASTA sequence %s", sequences[i].name);
    sequences[i].sequence = checked_realloc(sequences[i].sequence,
                                            (size_t)sequences[i].length + 1, 1);
    sequences[i].sequence[sequences[i].length] = '\0';
  }
  *count_out = count;
  return sequences;
}

/** Input: FeatureList *list, const Feature *feature. Output: none; output parameters are updated or released. */
static void append_feature(FeatureList *list, const Feature *feature) {
  if (list->count == list->capacity) {
    list->capacity = list->capacity ? list->capacity * 2 : 128;
    list->items = checked_realloc(list->items, list->capacity,
                                  sizeof(*list->items));
  }
  list->items[list->count++] = *feature;
}

/** Input: FeatureList *list. Output: none; output parameters are updated or released. */
static void feature_free(FeatureList *list) {
  for (size_t i = 0; i < list->count; i++) {
    free(list->items[i].parent);
  }
  free(list->items);
  *list = (FeatureList){0};
}

/** Input: const char *text. Output: int result. */
static int feature_kind(const char *text) {
  if (!strcasecmp(text, "exon")) return FEATURE_EXON;
  if (!strcasecmp(text, "cds")) return FEATURE_CDS;
  if (!strcasecmp(text, "five_prime_utr") || !strcasecmp(text, "5utr"))
    return FEATURE_5UTR;
  if (!strcasecmp(text, "three_prime_utr") || !strcasecmp(text, "3utr"))
    return FEATURE_3UTR;
  if (!strcasecmp(text, "intron")) return FEATURE_INTRON;
  return -1;
}

/** Input: char *attributes, const char *wanted. Output: char * result. */
static char *attribute_value(char *attributes, const char *wanted) {
  size_t wanted_length = strlen(wanted);
  for (char *field = attributes; field && *field;) {
    char *next = strchr(field, ';');
    if (next) *next++ = '\0';
    while (*field == ' ' || *field == '\t') field++;
    size_t key_length = strcspn(field, "= \t");
    if (key_length == wanted_length && !strncmp(field, wanted, wanted_length)) {
      char *value = field + key_length;
      while (*value == '=' || *value == ' ' || *value == '\t') value++;
      if (*value == '"') value++;
      char *quote = strrchr(value, '"');
      if (quote) *quote = '\0';
      return value;
    }
    field = next;
  }
  return NULL;
}

/** Input: char *attributes. Output: char * result. */
static char *feature_parent(char *attributes) {
  char *parent = attribute_value(attributes, "Parent");
  if (!parent) parent = attribute_value(attributes, "transcript_id");
  if (!parent) parent = attribute_value(attributes, "gene_id");
  if (!parent) return NULL;
  char *separator = strpbrk(parent, ",;");
  if (separator) *separator = '\0';
  return copy_string(parent);
}

/** Input: const char *text, int64_t *value. Output: int result. */
static int parse_integer(const char *text, int64_t *value) {
  char *end = NULL;
  errno = 0;
  long long parsed = strtoll(text, &end, 10);
  if (errno || end == text || *end || parsed < 0) return 0;
  *value = parsed;
  return 1;
}

/** Input: const char *path, const Sequence *sequences,                                  int sequence_count. Output: FeatureList result. */
static FeatureList read_features(const char *path, const Sequence *sequences,
                                 int sequence_count) {
  gzFile input = gzopen(path, "rb");
  if (!input) fail("cannot open annotation %s: %s", path, strerror(errno));
  FeatureList list = {0};
  char line[65536];
  while (gzgets(input, line, sizeof(line))) {
    if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;
    char *fields[9] = {0};
    char *cursor = line;
    int field_count = 0;
    while (field_count < 9 && cursor) {
      fields[field_count++] = cursor;
      cursor = strchr(cursor, '\t');
      if (cursor) *cursor++ = '\0';
    }
    if (field_count != 9) fail("annotation record does not have 9 columns");
    int kind = feature_kind(fields[2]);
    if (kind < 0) continue;
    int sequence = sequence_index(sequences, sequence_count, fields[0]);
    if (sequence < 0) fail("annotation sequence is absent from FASTA: %s", fields[0]);
    int64_t start = 0, end = 0;
    if (!parse_integer(fields[3], &start) || !parse_integer(fields[4], &end) ||
        start < 1 || end < start || end > sequences[sequence].length)
      fail("invalid annotation interval for %s", fields[0]);
    char strand = fields[6][0];
    if (strand != '+' && strand != '-')
      fail("unsupported annotation strand for %s", fields[0]);
    char *parent = feature_parent(fields[8]);
    if (kind >= FEATURE_EXON && !parent)
      fail("annotation feature has no transcript parent");
    Feature feature = {.parent = parent, .sequence = sequence, .kind = kind,
                       .strand = strand, .start = start - 1, .end = end};
    append_feature(&list, &feature);
  }
  if (!gzeof(input)) fail("cannot read annotation %s", path);
  if (gzclose(input) != Z_OK) fail("cannot close annotation %s", path);
  return list;
}

/** Input: const Feature *feature. Output: uint8_t result. */
static uint8_t explicit_label(const Feature *feature) {
  if (feature->kind == FEATURE_INTRON)
    return feature->strand == '+' ? ANNO_PLUS_INTRON : ANNO_MINUS_INTRON;
  if (feature->kind == FEATURE_5UTR)
    return feature->strand == '+' ? ANNO_PLUS_5UTR : ANNO_MINUS_5UTR;
  if (feature->kind == FEATURE_3UTR)
    return feature->strand == '+' ? ANNO_PLUS_3UTR : ANNO_MINUS_3UTR;
  return ANNO_INTERGENIC;
}

/** Input: const Feature *left, const Feature *right. Output: int result. */
static int same_transcript(const Feature *left, const Feature *right) {
  return left->parent && right->parent &&
         left->sequence == right->sequence &&
         !strcmp(left->parent, right->parent) &&
         left->strand == right->strand;
}

/** Input: const FeatureList *features,                                const Feature *target. Output: int result. */
static int transcript_has_exon(const FeatureList *features,
                               const Feature *target) {
  for (size_t i = 0; i < features->count; i++)
    if (features->items[i].kind == FEATURE_EXON &&
        same_transcript(&features->items[i], target))
      return 1;
  return 0;
}

/** Input: const FeatureList *features,                                     const Feature *target, int64_t *start,                                     int64_t *end. Output: int result. */
static int transcript_coding_bounds(const FeatureList *features,
                                    const Feature *target, int64_t *start,
                                    int64_t *end) {
  int found = 0;
  for (size_t i = 0; i < features->count; i++) {
    const Feature *feature = &features->items[i];
    if (feature->kind != FEATURE_CDS || !same_transcript(feature, target))
      continue;
    if (!found || feature->start < *start) *start = feature->start;
    if (!found || feature->end > *end) *end = feature->end;
    found = 1;
  }
  return found;
}

/** Input: Sequence *sequences, int sequence_count,                         const FeatureList *features. Output: none; output parameters are updated or released. */
static void mark_labels(Sequence *sequences, int sequence_count,
                        const FeatureList *features) {
  for (int i = 0; i < sequence_count; i++)
    sequences[i].labels = checked_calloc((size_t)sequences[i].length, 1);
  for (size_t i = 0; i < features->count; i++) {
    const Feature *feature = &features->items[i];
    uint8_t label = explicit_label(feature);
    if (!label) continue;
    for (int64_t position = feature->start; position < feature->end; position++)
      sequences[feature->sequence].labels[position] = label;
  }
  for (size_t i = 0; i < features->count; i++) {
    const Feature *feature = &features->items[i];
    if (feature->kind != FEATURE_EXON &&
        (feature->kind != FEATURE_CDS || transcript_has_exon(features, feature)))
      continue;
    const Feature *next = NULL;
    for (size_t j = 0; j < features->count; j++) {
      const Feature *candidate = &features->items[j];
      if (candidate->kind != feature->kind ||
          !same_transcript(feature, candidate) ||
          candidate->start < feature->end)
        continue;
      if (!next || candidate->start < next->start) next = candidate;
    }
    if (!next) continue;
    uint8_t label = feature->strand == '+'
      ? ANNO_PLUS_INTRON : ANNO_MINUS_INTRON;
    for (int64_t position = feature->end; position < next->start; position++)
      if (!sequences[feature->sequence].labels[position])
        sequences[feature->sequence].labels[position] = label;
  }
  for (size_t i = 0; i < features->count; i++) {
    const Feature *feature = &features->items[i];
    if (feature->kind != FEATURE_EXON) continue;
    int64_t coding_start = 0, coding_end = 0;
    if (!transcript_coding_bounds(features, feature, &coding_start, &coding_end))
      continue;
    uint8_t five_prime = feature->strand == '+'
      ? ANNO_PLUS_5UTR : ANNO_MINUS_5UTR;
    uint8_t three_prime = feature->strand == '+'
      ? ANNO_PLUS_3UTR : ANNO_MINUS_3UTR;
    for (int64_t position = feature->start; position < feature->end; position++) {
      uint8_t label = position < coding_start ? five_prime
        : position >= coding_end ? three_prime : ANNO_INTERGENIC;
      if (label && !sequences[feature->sequence].labels[position])
        sequences[feature->sequence].labels[position] = label;
    }
  }
}

/** Input: const char *sequence, int64_t length, float *complex_values. Output: none; output parameters are updated or released. */
static void cwt_values(const char *sequence, int64_t length, float *complex_values) {
  static const double signal[4][2] = {
    {1.0, 1.0}, {-1.0, 1.0}, {-1.0, -1.0}, {1.0, -1.0}
  };
  for (int64_t position = 0; position < length; position++) {
    for (int scale_index = 0; scale_index < ANNO_CWT_CHANNELS; scale_index++) {
      int radius = (int)cwt_scales[scale_index] * 4;
      double real = 0.0, imaginary = 0.0, normalization = 0.0;
      for (int offset = -radius; offset <= radius; offset++) {
        int64_t source = position + offset;
        if (source < 0 || source >= length) continue;
        int base = base_index(sequence[source]);
        if (base < 0) continue;
        double x = (double)offset / (double)cwt_scales[scale_index];
        double weight = exp(-0.5 * x * x);
        double phase = 6.0 * x;
        real += weight * (signal[base][0] * cos(phase) -
                          signal[base][1] * sin(phase));
        imaginary += weight * (signal[base][0] * sin(phase) +
                               signal[base][1] * cos(phase));
        normalization += weight;
      }
      size_t index = ((size_t)position * ANNO_CWT_CHANNELS +
                      (size_t)scale_index) * 2;
      complex_values[index] = normalization > 0.0
        ? (float)(real / normalization) : 0.0f;
      complex_values[index + 1] =
        normalization > 0.0 ? (float)(imaginary / normalization) : 0.0f;
    }
  }
}

/** Input: const char *sequence, int64_t length, float *power. Output: none; output parameters are updated or released. */
static void cwt_power(const char *sequence, int64_t length, float *power) {
  float *complex_values = checked_calloc((size_t)length * ANNO_CWT_CHANNELS * 2,
                                         sizeof(*complex_values));
  cwt_values(sequence, length, complex_values);
  for (int64_t position = 0; position < length; position++)
    for (int scale = 0; scale < ANNO_CWT_CHANNELS; scale++) {
      size_t complex_index = ((size_t)position * ANNO_CWT_CHANNELS +
                              (size_t)scale) * 2;
      power[(size_t)position * ANNO_CWT_CHANNELS + (size_t)scale] =
        complex_values[complex_index] * complex_values[complex_index] +
        complex_values[complex_index + 1] * complex_values[complex_index + 1];
    }
  free(complex_values);
}

/** Input: int input_channels, int output_channels, int kernel. Output: size_t result. */
static size_t conv_parameters(int input_channels, int output_channels, int kernel) {
  return (size_t)output_channels *
         ((size_t)input_channels * (size_t)kernel + 1);
}

/** Input: const AnnoUnetConfig *config. Output: size_t result. */
static size_t model_head_offset(const AnnoUnetConfig *config) {
  size_t offset = conv_parameters(config->input_channels, config->base_channels,
                                  config->kernel_size);
  for (int level = 1; level < config->depth; level++)
    offset += conv_parameters(config->base_channels, config->base_channels,
                              config->kernel_size);
  for (int level = 1; level < config->depth; level++)
    offset += conv_parameters(config->base_channels, config->base_channels,
                              config->kernel_size);
  return offset;
}

/** Input: const AnnoUnetConfig *config. Output: size_t result. */
static size_t model_parameters(const AnnoUnetConfig *config) {
  return model_head_offset(config) +
         conv_parameters(config->base_channels, config->classes, 1);
}

/** Input: const AnnoModelMetadata *metadata. Output: int result. */
static int valid_metadata(const AnnoModelMetadata *metadata) {
  if (!metadata || !metadata->window_size) return 0;
  for (int channel = 0; channel < 7; channel++)
    if (!isfinite(metadata->feature_mean[channel]) ||
        !isfinite(metadata->feature_std[channel]) ||
        metadata->feature_std[channel] <= 0.0f)
      return 0;
  for (int scale = 0; scale < ANNO_CWT_CHANNELS; scale++)
    if (!metadata->cwt_scales[scale]) return 0;
  for (int class = 0; class < ANNO_CLASSES; class++)
    if (!metadata->class_names[class][0]) return 0;
  return 1;
}

/** Input: const AnnoUnetConfig *config. Output: int result. */
static int valid_config(const AnnoUnetConfig *config) {
  return config && config->input_channels == 4 + ANNO_CWT_CHANNELS &&
         config->classes == ANNO_CLASSES && config->base_channels > 0 &&
         config->base_channels <= 128 && config->depth > 0 &&
         config->depth <= 8 && config->kernel_size > 0 &&
         config->kernel_size <= 15 && (config->kernel_size & 1);
}

/** Input: uint64_t *state. Output: uint32_t result. */
static uint32_t random_u32(uint64_t *state) {
  uint64_t value = *state;
  value ^= value >> 12;
  value ^= value << 25;
  value ^= value >> 27;
  *state = value;
  return (uint32_t)(value * 2685821657736338717ULL >> 32);
}

/** Input: uint64_t *state, float lower, float upper. Output: float result. */
static float random_uniform(uint64_t *state, float lower, float upper) {
  float unit = (float)random_u32(state) / (float)UINT32_MAX;
  return lower + (upper - lower) * unit;
}

/** Input: float *weights, int input_channels,                                    int output_channels, int kernel,                                    uint64_t *seed. Output: none; output parameters are updated or released. */
static void initialize_convolution(float *weights, int input_channels,
                                   int output_channels, int kernel,
                                   uint64_t *seed) {
  size_t per_output = (size_t)input_channels * (size_t)kernel + 1;
  float bound = sqrtf(6.0f / (float)(input_channels * kernel));
  for (int output = 0; output < output_channels; output++) {
    float *filter = weights + (size_t)output * per_output;
    for (int parameter = 0; parameter < input_channels * kernel; parameter++)
      filter[parameter] = random_uniform(seed, -bound, bound);
    filter[input_channels * kernel] = 0.0f;
  }
}

/** Input: const AnnoUnetConfig *config. Output: AnnoModel * result. */
AnnoModel *anno_model_create(const AnnoUnetConfig *config) {
  if (!valid_config(config)) return NULL;
  AnnoModel *model = checked_calloc(1, sizeof(*model));
  model->config = *config;
  model->metadata.window_size = ANNO_CWT_WINDOW;
  for (int channel = 0; channel < config->input_channels; channel++)
    model->metadata.feature_std[channel] = 1.0f;
  memcpy(model->metadata.cwt_scales, cwt_scales, sizeof(cwt_scales));
  for (int class = 0; class < ANNO_CLASSES; class++)
    snprintf(model->metadata.class_names[class],
             sizeof(model->metadata.class_names[class]), "%s", class_names[class]);
  size_t parameters = model_parameters(config);
  model->weights = checked_calloc(parameters, sizeof(*model->weights));
  uint64_t seed = 0x853c49e6748fea9bULL;
  size_t offset = 0;
  initialize_convolution(model->weights + offset, config->input_channels,
                         config->base_channels, config->kernel_size, &seed);
  offset += conv_parameters(config->input_channels, config->base_channels,
                            config->kernel_size);
  for (int level = 1; level < config->depth; level++) {
    initialize_convolution(model->weights + offset, config->base_channels,
                           config->base_channels, config->kernel_size, &seed);
    offset += conv_parameters(config->base_channels, config->base_channels,
                              config->kernel_size);
  }
  for (int level = 1; level < config->depth; level++) {
    initialize_convolution(model->weights + offset, config->base_channels,
                           config->base_channels, config->kernel_size, &seed);
    offset += conv_parameters(config->base_channels, config->base_channels,
                              config->kernel_size);
  }
  initialize_convolution(model->weights + offset, config->base_channels,
                         config->classes, 1, &seed);
  return model;
}

/** Input: AnnoModel *model. Output: none; output parameters are updated or released. */
void anno_model_destroy(AnnoModel *model) {
  if (!model) return;
  free(model->weights);
  free(model);
}

/** Input: const AnnoModel *model. Output: const AnnoModelMetadata * result. */
const AnnoModelMetadata *anno_model_metadata(const AnnoModel *model) {
  return model ? &model->metadata : NULL;
}

/** Input: const float *input, int64_t length, int channels,                            int64_t position, int channel. Output: float result. */
static float feature_value(const float *input, int64_t length, int channels,
                           int64_t position, int channel) {
  if (position < 0 || position >= length) return 0.0f;
  return input[(size_t)position * (size_t)channels + (size_t)channel];
}

/** Input: const float *input, int64_t length, int input_channels,                         int output_channels, int kernel, const float *weights,                         float *output, int relu. Output: none; output parameters are updated or released. */
static void convolution(const float *input, int64_t length, int input_channels,
                        int output_channels, int kernel, const float *weights,
                        float *output, int relu) {
  int radius = kernel / 2;
  size_t per_output = (size_t)input_channels * (size_t)kernel + 1;
  for (int64_t position = 0; position < length; position++)
    for (int output_channel = 0; output_channel < output_channels;
         output_channel++) {
      const float *filter = weights + (size_t)output_channel * per_output;
      float value = filter[input_channels * kernel];
      size_t index = 0;
      for (int offset = -radius; offset <= radius; offset++)
        for (int input_channel = 0; input_channel < input_channels;
             input_channel++)
          value += feature_value(input, length, input_channels,
                                 position + offset, input_channel) * filter[index++];
      output[(size_t)position * (size_t)output_channels +
             (size_t)output_channel] = relu && value < 0.0f ? 0.0f : value;
    }
}

static void softmax(const float *logits, int classes, float *probability);

/** Input: const AnnoUnetConfig *config, int64_t length,                          int64_t *padded_out. Output: int result. */
static int padded_length(const AnnoUnetConfig *config, int64_t length,
                         int64_t *padded_out) {
  int factor = 1 << (config->depth - 1);
  if (length <= 0 || length > INT64_MAX - factor + 1) return EOVERFLOW;
  if ((length + factor - 1) / factor > INT64_MAX / factor) return EOVERFLOW;
  *padded_out = (length + factor - 1) / factor * factor;
  return 0;
}

/** Input: AnnoForwardCache *cache, int depth. Output: none; output parameters are updated or released. */
static void forward_cache_free(AnnoForwardCache *cache, int depth) {
  if (!cache) return;
  for (int level = 0; level < depth; level++) {
    free(cache->enc[level]);
    free(cache->down[level]);
    if (level + 1 < depth) free(cache->dec[level]);
    if (level + 1 < depth) free(cache->up[level]);
  }
  free(cache->enc);
  free(cache->down);
  free(cache->up);
  free(cache->dec);
  free(cache->sizes);
  free(cache->padded_input);
  *cache = (AnnoForwardCache){0};
}

/** Input: const AnnoModel *model, const float *input,                                 int64_t length, AnnoForwardCache *cache. Output: int result. */
static int forward_cache_create(const AnnoModel *model, const float *input,
                                int64_t length, AnnoForwardCache *cache) {
  const AnnoUnetConfig *config = &model->config;
  int status = padded_length(config, length, &cache->padded);
  if (status) return status;
  cache->length = length;
  cache->sizes = checked_calloc((size_t)config->depth, sizeof(*cache->sizes));
  cache->enc = checked_calloc((size_t)config->depth, sizeof(*cache->enc));
  cache->down = checked_calloc((size_t)config->depth, sizeof(*cache->down));
  cache->up = checked_calloc((size_t)config->depth, sizeof(*cache->up));
  cache->dec = checked_calloc((size_t)config->depth, sizeof(*cache->dec));
  cache->padded_input = checked_calloc(
    (size_t)cache->padded * (size_t)config->input_channels, sizeof(float));
  memcpy(cache->padded_input, input,
         (size_t)length * (size_t)config->input_channels * sizeof(float));
  int base = config->base_channels;
  for (int level = 0; level < config->depth; level++)
    cache->sizes[level] = cache->padded >> level;
  size_t offset = 0;
  for (int level = 0; level < config->depth; level++) {
    int64_t size = cache->sizes[level];
    cache->enc[level] = checked_calloc((size_t)size * (size_t)base,
                                       sizeof(**cache->enc));
    if (!level) {
      convolution(cache->padded_input, size, config->input_channels, base,
                  config->kernel_size, model->weights + offset,
                  cache->enc[level], 1);
      offset += conv_parameters(config->input_channels, base,
                                config->kernel_size);
      continue;
    }
    cache->down[level] = checked_calloc((size_t)size * (size_t)base,
                                        sizeof(**cache->down));
    for (int64_t position = 0; position < size; position++)
      for (int channel = 0; channel < base; channel++)
        cache->down[level][(size_t)position * (size_t)base + (size_t)channel] =
          (cache->enc[level - 1][(size_t)(position * 2) * (size_t)base +
                                 (size_t)channel] +
           cache->enc[level - 1][(size_t)(position * 2 + 1) * (size_t)base +
                                 (size_t)channel]) * 0.5f;
    convolution(cache->down[level], size, base, base, config->kernel_size,
                model->weights + offset, cache->enc[level], 1);
    offset += conv_parameters(base, base, config->kernel_size);
  }
  cache->dec[config->depth - 1] = cache->enc[config->depth - 1];
  for (int level = config->depth - 2; level >= 0; level--) {
    int64_t size = cache->sizes[level];
    cache->up[level] = checked_calloc((size_t)size * (size_t)base,
                                      sizeof(**cache->up));
    cache->dec[level] = checked_calloc((size_t)size * (size_t)base,
                                       sizeof(**cache->dec));
    for (int64_t position = 0; position < size; position++)
      for (int channel = 0; channel < base; channel++)
        cache->up[level][(size_t)position * (size_t)base + (size_t)channel] =
          cache->enc[level][(size_t)position * (size_t)base + (size_t)channel] +
          cache->dec[level + 1][(size_t)(position / 2) * (size_t)base +
                                (size_t)channel];
    convolution(cache->up[level], size, base, base, config->kernel_size,
                model->weights + offset, cache->dec[level], 1);
    offset += conv_parameters(base, base, config->kernel_size);
  }
  return 0;
}

/** Input: const float *input, const float *output,                                  const float *gradient_output, int64_t length,                                  int input_channels, int output_channels,                                  int kernel, float *weights, float *gradient_input,                                  float *gradient_weights. Output: none; output parameters are updated or released. */
static void convolution_backward(const float *input, const float *output,
                                 const float *gradient_output, int64_t length,
                                 int input_channels, int output_channels,
                                 int kernel, float *weights, float *gradient_input,
                                 float *gradient_weights) {
  int radius = kernel / 2;
  size_t per_output = (size_t)input_channels * (size_t)kernel + 1;
  for (int64_t position = 0; position < length; position++)
    for (int output_channel = 0; output_channel < output_channels;
         output_channel++) {
      size_t output_index = (size_t)position * (size_t)output_channels +
                            (size_t)output_channel;
      float delta = output[output_index] > 0.0f
        ? gradient_output[output_index] : 0.0f;
      if (delta == 0.0f) continue;
      float *filter = weights + (size_t)output_channel * per_output;
      size_t index = 0;
      for (int offset = -radius; offset <= radius; offset++)
        for (int input_channel = 0; input_channel < input_channels;
             input_channel++) {
          int64_t source = position + offset;
          float old_weight = filter[index];
          if (gradient_input && source >= 0 && source < length)
            gradient_input[(size_t)source * (size_t)input_channels +
                           (size_t)input_channel] += delta * old_weight;
          if (source >= 0 && source < length)
            gradient_weights[(size_t)output_channel * per_output + index] +=
              delta * input[(size_t)source * (size_t)input_channels +
                            (size_t)input_channel];
          index++;
        }
      gradient_weights[(size_t)output_channel * per_output +
                       (size_t)input_channels * (size_t)kernel] += delta;
    }
}

/** Input: float *weights,                                        const float *gradient_weights,                                        int input_channels, int output_channels,                                        int kernel, float learning_rate. Output: none; output parameters are updated or released. */
static void apply_convolution_gradient(float *weights,
                                       const float *gradient_weights,
                                       int input_channels, int output_channels,
                                       int kernel, float learning_rate) {
  size_t parameters = conv_parameters(input_channels, output_channels, kernel);
  for (size_t parameter = 0; parameter < parameters; parameter++)
    weights[parameter] -= learning_rate * gradient_weights[parameter];
}

/** Input: const AnnoUnetConfig *config, int level. Output: size_t result. */
static size_t encoder_offset(const AnnoUnetConfig *config, int level) {
  size_t offset = 0;
  offset += conv_parameters(config->input_channels, config->base_channels,
                            config->kernel_size);
  for (int current = 1; current < level; current++)
    offset += conv_parameters(config->base_channels, config->base_channels,
                              config->kernel_size);
  return offset;
}

/** Input: AnnoModel *model,                                const AnnoForwardCache *cache,                                const int64_t *order, int64_t batch_start,                                int64_t batch_end, const uint8_t *label,                                const float *class_weight, float batch_weight,                                float learning_rate. Output: none; output parameters are updated or released. */
static void train_cached_batch(AnnoModel *model,
                               const AnnoForwardCache *cache,
                               const int64_t *order, int64_t batch_start,
                               int64_t batch_end, const uint8_t *label,
                               const float *class_weight, float batch_weight,
                               float learning_rate) {
  const AnnoUnetConfig *config = &model->config;
  int depth = config->depth;
  int base = config->base_channels;
  size_t per_class = (size_t)base + 1;
  size_t head_offset = model_head_offset(config);
  float **gradient_dec = checked_calloc((size_t)depth, sizeof(*gradient_dec));
  float **gradient_enc = checked_calloc((size_t)depth, sizeof(*gradient_enc));
  float *head_gradient = checked_calloc(
    (size_t)config->classes * per_class, sizeof(*head_gradient));
  for (int level = 0; level < depth; level++) {
    size_t cells = (size_t)cache->sizes[level] * (size_t)base;
    gradient_dec[level] = checked_calloc(cells, sizeof(**gradient_dec));
    gradient_enc[level] = checked_calloc(cells, sizeof(**gradient_enc));
  }
  for (int64_t batch_position = batch_start; batch_position < batch_end;
       batch_position++) {
    int64_t position = order[batch_position];
    float logits[ANNO_CLASSES], probability[ANNO_CLASSES];
    const float *features = cache->dec[0] +
      (size_t)position * (size_t)base;
    for (int class = 0; class < config->classes; class++) {
      const float *filter = model->weights + head_offset +
        (size_t)class * per_class;
      float value = filter[base];
      for (int channel = 0; channel < base; channel++)
        value += features[channel] * filter[channel];
      logits[class] = value;
    }
    softmax(logits, config->classes, probability);
    float sample_scale = class_weight[label[position]] / batch_weight;
    for (int class = 0; class < config->classes; class++) {
      float error = (probability[class] -
                     (label[position] == class ? 1.0f : 0.0f)) * sample_scale;
      const float *filter = model->weights + head_offset +
        (size_t)class * per_class;
      float *head_delta = head_gradient + (size_t)class * per_class;
      for (int channel = 0; channel < base; channel++) {
        head_delta[channel] += error * features[channel];
        gradient_dec[0][(size_t)position * (size_t)base + (size_t)channel] +=
          error * filter[channel];
      }
      head_delta[base] += error;
    }
  }
  size_t encoder_end = model_head_offset(config) -
    (size_t)(depth - 1) * conv_parameters(base, base, config->kernel_size);
  size_t decoder_offset = encoder_end;
  for (int level = depth - 2; level >= 0; level--) {
    size_t cells = (size_t)cache->sizes[level] * (size_t)base;
    float *gradient_up = checked_calloc(cells, sizeof(*gradient_up));
    size_t parameters = conv_parameters(base, base, config->kernel_size);
    float *gradient_weights = checked_calloc(parameters, sizeof(*gradient_weights));
    convolution_backward(cache->up[level], cache->dec[level],
                         gradient_dec[level], cache->sizes[level], base, base,
                         config->kernel_size,
                         model->weights + decoder_offset, gradient_up,
                         gradient_weights);
    apply_convolution_gradient(model->weights + decoder_offset, gradient_weights,
                               base, base, config->kernel_size, learning_rate);
    free(gradient_weights);
    for (int64_t position_at_level = 0;
         position_at_level < cache->sizes[level]; position_at_level++)
      for (int channel = 0; channel < base; channel++) {
        size_t index = (size_t)position_at_level * (size_t)base +
                      (size_t)channel;
        gradient_enc[level][index] += gradient_up[index];
        gradient_dec[level + 1][(size_t)(position_at_level / 2) *
                                (size_t)base + (size_t)channel] +=
          gradient_up[index];
      }
    free(gradient_up);
    decoder_offset += conv_parameters(base, base, config->kernel_size);
  }
  for (int channel = 0; channel < base; channel++)
    for (int64_t position_at_level = 0;
         position_at_level < cache->sizes[depth - 1]; position_at_level++)
      gradient_enc[depth - 1][(size_t)position_at_level * (size_t)base +
                              (size_t)channel] +=
        gradient_dec[depth - 1][(size_t)position_at_level * (size_t)base +
                                (size_t)channel];
  for (int level = depth - 1; level >= 0; level--) {
    float *gradient_down = NULL;
    if (level > 0)
      gradient_down = checked_calloc(
        (size_t)cache->sizes[level] * (size_t)base, sizeof(*gradient_down));
    int input_channels = level ? base : config->input_channels;
    size_t parameters = conv_parameters(input_channels, base,
                                        config->kernel_size);
    float *gradient_weights = checked_calloc(parameters, sizeof(*gradient_weights));
    convolution_backward(
      level ? cache->down[level] : cache->padded_input, cache->enc[level],
      gradient_enc[level], cache->sizes[level],
      input_channels, base, config->kernel_size,
      model->weights + encoder_offset(config, level), gradient_down,
      gradient_weights);
    apply_convolution_gradient(model->weights + encoder_offset(config, level),
                               gradient_weights, input_channels, base,
                               config->kernel_size, learning_rate);
    free(gradient_weights);
    if (level > 0) {
      for (int64_t position_at_level = 0;
           position_at_level < cache->sizes[level]; position_at_level++)
        for (int channel = 0; channel < base; channel++) {
          float gradient = gradient_down[(size_t)position_at_level *
                                         (size_t)base + (size_t)channel] * 0.5f;
          gradient_enc[level - 1][(size_t)(position_at_level * 2) *
                                  (size_t)base + (size_t)channel] += gradient;
          gradient_enc[level - 1][(size_t)(position_at_level * 2 + 1) *
                                  (size_t)base + (size_t)channel] += gradient;
        }
      free(gradient_down);
    }
  }
  for (int level = 0; level < depth; level++) {
    free(gradient_dec[level]);
    free(gradient_enc[level]);
  }
  for (int class = 0; class < config->classes; class++) {
    float *filter = model->weights + head_offset + (size_t)class * per_class;
    float *head_delta = head_gradient + (size_t)class * per_class;
    for (size_t parameter = 0; parameter < per_class; parameter++)
      filter[parameter] -= learning_rate * head_delta[parameter];
  }
  free(head_gradient);
  free(gradient_dec);
  free(gradient_enc);
}

/** Input: const AnnoModel *model, const float *input,                           int64_t length, float **features_out. Output: int result. */
static int model_features(const AnnoModel *model, const float *input,
                          int64_t length, float **features_out) {
  AnnoForwardCache cache = {0};
  int status = forward_cache_create(model, input, length, &cache);
  if (status) return status;
  size_t values = (size_t)length * (size_t)model->config.base_channels;
  float *features = checked_calloc(values, sizeof(*features));
  memcpy(features, cache.dec[0], values * sizeof(*features));
  forward_cache_free(&cache, model->config.depth);
  *features_out = features;
  return 0;
}

/** Input: const float *logits, int classes, float *probability. Output: none; output parameters are updated or released. */
static void softmax(const float *logits, int classes, float *probability) {
  float maximum = -INFINITY;
  for (int class = 0; class < classes; class++)
    if (logits[class] > maximum) maximum = logits[class];
  double total = 0.0;
  for (int class = 0; class < classes; class++) {
    probability[class] = expf(logits[class] - maximum);
    total += probability[class];
  }
  if (!(total > 0.0) || !isfinite(total)) {
    for (int class = 0; class < classes; class++)
      probability[class] = 1.0f / (float)classes;
    return;
  }
  for (int class = 0; class < classes; class++)
    probability[class] /= (float)total;
}

/** Input: const AnnoModel *model, const float *input,                        int64_t length, float *probability. Output: int result. */
int anno_model_predict(const AnnoModel *model, const float *input,
                       int64_t length, float *probability) {
  if (!model || !input || !probability || length <= 0 ||
      !valid_config(&model->config)) return EINVAL;
  float *features = NULL;
  int status = model_features(model, input, length, &features);
  if (status) return status;
  size_t offset = model_head_offset(&model->config);
  size_t per_class = (size_t)model->config.base_channels + 1;
  float logits[ANNO_CLASSES], values[ANNO_CLASSES];
  for (int64_t position = 0; position < length; position++) {
    for (int class = 0; class < model->config.classes; class++) {
      const float *filter = model->weights + offset + (size_t)class * per_class;
      float value = filter[model->config.base_channels];
      for (int channel = 0; channel < model->config.base_channels; channel++)
        value += features[(size_t)position * (size_t)model->config.base_channels +
                          (size_t)channel] * filter[channel];
      logits[class] = value;
    }
    softmax(logits, model->config.classes, values);
    memcpy(probability + (size_t)position * ANNO_CLASSES, values,
           sizeof(values));
  }
  free(features);
  return 0;
}

/** Input: AnnoModel *model, const float *input, const uint8_t *label,                      int64_t length, int iterations, int batch_size,                      float learning_rate, uint64_t seed. Output: int result. */
int anno_model_train(AnnoModel *model, const float *input, const uint8_t *label,
                     int64_t length, int iterations, int batch_size,
                     float learning_rate, uint64_t seed) {
  if (!model || !input || !label || length <= 0 || iterations < 1 ||
      batch_size < 1 || !(learning_rate > 0.0f) || !isfinite(learning_rate))
    return EINVAL;
  int64_t *order = checked_calloc((size_t)length, sizeof(*order));
  int64_t class_counts[ANNO_CLASSES] = {0};
  for (int64_t position = 0; position < length; position++) {
    if (label[position] >= ANNO_CLASSES) {
      free(order);
      return EINVAL;
    }
    class_counts[label[position]]++;
  }
  float class_weight[ANNO_CLASSES];
  for (int class = 0; class < ANNO_CLASSES; class++)
    class_weight[class] = class_counts[class]
      ? (float)length / ((float)ANNO_CLASSES * (float)class_counts[class])
      : 0.0f;
  for (int iteration = 0; iteration < iterations; iteration++) {
    AnnoForwardCache cache = {0};
    int status = forward_cache_create(model, input, length, &cache);
    if (status) {
      free(order);
      return status;
    }
    for (int64_t position = 0; position < length; position++)
      order[position] = position;
    for (int64_t position = length - 1; position > 0; position--) {
      int64_t swap = (int64_t)(random_u32(&seed) % (uint64_t)(position + 1));
      int64_t temporary = order[position];
      order[position] = order[swap];
      order[swap] = temporary;
    }
    for (int64_t batch_start = 0; batch_start < length;
         batch_start += batch_size) {
      int64_t batch_end = batch_start + batch_size;
      if (batch_end > length) batch_end = length;
      float batch_weight = 0.0f;
      for (int64_t batch_position = batch_start; batch_position < batch_end;
           batch_position++)
        batch_weight += class_weight[label[order[batch_position]]];
      if (!(batch_weight > 0.0f)) continue;
      train_cached_batch(model, &cache, order, batch_start, batch_end, label,
                         class_weight, batch_weight, learning_rate);
    }
    forward_cache_free(&cache, model->config.depth);
  }
  free(order);
  return 0;
}

/** Input: uint32_t checksum, const void *data, size_t size. Output: uint32_t result. */
static uint32_t checksum_update(uint32_t checksum, const void *data, size_t size) {
  const unsigned char *bytes = data;
  for (size_t i = 0; i < size; i++)
    checksum = (checksum ^ bytes[i]) * 16777619u;
  return checksum;
}

/** Input: const AnnoModel *model. Output: uint32_t result. */
static uint32_t model_checksum(const AnnoModel *model) {
  uint32_t checksum = 2166136261u;
  checksum = checksum_update(checksum, &model->metadata, sizeof(model->metadata));
  checksum = checksum_update(checksum, model->weights,
                             model_parameters(&model->config) * sizeof(float));
  return checksum;
}

/** Input: const AnnoModel *model, const char *path. Output: int result. */
int anno_model_save(const AnnoModel *model, const char *path) {
  if (!model || !path || !valid_config(&model->config) ||
      !valid_metadata(&model->metadata))
    return EINVAL;
  size_t parameters = model_parameters(&model->config);
  if (parameters > UINT32_MAX) return EOVERFLOW;
  FILE *file = fopen(path, "wb");
  if (!file) return errno;
  uint32_t header[10] = {
    model_version, model_endian, (uint32_t)model->config.input_channels,
    (uint32_t)model->config.classes, (uint32_t)model->config.base_channels,
    (uint32_t)model->config.depth, (uint32_t)model->config.kernel_size,
    (uint32_t)parameters, model->metadata.window_size, model_checksum(model)
  };
  int ok = fwrite(model_magic, 1, sizeof(model_magic), file) == sizeof(model_magic) &&
           fwrite(header, sizeof(*header), 10, file) == 10 &&
           fwrite(&model->metadata, sizeof(model->metadata), 1, file) == 1 &&
           fwrite(model->weights, sizeof(float), parameters, file) == parameters;
  int close_status = fclose(file);
  return ok && close_status == 0 ? 0 : EIO;
}

/** Input: AnnoModel **model_out, const char *path. Output: int result. */
int anno_model_load(AnnoModel **model_out, const char *path) {
  if (!model_out || !path) return EINVAL;
  *model_out = NULL;
  FILE *file = fopen(path, "rb");
  if (!file) return errno;
  char magic[sizeof(model_magic)];
  uint32_t header[10] = {0};
  int ok = fread(magic, 1, sizeof(magic), file) == sizeof(magic) &&
           fread(header, sizeof(*header), 10, file) == 10 &&
           !memcmp(magic, model_magic, sizeof(magic));
  if (!ok || header[0] != model_version || header[1] != model_endian) {
    fclose(file);
    return EINVAL;
  }
  AnnoUnetConfig config = {
    (int)header[2], (int)header[3], (int)header[4], (int)header[5], (int)header[6]
  };
  if (!valid_config(&config) || model_parameters(&config) != header[7]) {
    fclose(file);
    return EINVAL;
  }
  AnnoModel *model = anno_model_create(&config);
  if (!model) {
    fclose(file);
    return ENOMEM;
  }
  ok = fread(&model->metadata, sizeof(model->metadata), 1, file) == 1 &&
       fread(model->weights, sizeof(float), header[7], file) == header[7] &&
       fgetc(file) == EOF;
  int close_status = fclose(file);
  if (!ok || close_status || model->metadata.window_size != header[8] ||
      !valid_metadata(&model->metadata) ||
      model_checksum(model) != header[9]) {
    anno_model_destroy(model);
    return EINVAL;
  }
  for (size_t i = 0; i < model_parameters(&config); i++)
    if (!isfinite(model->weights[i])) {
      anno_model_destroy(model);
      return ERANGE;
    }
  *model_out = model;
  return 0;
}

/** Input: uint8_t state. Output: int result. */
int anno_class_for_state(uint8_t state) {
  static const int classes[ANNO_DECODED_STATES] = {0, 2, 0, 1, 3, 5, 0, 4, 6};
  return state < ANNO_DECODED_STATES ? classes[state] : -1;
}

/** Input: uint8_t previous, uint8_t current. Output: int result. */
int anno_grammar_transition(uint8_t previous, uint8_t current) {
  if (previous >= ANNO_DECODED_STATES || current >= ANNO_DECODED_STATES)
    return 0;
  if (previous == current) return 1;
  if (previous == ANNO_STATE_INTERGENIC)
    return current == ANNO_STATE_PLUS_5UTR || current == ANNO_STATE_PLUS_EXON ||
           current == ANNO_STATE_MINUS_5UTR || current == ANNO_STATE_MINUS_EXON;
  if (current == ANNO_STATE_INTERGENIC)
    return previous == ANNO_STATE_PLUS_3UTR ||
           previous == ANNO_STATE_MINUS_3UTR ||
           previous == ANNO_STATE_PLUS_EXON ||
           previous == ANNO_STATE_MINUS_EXON;
  if (previous == ANNO_STATE_PLUS_5UTR) return current == ANNO_STATE_PLUS_EXON;
  if (previous == ANNO_STATE_PLUS_EXON)
    return current == ANNO_STATE_PLUS_INTRON || current == ANNO_STATE_PLUS_3UTR;
  if (previous == ANNO_STATE_PLUS_INTRON) return current == ANNO_STATE_PLUS_EXON;
  if (previous == ANNO_STATE_MINUS_5UTR) return current == ANNO_STATE_MINUS_EXON;
  if (previous == ANNO_STATE_MINUS_EXON)
    return current == ANNO_STATE_MINUS_INTRON || current == ANNO_STATE_MINUS_3UTR;
  if (previous == ANNO_STATE_MINUS_INTRON) return current == ANNO_STATE_MINUS_EXON;
  return 0;
}

/** Input: uint8_t state. Output: int result. */
static int state_is_intron(uint8_t state) {
  return state == ANNO_STATE_PLUS_INTRON || state == ANNO_STATE_MINUS_INTRON;
}

/* Input: one DP block and its state buffers. Output: updated DP rows and backpointers. */
/** Input: const float *probability, int64_t block,                          int64_t block_size, int64_t length, int max_run,                          int min_exon, int min_intron, size_t row_cells,                          size_t row_bytes, float *checkpoints,                          float **previous_io, float **current_io,                          uint16_t *backs, int clamp_initial. Output: none; output parameters are updated or released. */
static void decode_block(const float *probability, int64_t block,
                         int64_t block_size, int64_t length, int max_run,
                         int min_exon, int min_intron, size_t row_cells,
                         size_t row_bytes, float *checkpoints,
                         float **previous_io, float **current_io,
                         uint16_t *backs, int clamp_initial) {
  int64_t start = block * block_size;
  int64_t end = start + block_size;
  float *previous = *previous_io;
  float *current = *current_io;
  if (end > length) end = length;
  if (!block) {
    for (size_t cell = 0; cell < row_cells; cell++) previous[cell] = -INFINITY;
    for (int state = 0; state < ANNO_DECODED_STATES; state++) {
      if (state != ANNO_STATE_INTERGENIC &&
          state != ANNO_STATE_PLUS_5UTR &&
          state != ANNO_STATE_PLUS_EXON &&
          state != ANNO_STATE_MINUS_5UTR &&
          state != ANNO_STATE_MINUS_EXON)
        continue;
      int class = anno_class_for_state((uint8_t)state);
      float value = probability[class];
      if (clamp_initial && !(value > FLT_MIN && isfinite(value)))
        value = FLT_MIN;
      previous[(size_t)state * (size_t)max_run] = logf(value);
    }
    start++;
  } else {
    memcpy(previous, checkpoints + (size_t)(block - 1) * row_cells, row_bytes);
  }
  for (int64_t position = start; position < end; position++) {
    for (size_t cell = 0; cell < row_cells; cell++) {
      current[cell] = -INFINITY;
      backs[(size_t)(position - block * block_size) * row_cells + cell] =
        UINT16_MAX;
    }
    for (int state = 0; state < ANNO_DECODED_STATES; state++) {
      int class = anno_class_for_state((uint8_t)state);
      float value = probability[(size_t)position * ANNO_CLASSES +
                                (size_t)class];
      float node = logf(value > FLT_MIN && isfinite(value) ? value : FLT_MIN);
      for (int previous_state = 0;
           previous_state < ANNO_DECODED_STATES; previous_state++)
        for (int run = 0; run < max_run; run++) {
          size_t previous_cell = (size_t)previous_state * (size_t)max_run +
                                 (size_t)run;
          float previous_score = previous[previous_cell];
          if (!isfinite(previous_score)) continue;
          int next_run = -1;
          if (state == previous_state)
            next_run = run + 1 < max_run ? run + 1 : max_run - 1;
          else if (anno_grammar_transition((uint8_t)previous_state,
                                           (uint8_t)state) &&
                   (previous_state == ANNO_STATE_INTERGENIC ||
                    run + 1 >=
                      (state_is_intron((uint8_t)previous_state)
                         ? min_intron : min_exon)))
            next_run = 0;
          if (next_run < 0) continue;
          size_t cell = (size_t)state * (size_t)max_run +
                        (size_t)next_run;
          if (previous_score + node > current[cell]) {
            current[cell] = previous_score + node;
            backs[(size_t)(position - block * block_size) * row_cells + cell] =
              (uint16_t)(previous_state * max_run + run);
          }
        }
    }
    float *swap = previous;
    previous = current;
    current = swap;
  }
  *previous_io = previous;
  *current_io = current;
  memcpy(checkpoints + (size_t)block * row_cells, previous, row_bytes);
}

/** Input: const float *probability, int64_t length,                         int min_exon, int min_intron, uint8_t *path. Output: int result. */
int anno_grammar_decode(const float *probability, int64_t length,
                        int min_exon, int min_intron, uint8_t *path) {
  if (!probability || !path || length < 1 || min_exon < 1 || min_intron < 1)
    return 0;
  int max_run = min_exon > min_intron ? min_exon : min_intron;
  size_t row_cells = (size_t)ANNO_DECODED_STATES * (size_t)max_run;
  if (max_run > UINT16_MAX / ANNO_DECODED_STATES ||
      row_cells > SIZE_MAX / sizeof(float))
    return 0;
  const int64_t block_size = 8192;
  int64_t block_count = (length + block_size - 1) / block_size;
  if ((size_t)block_count > SIZE_MAX / row_cells ||
      (size_t)block_count * row_cells > SIZE_MAX / sizeof(float))
    return 0;
  float *checkpoints = malloc((size_t)block_count * row_cells *
                              sizeof(*checkpoints));
  float *previous = malloc(row_cells * sizeof(*previous));
  float *current = malloc(row_cells * sizeof(*current));
  if (block_size > SIZE_MAX / row_cells ||
      (size_t)block_size * row_cells > SIZE_MAX / sizeof(uint16_t)) {
    free(checkpoints);
    free(previous);
    free(current);
    return 0;
  }
  size_t block_cells = (size_t)block_size * row_cells;
  uint16_t *backs = malloc(block_cells * sizeof(*backs));
  if (!checkpoints || !previous || !current || !backs) {
    free(checkpoints);
    free(previous);
    free(current);
    free(backs);
    return 0;
  }
  size_t row_bytes = row_cells * sizeof(*previous);
  for (int64_t block = 0; block < block_count; block++)
    decode_block(probability, block, block_size, length, max_run, min_exon,
                 min_intron, row_cells, row_bytes, checkpoints, &previous,
                 &current, backs, 1);
  int best_state = ANNO_STATE_INTERGENIC, best_run = 0;
  float best = -INFINITY;
  for (int state = 0; state < ANNO_DECODED_STATES; state++)
    for (int run = 0; run < max_run; run++) {
      if (state != ANNO_STATE_INTERGENIC &&
          run + 1 < (state_is_intron((uint8_t)state) ? min_intron : min_exon))
        continue;
      if (state == ANNO_STATE_PLUS_5UTR ||
          state == ANNO_STATE_PLUS_INTRON ||
          state == ANNO_STATE_MINUS_5UTR ||
          state == ANNO_STATE_MINUS_INTRON)
        continue;
      float score = checkpoints[(size_t)(block_count - 1) * row_cells +
                               (size_t)state * (size_t)max_run + (size_t)run];
      if (score > best) {
        best = score;
        best_state = state;
        best_run = run;
      }
    }
  if (!isfinite(best)) {
    free(checkpoints);
    free(previous);
    free(current);
    free(backs);
    return 0;
  }
  for (int64_t block = block_count - 1; block >= 0; block--) {
    int64_t end = block * block_size + block_size;
    if (end > length) end = length;
    decode_block(probability, block, block_size, length, max_run, min_exon,
                 min_intron, row_cells, row_bytes, checkpoints, &previous,
                 &current, backs, 0);
    for (int64_t position = end - 1; position >= block * block_size;
         position--) {
      path[position] = (uint8_t)best_state;
      if (!position) break;
      uint16_t back = backs[(size_t)(position - block * block_size) * row_cells +
                            (size_t)best_state * (size_t)max_run +
                            (size_t)best_run];
      if (back == UINT16_MAX) {
        free(checkpoints);
        free(previous);
        free(current);
        free(backs);
        return 0;
      }
      best_state = back / (uint16_t)max_run;
      best_run = back % (uint16_t)max_run;
    }
  }
  free(checkpoints);
  free(previous);
  free(current);
  free(backs);
  return 1;
}

/** Input: uint8_t state, char strand. Output: int result. */
static int path_is_transcript(uint8_t state, char strand) {
  return strand == '+'
    ? state >= ANNO_STATE_PLUS_5UTR && state <= ANNO_STATE_PLUS_3UTR
    : state >= ANNO_STATE_MINUS_5UTR && state <= ANNO_STATE_MINUS_3UTR;
}

/** Input: uint8_t state. Output: int result. */
static int path_is_intron(uint8_t state) {
  return state == ANNO_STATE_PLUS_INTRON || state == ANNO_STATE_MINUS_INTRON;
}

/** Input: Segment **segments, int *count,                            int64_t start, int64_t end. Output: none; output parameters are updated or released. */
static void append_segment(Segment **segments, int *count,
                           int64_t start, int64_t end) {
  *segments = checked_realloc(*segments, (size_t)*count + 1, sizeof(**segments));
  (*segments)[(*count)++] = (Segment){start, end};
}

/** Input: const char *sequence, int64_t sequence_length,                     const Segment *exons, int exon_count, char strand,                     Orf *orf. Output: int result. */
static int find_orf(const char *sequence, int64_t sequence_length,
                    const Segment *exons, int exon_count, char strand,
                    Orf *orf) {
  (void)sequence_length;
  int64_t transcript_length = 0;
  for (int i = 0; i < exon_count; i++)
    transcript_length += exons[i].end - exons[i].start;
  if (transcript_length < 6) return 0;
  char *transcript = checked_calloc((size_t)transcript_length + 1, 1);
  int64_t at = 0;
  for (int order = 0; order < exon_count; order++) {
    int i = strand == '+' ? order : exon_count - 1 - order;
    if (strand == '+')
      for (int64_t position = exons[i].start; position < exons[i].end; position++)
        transcript[at++] = sequence[position];
    else
      for (int64_t position = exons[i].end; position > exons[i].start;)
        transcript[at++] = complement(sequence[--position]);
  }
  int found = 0;
  int64_t best_start = 0, best_end = 0;
  for (int64_t start = 0; start + 2 < transcript_length; start++)
    if (is_atg(transcript + start))
      for (int64_t end = start + 3; end + 2 < transcript_length; end += 3)
        if (is_stop(transcript + end) &&
            (!found || end + 3 - start > best_end - best_start)) {
          found = 1;
          best_start = start;
          best_end = end + 3;
          break;
        }
  free(transcript);
  if (!found) return 0;
  int64_t offset = 0;
  for (int order = 0; order < exon_count; order++) {
    int i = strand == '+' ? order : exon_count - 1 - order;
    int64_t length = exons[i].end - exons[i].start;
    if (best_start < offset + length && best_end > offset) {
      int64_t local_start = best_start > offset ? best_start - offset : 0;
      int64_t local_end = best_end < offset + length ? best_end - offset : length;
      orf->start = strand == '+' ? exons[i].start + local_start
                                  : exons[i].end - local_end;
      orf->end = strand == '+' ? exons[i].start + local_end
                                : exons[i].end - local_start;
      break;
    }
    offset += length;
  }
  orf->transcript_start = best_start;
  orf->transcript_end = best_end;
  orf->length = best_end - best_start;
  orf->strand = strand;
  return 1;
}

/** Input: const Segment *exons, int exon_count,                                 int64_t from, int64_t to, char strand,                                 Segment **mapped, int *mapped_count. Output: int result. */
static int map_transcript_range(const Segment *exons, int exon_count,
                                int64_t from, int64_t to, char strand,
                                Segment **mapped, int *mapped_count) {
  int64_t offset = 0;
  for (int order = 0; order < exon_count; order++) {
    int i = strand == '+' ? order : exon_count - 1 - order;
    int64_t length = exons[i].end - exons[i].start;
    int64_t local_from = from > offset ? from - offset : 0;
    int64_t local_to = to < offset + length ? to - offset : length;
    if (local_from < local_to) {
      int64_t start = strand == '+' ? exons[i].start + local_from
                                    : exons[i].end - local_to;
      int64_t end = strand == '+' ? exons[i].start + local_to
                                  : exons[i].end - local_from;
      append_segment(mapped, mapped_count, start, end);
    }
    offset += length;
  }
  return *mapped_count > 0;
}

/** Input: FILE *output, const char *seqid, const char *type,                     int64_t start, int64_t end, char strand,                     int phase, const char *attributes. Output: int result. */
static int emit_gff(FILE *output, const char *seqid, const char *type,
                    int64_t start, int64_t end, char strand,
                    int phase, const char *attributes) {
  if (start < 0 || end <= start || (strand != '+' && strand != '-')) return 0;
  if (fprintf(output, "%s\tanno\t%s\t%lld\t%lld\t.\t%c\t",
              seqid, type, (long long)start + 1, (long long)end, strand) < 0)
    return 0;
  if (phase < 0) {
    if (fputs(".\t", output) == EOF) return 0;
  } else if (fprintf(output, "%d\t", phase) < 0) return 0;
  return fprintf(output, "%s\n", attributes) >= 0;
}

static float *normalized_features(const Sequence *sequence,
                                  const AnnoModelMetadata *metadata);

/** Input: FILE *output, const Sequence *sequence,                                const uint8_t *path, int64_t block_start,                                int64_t block_end, char strand,                                unsigned long gene_number. Output: int result. */
static int write_decoded_block(FILE *output, const Sequence *sequence,
                               const uint8_t *path, int64_t block_start,
                               int64_t block_end, char strand,
                               unsigned long gene_number) {
  Segment *exons = NULL;
  int exon_count = 0;
  for (int64_t position = block_start; position < block_end;) {
    while (position < block_end && path_is_intron(path[position])) position++;
    if (position == block_end) break;
    int64_t start = position;
    while (position < block_end && !path_is_intron(path[position])) position++;
    append_segment(&exons, &exon_count, start, position);
  }
  Orf orf;
  if (!exon_count || !find_orf(sequence->sequence, sequence->length, exons,
                                exon_count, strand, &orf)) {
    free(exons);
    return 0;
  }
  Segment *cds = NULL;
  int cds_count = 0;
  map_transcript_range(exons, exon_count, orf.transcript_start,
                       orf.transcript_end, strand, &cds, &cds_count);
  char gene_id[64], transcript_id[64], attributes[256];
  snprintf(gene_id, sizeof(gene_id), "gene%lu", gene_number);
  snprintf(transcript_id, sizeof(transcript_id), "mRNA%lu", gene_number);
  int64_t gene_start = exons[0].start, gene_end = exons[0].end;
  for (int i = 1; i < exon_count; i++) {
    if (exons[i].start < gene_start) gene_start = exons[i].start;
    if (exons[i].end > gene_end) gene_end = exons[i].end;
  }
  snprintf(attributes, sizeof(attributes), "ID=%s", gene_id);
  int ok = emit_gff(output, sequence->name, "gene", gene_start, gene_end,
                    strand, -1, attributes);
  snprintf(attributes, sizeof(attributes), "ID=%s;Parent=%s", transcript_id,
           gene_id);
  ok = ok && emit_gff(output, sequence->name, "mRNA", gene_start, gene_end,
                      strand, -1, attributes);
  for (int i = 0; ok && i < exon_count; i++) {
    snprintf(attributes, sizeof(attributes), "ID=%s.exon%d;Parent=%s",
             transcript_id, i + 1, transcript_id);
    ok = emit_gff(output, sequence->name, "exon", exons[i].start, exons[i].end,
                  strand, -1, attributes);
  }
  int64_t transcript_length = 0;
  for (int i = 0; i < exon_count; i++)
    transcript_length += exons[i].end - exons[i].start;
  int utr5 = 0, utr3 = 0;
  Segment *utr = NULL;
  int utr_count = 0;
  if (orf.transcript_start > 0 &&
      map_transcript_range(exons, exon_count, 0, orf.transcript_start, strand,
                           &utr, &utr_count))
    for (int i = 0; ok && i < utr_count; i++) {
      snprintf(attributes, sizeof(attributes), "ID=%s.five_prime_UTR%d;Parent=%s",
               transcript_id, ++utr5, transcript_id);
      ok = emit_gff(output, sequence->name, "five_prime_UTR", utr[i].start,
                    utr[i].end, strand, -1, attributes);
    }
  free(utr);
  for (int i = 0; ok && i < cds_count; i++) {
    int64_t coding_offset = 0;
    for (int j = 0; j < i; j++) coding_offset += cds[j].end - cds[j].start;
    int phase = (3 - (int)(coding_offset % 3)) % 3;
    snprintf(attributes, sizeof(attributes), "ID=%s.cds%d;Parent=%s",
             transcript_id, i + 1, transcript_id);
    ok = emit_gff(output, sequence->name, "CDS", cds[i].start, cds[i].end,
                  strand, phase, attributes);
  }
  utr = NULL;
  utr_count = 0;
  if (orf.transcript_end < transcript_length &&
      map_transcript_range(exons, exon_count, orf.transcript_end,
                           transcript_length, strand, &utr, &utr_count))
    for (int i = 0; ok && i < utr_count; i++) {
      snprintf(attributes, sizeof(attributes), "ID=%s.three_prime_UTR%d;Parent=%s",
               transcript_id, ++utr3, transcript_id);
      ok = emit_gff(output, sequence->name, "three_prime_UTR", utr[i].start,
                    utr[i].end, strand, -1, attributes);
    }
  free(utr);
  free(cds);
  free(exons);
  return ok ? 1 : -1;
}

/** Input: FILE *output, const Sequence *sequences,                                        int count, AnnoModel *model,                                        const char *bedgraph. Output: unsigned long result. */
static unsigned long write_decoded_gff(FILE *output, const Sequence *sequences,
                                       int count, AnnoModel *model,
                                       const char *bedgraph) {
  BedgraphWriter writer = {0};
  if (bedgraph) {
    if (mkdir(bedgraph, 0777) && errno != EEXIST)
      fail("cannot create bedGraph directory %s: %s", bedgraph, strerror(errno));
    struct stat status;
    if (stat(bedgraph, &status) || !S_ISDIR(status.st_mode))
      fail("bedGraph path is not a directory: %s", bedgraph);
    for (int class = 0; class < ANNO_CLASSES; class++) {
      size_t size = strlen(bedgraph) + strlen(class_names[class]) + 13;
      char *path = checked_calloc(size, 1);
      snprintf(path, size, "%s/%s.bedgraph", bedgraph, class_names[class]);
      writer.files[class] = fopen(path, "w");
      if (!writer.files[class]) fail("cannot write %s", path);
      fprintf(writer.files[class], "track type=bedGraph name=\"%s probability\"\n",
              class_names[class]);
      free(path);
    }
  }
  unsigned long genes = 0;
  for (int sequence_index_value = 0; sequence_index_value < count;
       sequence_index_value++) {
    const Sequence *sequence = &sequences[sequence_index_value];
    float *input = normalized_features(sequence, &model->metadata);
    float *probability = checked_calloc((size_t)sequence->length * ANNO_CLASSES,
                                          sizeof(*probability));
    uint8_t *path = checked_calloc((size_t)sequence->length, 1);
    if (anno_model_predict(model, input, sequence->length, probability) ||
        !anno_grammar_decode(probability, sequence->length, 9, 3, path))
      fail("UNet decoding failed for %s", sequence->name);
    if (writer.files[0]) {
      for (int class = 0; class < ANNO_CLASSES; class++) {
        int64_t start = 0;
        float value = -1.0f;
        for (int64_t position = 0; position <= sequence->length; position++) {
          float next = position < sequence->length
            ? roundf(probability[(size_t)position * ANNO_CLASSES +
                                 (size_t)class] *
                     100000000.0f) / 100000000.0f
            : -1.0f;
          if (position == 0) {
            value = next;
            continue;
          }
          if (position < sequence->length && next == value) continue;
          fprintf(writer.files[class], "%s\t%lld\t%lld\t%.8f\n",
                  sequence->name, (long long)start, (long long)position, value);
          start = position;
          value = next;
        }
      }
    }
    for (int64_t position = 0; position < sequence->length;) {
      if (!path_is_transcript(path[position], '+') &&
          !path_is_transcript(path[position], '-')) {
        position++;
        continue;
      }
      char strand = path_is_transcript(path[position], '+') ? '+' : '-';
      int64_t start = position;
      while (position < sequence->length &&
             path_is_transcript(path[position], strand)) position++;
      int result = write_decoded_block(output, sequence, path, start, position,
                                       strand, genes + 1);
      if (result < 0) fail("cannot write GFF for %s", sequence->name);
      if (result > 0) genes++;
    }
    free(path);
    free(probability);
    free(input);
  }
  for (int class = 0; class < ANNO_CLASSES; class++)
    if (writer.files[class] && fclose(writer.files[class]))
      fail("cannot close bedGraph output");
  return genes;
}

/** Input: const Sequence *sequence,                                   const AnnoModelMetadata *metadata. Output: float * result. */
static float *normalized_features(const Sequence *sequence,
                                  const AnnoModelMetadata *metadata) {
  const int channels = 4 + ANNO_CWT_CHANNELS;
  float *input = checked_calloc((size_t)sequence->length * channels,
                                sizeof(*input));
  for (int64_t position = 0; position < sequence->length; position++) {
    int base = base_index(sequence->sequence[position]);
    for (int channel = 0; channel < 4; channel++)
      input[(size_t)position * (size_t)channels + (size_t)channel] =
        base == channel;
    for (int channel = 0; channel < ANNO_CWT_CHANNELS; channel++)
      input[(size_t)position * (size_t)channels + 4 +
            (size_t)channel] =
        (sequence->cwt[(size_t)position * ANNO_CWT_CHANNELS + (size_t)channel] -
         metadata->feature_mean[4 + channel]) /
        metadata->feature_std[4 + channel];
  }
  return input;
}

/** Input: Sequence *sequences, int count,                              AnnoModelMetadata *metadata. Output: none; output parameters are updated or released. */
static void compute_metadata(Sequence *sequences, int count,
                             AnnoModelMetadata *metadata) {
  long double sums[7] = {0}, squares[7] = {0};
  uint64_t total = 0;
  for (int i = 0; i < count; i++) {
    for (int64_t position = 0; position < sequences[i].length; position++) {
      int base = base_index(sequences[i].sequence[position]);
      total++;
      for (int channel = 0; channel < 7; channel++) {
        float value = channel < 4
          ? (float)(base == channel)
          : sequences[i].cwt[(size_t)position * ANNO_CWT_CHANNELS +
                             (size_t)channel - 4];
        sums[channel] += value;
        squares[channel] += (long double)value * value;
      }
    }
  }
  for (int channel = 0; channel < 7; channel++) {
    long double samples = (long double)total;
    long double mean = total ? sums[channel] / samples : 0.0;
    long double variance =
      total ? squares[channel] / samples - mean * mean : 1.0;
    if (variance < 1e-12) variance = 1.0;
    metadata->feature_mean[channel] = (float)mean;
    metadata->feature_std[channel] = sqrtf((float)variance);
  }
}

/** Input: const Options *options, Sequence *sequences,                               int count, const char *annotation. Output: AnnoModel * result. */
static AnnoModel *train_model(const Options *options, Sequence *sequences,
                              int count, const char *annotation) {
  FeatureList features = read_features(annotation, sequences, count);
  mark_labels(sequences, count, &features);
  feature_free(&features);
  AnnoUnetConfig config = {7, ANNO_CLASSES, 8, 3, 3};
  AnnoModel *model = anno_model_create(&config);
  if (!model) fail("cannot create 7-state UNet");
  compute_metadata(sequences, count, &model->metadata);
  for (int i = 0; i < count; i++) {
    float *input = normalized_features(sequences + i, &model->metadata);
    int64_t training_length = sequences[i].length;
    int training_batch_size = training_length > INT_MAX
      ? INT_MAX : (int)training_length;
    if (anno_model_train(model, input, sequences[i].labels, training_length,
                         options->iterations, training_batch_size, 0.05f,
                         (uint64_t)i + 1))
      fail("training failed for %s", sequences[i].name);
    free(input);
  }
  if (anno_model_save(model, options->model))
    fail("cannot save model %s", options->model);
  fprintf(stderr, "Trained 7-state UNet on %d contigs\n", count);
  return model;
}

/** Input: const char *sequence, int64_t length, char *out. Output: none; output parameters are updated or released. */
static void reverse_complement(const char *sequence, int64_t length, char *out) {
  for (int64_t i = 0; i < length; i++)
    out[i] = complement(sequence[length - i - 1]);
  out[length] = '\0';
}

/** Input: const char *fasta, const char *directory. Output: none; output parameters are updated or released. */
static void run_cwt(const char *fasta, const char *directory) {
  int count = 0;
  Sequence *sequences = read_fasta(fasta, &count);
  if (mkdir(directory, 0777) && errno != EEXIST)
    fail("cannot create CWT directory %s: %s", directory, strerror(errno));
  FILE *matrix = NULL, *index = NULL;
  char *matrix_path = checked_calloc(strlen(directory) + 12, 1);
  char *index_path = checked_calloc(strlen(directory) + 14, 1);
  sprintf(matrix_path, "%s/matrix.bin", directory);
  sprintf(index_path, "%s/contigs.tsv", directory);
  matrix = fopen(matrix_path, "wb");
  index = fopen(index_path, "w");
  if (!matrix || !index) fail("cannot create CWT output in %s", directory);
  fprintf(index, "# anno-cwt-v2\n# scales\t4,5,6\nseqid\tlength\tplus_offset\tminus_offset\n");
  size_t rows = 0;
  for (int i = 0; i < count; i++) {
    size_t values = (size_t)sequences[i].length * ANNO_CWT_CHANNELS * 2;
    float *complex_values = checked_calloc(values, sizeof(*complex_values));
    char *reverse = checked_calloc((size_t)sequences[i].length + 1, 1);
    cwt_values(sequences[i].sequence, sequences[i].length, complex_values);
    reverse_complement(sequences[i].sequence, sequences[i].length, reverse);
    float *reverse_values = checked_calloc(values, sizeof(*reverse_values));
    cwt_values(reverse, sequences[i].length, reverse_values);
    fprintf(index, "%s\t%lld\t%zu\t%zu\n", sequences[i].name,
            (long long)sequences[i].length, rows, rows + (size_t)sequences[i].length);
    if (fwrite(complex_values, sizeof(float), values, matrix) != values ||
        fwrite(reverse_values, sizeof(float), values, matrix) != values)
      fail("cannot write CWT matrix");
    rows += (size_t)sequences[i].length * 2;
    free(reverse_values);
    free(reverse);
    free(complex_values);
  }
  if (fclose(matrix) || fclose(index)) fail("cannot close CWT output");
  free(matrix_path);
  free(index_path);
  sequence_free(sequences, count);
}

/** Input: const char *text, const char *name. Output: int result. */
static int positive_integer(const char *text, const char *name) {
  char *end = NULL;
  errno = 0;
  long value = strtol(text, &end, 10);
  if (errno || end == text || *end || value < 1 || value > INT_MAX)
    fail("invalid %s: %s", name, text);
  return (int)value;
}

/** Input: const char *program. Output: none; output parameters are updated or released. */
static void usage(const char *program) {
  fprintf(stderr,
    "Usage: %s [-m model] [-i iterations] [-p threads] [-b bedgraph_dir] "
    "<genome.fa[.gz]> [annotation.gff3[.gz]]\n"
    "       %s [-p threads] cwt <genome.fa[.gz]> <directory>\n",
    program, program);
}

/** Input: void. Output: int result. */
static int self_test(void) {
  AnnoUnetConfig config = {7, ANNO_CLASSES, 4, 2, 3};
  AnnoModel *model = anno_model_create(&config);
  if (!model) return 0;
  float input[17 * 7] = {0};
  float probability[17 * ANNO_CLASSES];
  uint8_t labels[17];
  for (int position = 0; position < 17; position++) {
    input[position * 7 + position % 4] = 1.0f;
    input[position * 7 + 4 + position % ANNO_CWT_CHANNELS] =
      (float)(position + 1) / 17.0f;
    labels[position] = (uint8_t)(position % 4);
  }
  if (anno_model_predict(model, input, 17, probability)) {
    anno_model_destroy(model);
    return 0;
  }
  int variable = 0;
  for (int i = 0; i < 17; i++) {
    float sum = 0.0f;
    for (int class = 0; class < ANNO_CLASSES; class++)
      sum += probability[i * ANNO_CLASSES + class];
    if (i && fabsf(probability[i * ANNO_CLASSES] - probability[0]) > 1e-7f)
      variable = 1;
    if (!isfinite(sum) || fabsf(sum - 1.0f) > 1e-5f) {
      anno_model_destroy(model);
      return 0;
    }
  }
  if (!variable ||
      anno_model_train(model, input, labels, 17, 10, 17, 0.1f, 1) ||
      anno_model_save(model, "/tmp/anno-7-state-test.model")) {
    anno_model_destroy(model);
    return 0;
  }
  AnnoModel *loaded = NULL;
  int ok = !anno_model_load(&loaded, "/tmp/anno-7-state-test.model");
  remove("/tmp/anno-7-state-test.model");
  anno_model_destroy(loaded);
  anno_model_destroy(model);
  if (!ok) return 0;
  uint8_t path[] = {0, 1, 2, 2, 3, 3, 2, 4, 0};
  if (!anno_grammar_transition(0, 2) ||
      !anno_grammar_transition(2, 3) ||
      !anno_grammar_transition(2, 0) ||
      anno_class_for_state(ANNO_STATE_MINUS_INTRON) != ANNO_MINUS_INTRON ||
      path[0] != ANNO_STATE_INTERGENIC)
    return 0;
  const char plus_sequence[] = "ATGAAATAGCCC";
  Segment plus_exon = {0, 12};
  Orf plus_orf;
  if (!find_orf(plus_sequence, 12, &plus_exon, 1, '+', &plus_orf) ||
      plus_orf.transcript_start != 0 || plus_orf.transcript_end != 9)
    return 0;
  const char minus_sequence[] = "GGGCTATTTCATGGG";
  Segment minus_exon = {3, 12};
  Orf minus_orf;
  if (!find_orf(minus_sequence, 15, &minus_exon, 1, '-', &minus_orf) ||
      minus_orf.transcript_start != 0 || minus_orf.transcript_end != 9)
    return 0;
  const char spliced_sequence[] = "ATGAAACCCCTAGCCC";
  Segment spliced_exons[] = {{0, 6}, {10, 16}};
  Orf spliced_orf;
  if (!find_orf(spliced_sequence, 16, spliced_exons, 2, '+', &spliced_orf) ||
      spliced_orf.transcript_end != 9)
    return 0;
  uint8_t decoded[12];
  for (int i = 0; i < 12; i++) decoded[i] = ANNO_STATE_PLUS_EXON;
  Sequence sequence = {.name = "self", .sequence = (char *)plus_sequence,
                       .length = 12};
  FILE *gff = tmpfile();
  if (!gff || write_decoded_block(gff, &sequence, decoded, 0, 12, '+', 1) != 1) {
    if (gff) fclose(gff);
    return 0;
  }
  rewind(gff);
  char line[1024];
  int records = 0;
  while (fgets(line, sizeof(line), gff)) records++;
  fclose(gff);
  return records >= 4;
}

/** Input: int argc, char **argv. Output: int result. */
int main(int argc, char **argv) {
  if (argc == 2 && !strcmp(argv[1], "--self-test"))
    return self_test() ? EXIT_SUCCESS : EXIT_FAILURE;
  Options options = {
    .model = "anno.model", .iterations = ANNO_DEFAULT_ITERATIONS, .threads = 1
  };
  int argument = 1;
  while (argument < argc && argv[argument][0] == '-') {
    const char *option = argv[argument++];
    if (!strcmp(option, "-h") || !strcmp(option, "--help")) {
      usage(argv[0]);
      return EXIT_SUCCESS;
    }
    if (argument >= argc) {
      usage(argv[0]);
      return EXIT_FAILURE;
    }
    if (!strcmp(option, "-m") || !strcmp(option, "--model")) {
      options.model = argv[argument++];
    } else if (!strcmp(option, "-b") || !strcmp(option, "--bedgraph")) {
      options.bedgraph = argv[argument++];
    } else if (!strcmp(option, "-i") || !strcmp(option, "--iterations")) {
      options.iterations = positive_integer(argv[argument++], "iterations");
    } else if (!strcmp(option, "-p") || !strcmp(option, "--threads")) {
      options.threads = positive_integer(argv[argument++], "threads");
    } else {
      usage(argv[0]);
      return EXIT_FAILURE;
    }
  }
  if (argument < argc && !strcmp(argv[argument], "cwt")) {
    if (options.bedgraph || argc - argument != 3) {
      usage(argv[0]);
      return EXIT_FAILURE;
    }
    run_cwt(argv[argument + 1], argv[argument + 2]);
    return EXIT_SUCCESS;
  }
  int remaining = argc - argument;
  if (remaining < 1 || remaining > 2) {
    usage(argv[0]);
    return EXIT_FAILURE;
  }
  (void)options.threads;
  int count = 0;
  Sequence *sequences = read_fasta(argv[argument], &count);
  for (int i = 0; i < count; i++) {
    sequences[i].cwt = checked_calloc((size_t)sequences[i].length *
                                      ANNO_CWT_CHANNELS, sizeof(float));
    cwt_power(sequences[i].sequence, sequences[i].length, sequences[i].cwt);
  }
  AnnoModel *model = remaining == 2
    ? train_model(&options, sequences, count, argv[argument + 1]) : NULL;
  if (!model && anno_model_load(&model, options.model))
    fail("cannot read 7-state model %s", options.model);
  puts("##gff-version 3");
  unsigned long genes = write_decoded_gff(stdout, sequences, count, model,
                                          options.bedgraph);
  if (fflush(stdout) || ferror(stdout)) fail("cannot write GFF3");
  fprintf(stderr, "Predicted %lu genes\n", genes);
  anno_model_destroy(model);
  sequence_free(sequences, count);
  return EXIT_SUCCESS;
}
