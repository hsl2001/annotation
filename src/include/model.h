#ifndef MODEL_H
#define MODEL_H

#include "cwt.h"
#include "surf.h"

enum { GS, IS, IE, GE, SITES, JUNCTION = SITES, SPAN, CLASSIFIERS };
extern const char *const classifier_names[CLASSIFIERS];
#ifndef SITE_RADIUS
#define SITE_RADIUS 15
#endif
#if SITE_RADIUS < 0
#error "SITE_RADIUS must be nonnegative"
#endif
#define SITE_WIDTH (2 * SITE_RADIUS + 1)
#define SITE_FEATURES (SITE_WIDTH * (NUCLEOTIDES + 1 + 64) + 1)
#define PAIR_FEATURES (SITE_WIDTH * (2 * NUCLEOTIDES + 1) + 3)
#define MODEL_FEATURES                                                         \
  (SITE_FEATURES > PAIR_FEATURES ? SITE_FEATURES : PAIR_FEATURES)

typedef struct {
  union {
    struct {
      uint8_t bases[SITE_WIDTH];
      float power[SITE_WIDTH];
    };
    struct {
      uint8_t ends[2][SITE_WIDTH];
      float distance[2];
    };
  };
  uint8_t type, target;
  float weight;
} Sample;

typedef struct {
  double weights[CLASSIFIERS][MODEL_FEATURES];
  int64_t max_intron, max_span;
} Model;

int site_candidate(const char *bases, int64_t length, int64_t p, int type);
void sample_site(const Sequence *sequence, const float *power, int64_t origin,
                 int64_t p, int type, Sample *sample);
void sample_pair(const Sequence *sequence, int64_t start, int64_t end, int type,
                 Sample *sample);
double site_logit(const Model *model, const Sample *sample);
double pair_logit(const Model *model, const Sequence *sequence, int64_t start,
                  int64_t end, int type);
Model *train_model(const Sequence *sequences, int count, int iterations,
                   int threads);
unsigned long predict(const Model *model, const Sequence *sequences, int count,
                      FILE *gff, const char *prefix);
int save_model(const char *path, const Model *model);
int load_model(const char *path, Model *model);

#endif
