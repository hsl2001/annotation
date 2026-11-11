#ifndef MODEL_H
#define MODEL_H

#include "cwt.h"
#include "surf.h"

enum { TSS, IS, IE, TTS, SITES, JUNCTION = SITES, CLASSIFIERS };
extern const char *const classifier_names[CLASSIFIERS];
#ifndef SITE_RADIUS
#define SITE_RADIUS 15
#endif
#if SITE_RADIUS < 0
#error "SITE_RADIUS must be nonnegative"
#endif
#define SITE_WIDTH (2 * SITE_RADIUS + 1)
#define SITE_FEATURES (SITE_WIDTH * (NUCLEOTIDES + 1 + 64) + 1)
#define JUNCTION_FEATURES (SITE_WIDTH * (2 * NUCLEOTIDES + 1) + 3)
#define MODEL_FEATURES SITE_FEATURES

typedef struct {
  union {
    struct {
      uint8_t bases[SITE_WIDTH];
      float power[SITE_WIDTH];
    };
    struct {
      uint8_t donor[SITE_WIDTH], acceptor[SITE_WIDTH];
      float distance[2];
    };
  };
  uint8_t type, target;
  float weight;
} Sample;

typedef struct {
  double weights[CLASSIFIERS][MODEL_FEATURES];
  int64_t max_intron, max_exon;
} Model;

int site_candidate(const char *bases, int64_t length, int64_t p, int type);
void sample_site(const Sequence *sequence, const float *power, int64_t origin,
                 int64_t p, int type, Sample *sample);
void sample_junction(const Sequence *sequence, int64_t donor, int64_t acceptor,
                     Sample *sample);
double site_logit(const Model *model, const Sample *sample);
double junction_logit(const Model *model, const Sequence *sequence,
                      int64_t donor, int64_t acceptor);
Model *train_model(const Sequence *sequences, int count, int iterations,
                   int threads);
unsigned long predict(const Model *model, const Sequence *sequences, int count,
                      FILE *gff, const char *prefix);
int save_model(const char *path, const Model *model);
int load_model(const char *path, Model *model);

#endif
