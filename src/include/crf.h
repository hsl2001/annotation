/**
 * @file crf.h
 * @brief Linear-chain CRF
 */
#ifndef CRF_H
#define CRF_H

#include "cwt.h"
#include "surf.h"

#include <stdint.h>

#define GRAMMAR_MAX_BACKPOINTERS 5

/** Parameters of CRF */
#define EMISSIONS (CLASSES * (INPUT_CHANNELS + 1))
#define TRANSITIONS (CLASSES * CLASSES)
#define BOUNDARY_RADIUS 64
#define BOUNDARY_FEATURES ((2 * BOUNDARY_RADIUS + 1) * NUCLEOTIDES)
#define BOUNDARY_WEIGHTS (2 * CLASSES * BOUNDARY_FEATURES)
#define PARAMETERS (EMISSIONS + TRANSITIONS + BOUNDARY_WEIGHTS)

/** Learned emission and transition parameters for the observed classes. */
typedef union {
  struct {
    double emission_weights[EMISSIONS];
    double transition_weights[TRANSITIONS];
    double boundary_weights[BOUNDARY_WEIGHTS];
  };
  double weights[PARAMETERS];
} Model;

/** Just one sampled, compact, and position-major feature window */
typedef struct {
  int64_t length;       // Number of labeled positions
  const uint8_t *bases; // A, C, G, T as 0..3, (4 for ambiguous); CAUTION: other
                        // implemations are 0, 1, 2, 3 and -1
  const float *cwt;     // Position-major CWT power array
  const uint8_t *classes;          // Observed class labels per position
  int64_t offset, sequence_length; // Preserve DNA context across training cuts.
} Sample;

/** Grammar for genes */
typedef struct {
  int count;
  int backpointers[CLASSES];
  int from[CLASSES][GRAMMAR_MAX_BACKPOINTERS];
} Grammar;

typedef struct Decoder Decoder;

void compute_emissions(const Model *model, const float *features,
                       double scores[CLASSES]);

/** Score class transitions at the boundary before position, using DNA only. */
int compute_transitions(const double *bias, const double *pwm,
                        const Sample *sample, size_t position,
                        double scores[TRANSITIONS]);

/**
 * @brief Optimize a CRF model with L-BFGS over labeled windows.
 * @param[in,out] model: Initial model parameters; replaced with the trained
 * parameters on success.
 * @param[in] grammar: State grammar.
 * @param[in] samples: Labeled training windows.
 * @param[in] sample_count Number of training windows.
 * @param[in] iterations: Maximum L-BFGS iterations.
 * @param[in] threads: Number of objective workers.
 * @return An LbfgsStatus value describing convergence or failure.
 */
int train(Model *model, const Grammar *grammar, const Sample *samples,
          size_t sample_count, int iterations, int threads);

/**
 * @brief Train a CRF from labeled sequences using L-BFGS.
 * @param[in] sequences: Sequences with per-base training labels.
 * @param[in] count: Number of sequences.
 * @param[in] iterations: Maximum optimizer iterations.
 * @param[in] threads: Number of parallel objective workers.
 * @return Newly allocated trained CRF model.
 */
Model *train_pipeline(const Sequence *sequences, int count, int iterations,
                      int threads);

/**
 * @brief Check that every CRF parameter is finite.
 * @param[in] model: Model to inspect.
 * @return Nonzero when all parameters are finite.
 */
int is_model_finite(const Model *model);

/** Context shared with the CRF objective callback. */
typedef struct {
  const Grammar *grammar;
  const Sample *samples;
  size_t sample_count;
  uint64_t total_length;
  int threads;
  int failure_status;
} ModelObjectiveContext;

/**
 * @brief Compute normalized conditional NLL and its full parameter gradient.
 * @param[in] x: Current parameter vector.
 * @param[out] gradient: Receives the objective gradient.
 * @param[in,out] opaque: CRF objective context; records any worker failure.
 * @return Objective value, or NaN if evaluation fails.
 */
double model_objective(const double *x, double *gradient, void *opaque);

/**
 * @brief Build the strand-specific grammar.
 * @param[out] grammar: Receives the initialized grammar.
 * @return None.
 */
void build_grammar(Grammar *grammar);

Decoder *decoder_create(int64_t length);
void decoder_destroy(Decoder *decoder);
void decoder_step(Decoder *decoder, const double *emission,
                  const double *transition);
int decoder_finish(Decoder *decoder, uint8_t *path);

/**
 * @brief Predict sequences, writing GFF3 and optional class bedGraphs.
 */
unsigned long predict(const Model *model, const Sequence *sequences, int count,
                      FILE *gff, const char *prefix);

/**
 * @brief Write a versioned CRF model file.
 * @param[in] path: Destination model path.
 * @param[in] model: Finite model parameters to serialize.
 * @return Zero on success, or an errno-style error code.
 */
int save_model(const char *path, const Model *model);

int load_model(const char *path, Model *model);

#endif /* CRF_H */
