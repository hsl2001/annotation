/**
 * @file cwt.h
 */
#ifndef CWT_H
#define CWT_H

#include "surf.h"

#define CWT_CHANNELS 3 // Number of CWT scales

#define CWT_RADIUS 4       // Radius of scale units of Morlet complex wavelet
#define MORLET_CENTRAL 6.0 // Central frequency of Morlet complex wabvvelet
#define INPUT_CHANNELS                                                         \
  (NUCLEOTIDES +                                                               \
   CWT_CHANNELS) // Number of input channels

/**
 * @brief Compute CWT power in position-major layout.
 * @param[in] bases: Sequence bases.
 * @param[in] length: Total sequence length.
 * @param[in] start: Inclusive interval start.
 * @param[in] end: Exclusive interval end.
 * @param[out] out: Receives the interval's CWT power channels.
 * @return None.
 */
void cwt_power(const char *bases, int64_t length, int64_t start, int64_t end,
               float *out);

/**
 * @brief Fill the nucleotide and CWT power features for one position.
 */
void compute_features(char base, const float *power, float *features);

#endif /** CWT_H */
