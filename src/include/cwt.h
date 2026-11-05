/**
* @file cwt.h
*/
#ifndef CWT_H
#define CWT_H

#include <stdint.h>

#define CWT_CHANNELS 3 // Number of CWT scales
#define CWT_WINDOW 1024 // Length of parallely processed CWT window
#define CWT_RADIUS 4 // Radius of scale units of Morlet complex wavelet
#define MORLET_CENTRAL 6.0 // Central frequency of Morlet complex wabvvelet

/**
* @brief Morlet CWT of the complex base signal for positions [start,end).
* @param[in] bases: Input sequence bases.
* @param[in] length: Total sequence length.
* @param[in] start: First position, inclusive.
* @param[in] end: Last position, exclusive.
* @param[out] out: Receives tensor of (end-start) x CWT_CHANNELS x 2 floats (real, imaginary).
* @return None.
*/
void cwt(const char *bases, int64_t length, int64_t start, int64_t end,
    float *out);

/**
 * @brief Compute CWT power in position-major layout.
 * @param[in] bases Sequence bases.
 * @param[in] length Total sequence length.
 * @param[in] start Inclusive interval start.
 * @param[in] end Exclusive interval end.
 * @param[out] out Receives the interval's CWT power channels.
 * @return None.
 */
void cwt_power(const char *bases, int64_t length, int64_t start,
                    int64_t end, float *out);



#endif /** CWT_H */
