#ifndef RFFT_H
#define RFFT_H
// Reasonably Fast Fourier Transform (in public domain)

#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>

// Perform the FFT in place for an array of size 2^k.
// No normalization is done.
void fft_transform(double complex* vec, size_t n, bool inverse);

#endif // RFFT_H

#ifdef RFFT_IMPLEMENTATION
#ifndef M_PI
	#define M_PI 3.14159265358979323846
#endif

void fft_transform(double complex* vec, size_t n, bool inverse) {
	int levels = 0;	 // Compute levels = floor(log2(n))
	for (size_t k =	1; (k &	n) == 0; k <<= 1)
		levels++;

	// Permute vec by reversing the bits of addresses
	for (size_t i =	0; i < n; i++) {
		// Reverse the bits of i
		size_t j = 0, ii = i;
		for (int k = 0;	k < levels; k++) {
			j = (j << 1) | (ii & 1);
			ii >>= 1;
		}

		if (j >	i) {
			double complex tmp = vec[i];
			vec[i] = vec[j];
			vec[j] = tmp;
		}
	}

	// Cooley-Tukey	in place
	for (size_t half = 1; half < n;	half *=	2) {
		size_t size = 2	* half;
		size_t step = n	/ size;
		for (size_t j =	0, k = 0; j < half; j++, k += step) {
			double angle = (inverse	? 2 : -2) * M_PI * k / n;
			double complex omega = cos(angle) + I *	sin(angle);
        		for (size_t i =	0; i < n; i += size) {
				double complex tmp = vec[i + j + half] * omega;
				vec[i + j + half] = vec[i + j] - tmp;
				vec[i + j] += tmp;
			}
		}
	}
}

#endif // RFFT_IMPLEMENTATION