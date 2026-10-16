#ifndef RFFT_H
#define RFFT_H
// Reasonably Fast Fourier Transform (in public domain)

#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

// Perform the FFT in place for an array of size 2^k.
// No normalization is done.
void fft_transform_radix2(double _Complex* vec, size_t n, bool inverse);
void fft_transform_radix2_twiddle(double _Complex* restrict vec, size_t n, bool inverse,
                                  const double _Complex* restrict twiddle);

// Perform the FFT for an arbitrary array using the Bluestein's algorithm.
// The code needs to allocate two supplementary buffers (of size < 4 * n - 3).
// No normalization is done.
void fft_transform_bluestein(double _Complex* vec, size_t n, bool inverse);

// Perform the FFT, choosing the suitable algorithm from the two above.
void fft_transform(double _Complex* vec,	size_t n, bool inverse);

#endif // RFFT_H

#ifdef RFFT_IMPLEMENTATION
#ifndef RFFT_CALLOC
	#define RFFT_CALLOC(n, size)	calloc(n, size)
	#define RFFT_FREE(p)	free(p)
#endif

#ifndef M_PI
	#define M_PI 3.14159265358979323846
#endif

void fft_transform_radix2_twiddle(double complex* restrict vec, size_t n, bool inverse,
                                  const double complex* restrict twiddle) {
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
		for (size_t i =	0; i < n; i += size) {
			if (!inverse) {
				for (size_t j =	0, k = 0; j < half; j++, k += step) {
					double complex tmp = vec[i + j + half] * twiddle[k];
					vec[i + j + half] = vec[i + j] - tmp;
					vec[i + j] += tmp;
				}
			} else {
				for (size_t j =	0, k = 0; j < half; j++, k += step) {
					double complex tmp = vec[i + j + half] * conj(twiddle[k]);
					vec[i + j + half] = vec[i + j] - tmp;
					vec[i + j] += tmp;
				}
			}
		}
	}
}

void fft_transform_radix2(double complex* vec, size_t n, bool inverse) {
	if (n <= 1) return;
	double complex twiddle_stack[4096];
	double complex* twiddle = (n / 2 <= 4096)
		? twiddle_stack : RFFT_CALLOC(n / 2, sizeof(double complex));
	double base_angle = -2.0 * M_PI / n;
	for (size_t k = 0; k < n / 2; k++) {
		double angle = base_angle * k;
		twiddle[k] = cos(angle) + I * sin(angle);
	}
	fft_transform_radix2_twiddle(vec, n, inverse, twiddle);
	if (twiddle != twiddle_stack) RFFT_FREE(twiddle);
}

void fft_transform_bluestein(double complex* vec, size_t n, bool inverse) {
	// Find m = 2^k such that m >= 2 * n + 1
	size_t m = 1;
	while (m <= 2 * n) {
		m *= 2;
	}

	// Extended vectors of size m
	double complex*	avec = RFFT_CALLOC(m, sizeof(double complex));
	double complex*	bvec = RFFT_CALLOC(m, sizeof(double complex));
	avec[0]	= vec[0];
	bvec[0]	= 1.0;
	for (size_t i =	1; i < n; i++) {
		size_t k = (i * i) % (2 * n);
		double angle = (inverse	? M_PI : -M_PI)	* k / n;
		double complex omega = cos(angle) + I *	sin(angle);
		avec[i]	= vec[i] * omega;
		bvec[i]	= bvec[m - i] =	conj(omega);
	}

	// Convolution
	fft_transform_radix2(avec, m, false);
	fft_transform_radix2(bvec, m, false);
	for (size_t i =	0; i < m; i++) {
		avec[i]	*= bvec[i];
	}
	fft_transform_radix2(avec, m, true);

	for (size_t i =	0; i < n; i++) {
		size_t k = (i * i) % (2 * n);
		double angle = (inverse	? M_PI : -M_PI)	* k / n;
		double complex omega = cos(angle) + I *	sin(angle);
		vec[i] = avec[i] * omega / m;
	}

	RFFT_FREE(avec);
	RFFT_FREE(bvec);
}

void fft_transform(double complex* vec,	size_t n, bool inverse)	{
	if (n <= 1)
		return;
	else if	((n & (n - 1)) == 0)  // Power of 2
		fft_transform_radix2(vec, n, inverse);
	else
		fft_transform_bluestein(vec, n,	inverse);
}

#endif // RFFT_IMPLEMENTATION