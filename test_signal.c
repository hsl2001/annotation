#include "anno_signal.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(void) {
  float peaks[] = {0, 0, 0, 10, 10, 10, NAN};
  assert(fabs(otsu_threshold(peaks, 7) - 10.0 / 256) < 1e-12);
  float constant[] = {2, 2, NAN}, missing[] = {NAN, NAN};
  assert(isinf(otsu_threshold(constant, 3)));
  assert(isinf(otsu_threshold(missing, 2)));
  assert(isinf(otsu_threshold(NULL, 0)));
  float asymmetric[] = {0, 1, 2, 8, 9, 10};
  double threshold = otsu_threshold(asymmetric, 6);
  assert(threshold > 2 && threshold < 8);
  assert(base_signal('A') == 1.0 && base_signal('C') == I);
  assert(base_signal('G') == -I && base_signal('T') == -1.0);
  assert(base_signal('N') == 0.0 && base_signal('c') == I);
  Wavelets *wavelets = anno_alloc(1, sizeof(*wavelets));
  wavelets_init(wavelets);
  for (int scale = 0; scale < WAVE_COUNT; scale++) {
    double complex sum = 0.0;
    double energy = 0.0;
    for (int tap = 0; tap < wave_sizes[scale]; tap++) {
      double magnitude = cabs(wavelets->kernel[scale][tap]);
      energy += magnitude * magnitude;
      sum += wavelets->kernel[scale][tap];
    }
    assert(cabs(sum) < 1e-12 && fabs(energy - 1.0) < 1e-12);
  }
  double first[CWT_CHANNELS];
  cwt_extract(wavelets, "NNNANNNNN", 9, 3, 1, first);
  for (int scale = 0; scale < WAVE_COUNT; scale++) {
    double complex expected = wavelets->kernel[scale][wave_sizes[scale] / 2];
    assert(fabs(first[2 * scale] - creal(expected)) < 1e-12);
    assert(fabs(first[2 * scale + 1] - cimag(expected)) < 1e-12);
  }
  char *reverse = reverse_complement("aCGTN", 5);
  assert(strcmp(reverse, "NACGT") == 0);
  free(reverse);
  free(wavelets);
  puts("signal tests passed");
  return 0;
}