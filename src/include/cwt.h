#ifndef CWT_H
#define CWT_H

#include "surf.h"

#define CWT_CHANNELS 3
#define CWT_RADIUS 4
#define MORLET_CENTRAL 6.0

void cwt_power(const char *bases, int64_t length, int64_t start, int64_t end,
               float *out);

#endif
