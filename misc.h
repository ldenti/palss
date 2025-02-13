#ifndef PS_MISC_H
#define PS_MISC_H

#include <sys/time.h>

#ifndef MIN
#define MIN(a, b) ((a) < (b) ? a : b)
#endif

static inline double realtime();

#endif
