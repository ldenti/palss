#include "misc.h"

double realtime() {
  struct timeval tp;
  struct timezone tzp;
  gettimeofday(&tp, &tzp);
  return (double)tp.tv_sec + (double)tp.tv_usec * 1e-6;
}
