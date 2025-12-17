#include "misc.hpp"

// reverse&complement a 0123-encoded string
void rc3(char *s, int l) {
  int i;
  for (i = 0; i < (l >> 1); ++i) {
    int tmp = s[l - 1 - i];
    s[l - 1 - i] = 4 - s[i];
    s[i] = 4 - tmp;
  }
  if (l & 1)
    s[i] = 4 - s[i];
}

// reverse&complement a 1234-encoded string
void rc4(char *s, int l) {
  int i;
  for (i = 0; i < (l >> 1); ++i) {
    int tmp = s[l - 1 - i];
    s[l - 1 - i] = 5 - s[i];
    s[i] = 5 - tmp;
  }
  if (l & 1)
    s[i] = 5 - s[i];
}

double realtime() {
  struct timeval tp;
  struct timezone tzp;
  gettimeofday(&tp, &tzp);
  return (double)tp.tv_sec + (double)tp.tv_usec * 1e-6;
}