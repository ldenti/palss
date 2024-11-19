#ifndef PATH_HPP
#define PATH_HPP

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  int idx;   // identifier (as in gfa)
  char *seq; // sequence
  int l;     // actual length
  int c;     // capacity
} seg_t;

typedef struct {
  int idx1;
  int idx2;
} link_t;

seg_t *init_seg();
void destroy_seg(seg_t *seg);
void gfa_parse_S(char *s, seg_t *ret);

#endif
