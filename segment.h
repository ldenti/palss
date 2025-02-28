#ifndef PS_SEG_H
#define PS_SEG_H

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "kvec.h"

typedef struct {
  int idx;   // identifier (as in gfa) - assuming integer
  char *seq; // sequence
  int l;     // length
  int c;     // capacity
  kvec_t(int) paths;
} seg_t;

/* Initialize a segment */
seg_t *init_seg();

/* Destroy the segment */
void destroy_seg(seg_t *seg);

#endif
