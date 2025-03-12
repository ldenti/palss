#ifndef PS_SEG_H
#define PS_SEG_H

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ic.h"
#include "kvec.h"
#include "rle.h"

typedef struct {
  int idx;   /* identifier (as in gfa) - assuming integer */
  char *seq; /* sequence */
  int l;     /* length */
  int c;     /* capacity */

  /* Everything we need to build the rle bv of paths */
  uint8_t *paths;
  int paths_p;
  int paths_c;
  int paths_b; /* this is the same as *rle_nptr(paths) */
  int lasth;
  int lastp;
  int rl1;
  int64_t cnts[6];
  /* --- */

  kvec_t(int) starts;   /* Offset of each path over pord */
  kvec_t(int) pord;     /* Path-ordering */
  unsigned char *cpord; /* Compressed path-ordering */
} seg_t;

/* Initialize a segment */
seg_t *init_seg();

/* Update segment with path h and order o */
void update_seg(seg_t *seg, int h, int o);

/* Finalize rle and compress path-ordering */
void compress_seg(seg_t *seg);

/* Destroy the segment */
void destroy_seg(seg_t *seg);

#endif
