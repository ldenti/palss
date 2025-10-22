#ifndef PS_READS_HPP
#define PS_READS_HPP

#include <stdio.h>
#include <zlib.h>

#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

typedef struct {
  char *name;
  size_t name_l, name_m;
  char *seq;
  size_t seq_l, seq_m;
} read_t;

typedef struct {
  gzFile fp;
  read_t **reads;
  kseq_t *seq;
  size_t n;
  size_t m;
} rbatch_t;

read_t *rd_init();
void rd_load(read_t *r, kseq_t *seq);
void rd_destroy(read_t *r);

rbatch_t *rbx_init(const char *fn, int m);
int rbx_load(rbatch_t *rb);
void rbx_destroy(rbatch_t *rb);

#endif