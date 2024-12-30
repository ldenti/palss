#ifndef PSSFS_H
#define PSSFS_H

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"

typedef struct anchor_t {
  int64_t v = -1;    // vertex on graph
  int offset = -1;   // offset on vertex
  int p = -1;        // position on query
  uint64_t seq = -1; // kmer
} anchor_t;

typedef struct sfs_t {
  int qidx;        // read index
  int s;           // start on query
  int l;           // length
  anchor_t a = {}; // left anchor
  anchor_t b = {}; // right anchor
  int strand = 1;  // inferred strand
  uint64_t esk = -1,
           eek = -1; // expected starting and ending kmers (from cluster)
  int good = 1;      // is it good for calling step?
  char *rname;
  uint8_t *seq; // sequence
} sfs_t;

anchor_t parse_anchor(char *line);
sfs_t parse_sfs_line(char *line);

#endif
