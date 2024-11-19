#ifndef PSPATH_HPP
#define PSPATH_HPP

#include <stdint.h>

struct anchor_t {
  int64_t v = -1;    // vertex on graph
  int offset = -1;   // offset on vertex
  int p = -1;        // position on query
  uint64_t seq = -1; // kmer
};

struct sfs_t {
  int qidx;        // read index
  int s;           // start on query
  int l;           // length
  anchor_t a = {}; // left anchor
  anchor_t b = {}; // right anchor
  int strand = 1;  // inferred strand
  uint64_t esk = -1,
           eek = -1; // expected starting and ending kmers (from cluster)
  int good = 1;      // is it good for calling step?
  char *seq = NULL;  // sequence
};

#endif
