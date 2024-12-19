#ifndef PSPATH_HPP
#define PSPATH_HPP

#include <stdint.h>
#include <vector>

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
  uint8_t *seq;      // sequence
};

struct cluster_t {
  std::vector<sfs_t> specifics; // TODO: kvec
  int va = -1, vb = -1;         // starting and ending vertices
  int offa = -1, offb = -1;     // offsets on the two vertices
  uint64_t ka = -1, kb = -1;    // starting and ending kmers
};

#endif
