#ifndef PS_ANCHOR_HPP
#define PS_ANCHOR_HPP

#include <stdint.h>
#include <vector>
#include <map>

typedef struct {
  // uint32_t id;      // path identifier w/ strand
  uint32_t offset1; // offset along path
  uint32_t offset2; // offset along path
  bool is_reference;
  // bool is_reverse;
} pp_t;

typedef struct {
  // 0123-encoded kmer
  uint64_t kmer;
  // gbwt identifier w/ strand
  uint32_t v1;
  uint32_t v2;
  // position along vertices (on "their" strand)
  uint32_t pos1;
  uint32_t pos2;
  //
  bool is_reference; // anchor is on reference path (can be also on other paths)
  bool is_valid;     // anchor is valid, so not repeated
  //
  int qp; // position on query
  //
  std::map<uint32_t, pp_t> paths; // paths, key is path identifier w/ strand
  //
  bool is_canonical; // canonical is this version (eg for AAA is true, for TTT
                     // is false)
  bool has_both;     // we saw this anchor on both "strand" (++/--, +-/-+)
} anchor_t;

static inline bool operator<(const anchor_t &x, const anchor_t &y) {
  return x.kmer < y.kmer;
}

typedef std::vector<anchor_t> anchors_t;

#endif