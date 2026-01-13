#ifndef PS_ANCHOR_HPP
#define PS_ANCHOR_HPP

#include <map>
#include <stdint.h>
#include <vector>

typedef struct {
  // uint32_t id;      // path identifier w/ strand
  uint32_t offset1;  // offset along path
  uint32_t offset2;  // offset along path
  bool strand;       // is the anchor consistent with the path? (depending on
                     // canonical)
  bool reversed;     // we reversed the vertices of the anchor to find this path
  bool is_reference; // this path is the reference path
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
  int qp;     // position on query
  int qp_rev; // position on query (reverse)
  //
  std::map<uint32_t, pp_t> paths; // paths, key is path identifier w/ strand
  //
  bool is_canonical;
  bool has_both;     // we saw this anchor on both "strand" (++/--, +-/-+)
} anchor_t;

typedef std::vector<anchor_t> anchors_t;

#endif