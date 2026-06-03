#ifndef PS_ANCHOR_HPP
#define PS_ANCHOR_HPP

#include <map>
#include <stdint.h>
#include <vector>

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
  std::map<uint32_t, uint64_t>
      paths; // paths, key is path identifier w/ strand bit + referece bit,
             // consistency bit, and inverted bit. (consistency: "is the kmer
             // along the path consistent with the kmer extracted from the
             // read?", to model strand. "inverted": anchor on this path comes
             // from B-A- and not A+B+). Values are encoded path offsets
  //
  bool inverted; // this anchor is inverted along the path (XXX: I think this is
                 // bad since this is used only when anchor is assigned to a
                 // path during chaining)
  bool is_canonical;
  bool has_both; // we saw this anchor on both "strand" (++/--, +-/-+)
} anchor_t;

typedef std::vector<anchor_t> anchors_t;

#endif