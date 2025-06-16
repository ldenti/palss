#ifndef PS_SFS_HPP
#define PS_SFS_HPP

#include <string>

typedef struct {
  int64_t v;    // vertex on graph
  int offset;   // offset on vertex
  int p;        // position on query
  uint64_t seq; // kmer
  int closest;  // is this anchor the closest one to the specific string?
} anchor_t;

typedef struct {
  std::string qname; // read name
  int s;             // start on read
  int l;             // length
  anchor_t a;        // left anchor
  anchor_t b;        // right anchor
  int strand;        // inferred strand
  bool good = true;  // is it good for calling step?
  uint8_t *seq;      // 1-4encoded sequence
} sfs_t;

// anchor_t parse_anchor(char *line);
// sfs_t parse_sfs_line(char *line);

#endif
