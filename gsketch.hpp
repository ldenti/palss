#ifndef GSKETCH_HPP
#define GSKETCH_HPP

#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <zlib.h>

#include "kseq.h"

#include "path.h"
#include "utils.h"

KSTREAM_INIT(gzFile, gzread, 65536)

using namespace std;

inline uint64_t encode(int v, int off, int unique) {
  uint64_t e = ((uint64_t)v << 17) | (off & 0xFFFF) << 1 | (unique & 1);
  return e;
}
inline int8_t decode_unique(uint64_t e) { return e & 1; }
inline int64_t decode_v(uint64_t e) { return (e >> 17); }
inline int16_t decode_off(uint64_t e) { return (e >> 1) & 0xFFFF; }

struct seg_t {
  int idx;   // identifier (as in gfa)
  char *seq; // sequence
  int l;     // actual length
  int c;     // capacity
};

struct link_t {
  int idx1;
  int idx2;
};

class GSK {
public:
  uint8_t klen;
  int nvertices;                  // number of vertices
  map<uint64_t, uint64_t> sketch; // sketch : (47[vertex, 16[offset, 1[unique)
  vector<vector<int>> graph;
  vector<path_t *> paths;    // TODO: sampling
  map<int, string> vertices; // make this an array assuming indices in [0,n]

  GSK(char *, uint8_t);
  int build_graph();
  void destroy_graph();
  int build_sketch();
  pair<int64_t, int16_t> get(uint64_t &);
  vector<int> adj(int);
  vector<path_t *> get_subpaths(int, int);
  int get_sequence(const path_t *, char **, int *);
  int compatible(int, int);
  int get_vl(int);
  int store(char *);

private:
  char *gfa_fn;
  void gfa_parse_S(char *, seg_t *);
  void gfa_parse_L(char *, link_t *);
  void gfa_parse_P(char *, path_t *);
  void gfa_parse_W(char *, path_t *);

  void get_s();
  void add_kmer(uint64_t, uint64_t, uint16_t);
};

#endif
