#ifndef GSKETCH_HPP
#define GSKETCH_HPP

/*#include <boost/graph/adjacency_list.hpp>*/
#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <zlib.h>

#include "kseq.h"

KSTREAM_INIT(gzFile, gzread, 65536)

using namespace std;

static const uint8_t to_int[128] = {0, 0, 1, 2, 3, 0, 0, 0, 0, 0, // 0
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 10
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 20
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 30
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 40
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 50
                                    0, 0, 0, 0, 0, 1, 0, 2, 0, 0, // 60
                                    0, 3, 0, 0, 0, 0, 0, 0, 0, 0, // 70
                                    0, 0, 0, 0, 4, 0, 0, 0, 0, 0, // 80
                                    0, 0, 0, 0, 0, 0, 0, 1, 0, 2, // 90
                                    0, 0, 0, 3, 0, 0, 0, 0, 0, 0, // 100
                                    0, 0, 0, 0, 0, 0, 4, 0, 0, 0, // 110
                                    0, 0, 0, 0, 0, 0, 0, 0};      // 120

inline uint8_t reverse_char(const uint8_t c) { return ((~c) & 3); }

inline uint64_t k2d(char *kmer, uint8_t k) {
  uint64_t kmer_d = 0;
  uint8_t x;
  for (uint8_t i = 0; i < k; ++i) {
    x = (kmer[i] < 6 ? kmer[i] : to_int[kmer[i]]) -
        1; // we assume sequence to be encoded as A:1 but in kmer A is 0
    kmer_d = (kmer_d << 2) | (x < 4 ? x : rand() % 4);
  }
  return kmer_d;
}

inline uint64_t rc(uint64_t kmer, const uint8_t k) {
  uint64_t rckmer = 0;
  kmer = ~kmer;
  for (uint8_t i = 0; i < k; ++i) {
    rckmer = (rckmer << 2) | (kmer & 3);
    kmer >>= 2;
  }
  return rckmer;
}

inline uint64_t lsappend(const uint64_t kmer, const uint64_t c,
                         const uint64_t k) { // left shift and append
  return ((kmer << 2) | c) & ((1UL << 2 * k) - 1);
}

inline uint64_t rsprepend(const uint64_t kmer, const uint64_t c,
                          const uint64_t k) { // right shift and prepend
  return (kmer >> 2) | (c << (2 * k - 2));
}

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

struct path_t {
  char *idx;
  int *vertices;
  int l;
  int capacity;
};

class GSK {
public:
  /*boost::adjacency_list<> graph;*/
  vector<vector<int>> graph;
  vector<path_t *> paths;    // TODO: sampling
  map<int, string> vertices; // make this an array assuming indices in [0,n]
  map<uint64_t, uint64_t> sketch;
  int k;
  int nvertices; // number of vertices
  GSK(char *);
  int build_graph();
  void destroy_graph();
  int build_sketch(int l);
  int get(uint64_t &);
  vector<int> adj(int);
  vector<path_t *> get_subpaths(int, int);
  int get_sequence(const path_t *, char **, int *);
  int compatible(int, int);

private:
  char *gfa_fn;
  set<uint64_t> multi;

  void gfa_parse_S(char *, seg_t *);
  void gfa_parse_L(char *, link_t *);
  void gfa_parse_P(char *, path_t *);
  void gfa_parse_W(char *, path_t *);

  void get_s();
  void add_kmer(uint64_t &, int &);
};

#endif
