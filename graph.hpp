#ifndef PS_GRAPH_HPP
#define PS_GRAPH_HPP

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include <map>
#include <set>
#include <string>
#include <vector>

#include "kseq.h"

#include "misc.hpp"
// #include "segments.hpp"

// Some assumptions I made:
//   - vertex identifiers in GFA are integers
//   - graph has less than INT_MAX nodes
//   - we do not need to store all paths here

typedef struct {
  int idx;         // identifier (as in gfa) - assuming integer
  std::string seq; // sequence
  int l;           // length
  int path = -1;   // reference path
  int pos = -1;    // position on reference path
} seg_t;

typedef struct {
  // identifier (as in gfa)
  std::string idx;
  // identifiers of the vertices along the path (graph space)
  std::vector<int> vertices;
} path_t;

class Graph {
public:                        // XXX: private
  std::string fn;              // gfa file name
  std::vector<seg_t> vertices; // labels in graph vertex order
  // segments_t *vertices;
  std::map<int, int>
      v_map; // mapping between gfa idx and internal idx (position on vertices)

  int ne;                                  // total number of edges
  std::vector<std::vector<int>> out_edges; // outgoing edges
  std::vector<std::vector<int>> in_edges;  // incoming edges

  std::vector<path_t> paths;

  std::map<int, uint64_t> in_distances;
  std::map<int, uint64_t> out_distances;

public:
  Graph(const std::string &fn);

  // Load all vertices from gfa
  int load_vertices();

  // Load all edges from gfa
  int load_edges();

  // Load all paths from gfa
  int load_paths(const std::string &ref);

  // Build distance index
  int build_distance_index();

  // iterative bfs from vertex v, direction out
  uint64_t bfs(int v, int out) const;

  // Get distance between two vertices
  uint distance(int v1, int v2);

  // Convert vertex from GFA space to graph space
  int get_iidx(int v) const;

  // int distance(gbwtgraph::nid_t a, gbwtgraph::nid_t b) const;
  // std::map<gbwt::size_type, std::vector<gbwt::size_type>>
  // locate(gbwtgraph::nid_t v) const;
  // void print_stats() const;

  // // Get subpaths going from v1 to v2
  // std::map<gbwt::size_type, gbwt::vector_type>
  // get_subpaths(gbwtgraph::nid_t v1, gbwtgraph::nid_t v2) const;
  // // Get positions of vertices on reference paths
  // positions_t get_positions() const;
  std::string get_sequence(int v) const;
  // std::string get_gfa_idx(gbwtgraph::nid_t v) const;
};

// /* Initialize a path with capacity c*/
// path_t *init_path();
// /* Add vertex v to path, reallocating it if needed */
// void add_vertex(path_t *path, int v);
// /* Destroy a path */
// void destroy_path(path_t *p);

// GFA reading utilities
void gfa_parse_S(char *s, seg_t &ret);
void gfa_parse_L(char *s, int *idx1, int *idx2);
void gfa_parse_P(char *s, path_t &p, const std::string &ref);
void gfa_parse_W(char *s, path_t &p, const std::string &ref);

#endif
