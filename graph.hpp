#ifndef PS_GRAPH_HPP
#define PS_GRAPH_HPP

#include <map>
#include <stdio.h>
#include <string>

#include "gbwt/fast_locate.h"
#include "gbwt/gbwt.h"
#include "gbwtgraph/gbz.h"
#include "gbwtgraph/utils.h"
#include "sdsl/simple_sds.hpp"

#include <misc.hpp>

// TODO: adapt this and add to index subcommand

typedef struct {
  std::map<std::string, int> references;
  std::map<int, std::string> v2ref;
  std::map<int, int> offsets;
} positions_t;

class Graph {
private:
  std::string fn;
  gbwt::FastLocate fl;
  sdsl::int_vector<0> plens;
  gbwtgraph::GBZ gbz;

public:
  Graph(const std::string &fn);
  int load();
  int distance(gbwtgraph::nid_t a, gbwtgraph::nid_t b) const;
  std::map<gbwt::size_type, std::vector<gbwt::size_type>>
  locate(gbwtgraph::nid_t v) const;
  void print_stats() const;

  // Get subpaths going from v1 to v2
  std::map<gbwt::size_type, gbwt::vector_type>
  get_subpaths(gbwtgraph::nid_t v1, gbwtgraph::nid_t v2) const;
  // Get positions of vertices on reference paths
  positions_t get_positions() const;
  std::string get_sequence(gbwtgraph::nid_t v) const;
  std::string get_gfa_idx(gbwtgraph::nid_t v) const;
};

#endif
