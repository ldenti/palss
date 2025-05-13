#ifndef PS_GRAPH_HPP
#define PS_GRAPH_HPP

#include <stdio.h>
#include <string>

#include "gbwt/fast_locate.h"
#include "gbwt/gbwt.h"
#include "gbwtgraph/gbz.h"
#include "gbwtgraph/utils.h"
#include "sdsl/simple_sds.hpp"

#include <misc.hpp>

// TODO: adapt this and add to index subcommand

class Graph {
private:
  std::string fn;
  gbwt::FastLocate fl;
  sdsl::int_vector<0> plens;

public:
  // XXX: make this private (breaks main_map)
  gbwtgraph::GBZ gbz;

  Graph(const std::string &fn);
  int load();
  int distance(gbwtgraph::nid_t a, gbwtgraph::nid_t b) const;
  std::map<gbwt::size_type, std::vector<gbwt::size_type>>
  locate(gbwtgraph::nid_t v) const;
  void print_stats() const;
};

#endif
