#ifndef PS_GRAPH_HPP
#define PS_GRAPH_HPP

#include <map>
#include <stdio.h>
#include <string>
#include <vector>

#include "gbwt/fast_locate.h"
#include "gbwt/gbwt.h"
#include "gbwtgraph/gbz.h"
#include "gbwtgraph/utils.h"
#include "sdsl/simple_sds.hpp"

#include <misc.hpp>

typedef struct {
  gbwt::size_type id;                    // with strand
  std::vector<gbwt::size_type> vertices; // with strand
  std::string sequence;
  gbwt::size_type dist_to_end;
  gbwt::size_type offset1;
  gbwt::size_type offset2;
  size_t cutpfx, cutsfx;
} path_t;
bool operator==(const path_t &path1, const path_t &path2);
bool operator<(const path_t &path1, const path_t &path2);
void p_reverse(path_t &path);

// All functions expect vertices without strand information

/*
 * if vertex is too long, it gets split. So we have different curr.first but
 * still same gfa_idx
 */

class Graph {
private:
  std::string fn;
  std::string reference;
  gbwtgraph::GBZ gbz;
  gbwt::FastLocate fl;

public:
  Graph(const std::string &fn, const std::string &reference = "");
  int load();
  int load_fl();
  int build_fl();
  const gbwt::Metadata &get_metadata() const;
  const gbwtgraph::GBZ &get_gbz() const;
  const gbwt::FastLocate &get_fl() const;
  //   int distance(gbwtgraph::nid_t a, gbwtgraph::nid_t b) const;
  //   std::map<gbwt::size_type, std::vector<gbwt::size_type>>
  //   locate(gbwtgraph::nid_t v) const;
  void print_stats() const;
  path_t get_path(gbwt::size_type path_id) const;
  // Get all paths from v1 to v2
  std::vector<path_t> get_paths(uint32_t v1, uint32_t v2,
                                bool ref_only = false) const;
  std::string get_gfa_name(uint32_t v) const;
  size_t get_vertex_len(gbwt::size_type v) const;
  std::string get_path_contig(gbwt::size_type pid) const;
  // std::string get_path_sample(gbwt::size_type pid) const;
  std::map<std::string, size_t> get_reference_paths() const;
  std::string
  get_path_sequence(const std::vector<gbwt::node_type> &vertices) const;
};

#endif
