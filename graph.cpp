#include "graph.hpp"

Graph::Graph(const std::string &fn) : fn(fn) {}

int Graph::load() {
  double rt = realtime();
  sdsl::simple_sds::load_from(gbz, fn);
  fprintf(stderr, "[M::%s] restored GBZ in %.3f sec\n", __func__,
          realtime() - rt);

  rt = realtime();

  std::ifstream in;
  in.open(fn + ".ri", std::ifstream::in | std::ifstream::app);
  fl.load(in);
  fl.setGBWT(gbz.index);
  in.close();
  fprintf(stderr, "[M::%s] restored R-index in %.3f sec\n", __func__,
          realtime() - rt);

  rt = realtime();
  uint npaths = gbz.index.metadata.paths();
  plens = sdsl::int_vector<0>(npaths, 0, 32);
  sdsl::simple_sds::load_from(plens, fn + ".pl");
  fprintf(stderr, "[M::%s] restored path lengths in %.3f sec\n", __func__,
          realtime() - rt);

  return 0;
}

void Graph::print_stats() const {
  gbwt::printStatistics(gbz.index, "gbwt");
  gbwt::printStatistics(fl, "fl");
}

int Graph::distance(gbwtgraph::nid_t a, gbwtgraph::nid_t b) const {
  auto Oa = locate(a);
  auto Ob = locate(b);
  uint mind = -1, d;
  for (auto &[key, valA] : Oa) {
    if (Ob.find(key) == Ob.end())
      continue;
    for (const gbwt::size_type &oa : valA) {
      for (const gbwt::size_type &ob : Ob[key]) {
        d = std::max(oa, ob) - std::min(oa, ob);
        if (d < mind)
          mind = d;
      }
    }
  }
  return mind;
}

std::map<gbwt::size_type, std::vector<gbwt::size_type>>
Graph::locate(gbwtgraph::nid_t v) const {
  std::map<gbwt::size_type, std::vector<gbwt::size_type>> result;

  int vv = gbwt::Node::encode(v, 0);
  std::vector<gbwt::size_type> xxx = fl.decompressSA(vv);
  for (const gbwt::size_type &x : xxx) {
    gbwt::size_type pidx = fl.seqId(x);
    int is_reverse = gbwt::Path::is_reverse(pidx);
    gbwt::size_type p = pidx >> 1;
    int pl = plens[p];
    int o = is_reverse ? fl.seqOffset(x) : pl - fl.seqOffset(x) - 1;
    result[p].push_back(o);
  }

  for (auto &[_, val] : result) {
    std::sort(val.begin(), val.end());
  }
  return result;
}

// size_t Graph::get_nsamples() const { return gbz.index.da_samples.size(); }

// std::string Graph::get_gfa_idx(gbwt::node_type v) const {
//   return gbz.graph.get_segment_name(gbz.graph.get_handle(v));
// }
