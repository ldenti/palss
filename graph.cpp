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

  // rt = realtime();
  // uint npaths = gbz.index.metadata.paths();
  // plens = sdsl::int_vector<0>(npaths, 0, 32);
  // sdsl::simple_sds::load_from(plens, fn + ".pl");
  // fprintf(stderr, "[M::%s] restored path lengths in %.3f sec\n", __func__,
  //         realtime() - rt);

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

positions_t Graph::get_positions() const {
  positions_t output;
  uint npaths = gbz.index.metadata.paths();
  for (uint pp = 0; pp < npaths; ++pp) {
    int p = gbwt::Path::encode(pp, 0);
    // here we need the path identifier without strand
    std::string sample_name = gbz.index.metadata.fullPath(pp).sample_name;
    std::string contig_name = gbz.index.metadata.fullPath(pp).contig_name;
    // size_t haplotype = gbz.index.metadata.fullPath(pp).haplotype;
    // std::cerr << sample_name << " " << contig_name << " " << haplotype
    //           << std::endl;
    int cp = 0;
    // here we need the encoded path identifier (with strand)
    if (sample_name.compare("_gbwt_ref") == 0) {
      gbwt::vector_type path = gbz.index.extract(p);
      // std::cerr << path.size() << std::endl;
      for (const gbwt::node_type v : path) {
        // v has strand information
        gbwtgraph::nid_t vv = gbwt::Node::id(v);
        // vv is internal node_id
        gbwtgraph::handle_t hh = gbz.graph.get_handle(vv);
        // string name = gbz.graph.get_segment_name(hh);
        int l = gbz.graph.get_length(hh);
        output.offsets[vv] = cp;
        output.v2ref[vv] = contig_name;
        cp += l;
      }
      output.references[contig_name] = cp;
    }
  }
  return output;
}
