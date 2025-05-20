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
  // std::cerr << " === " << v << " ===" << std::endl;
  // gbwtgraph::handle_t hh = gbz.graph.get_handle(v);
  // std::string name = gbz.graph.get_segment_name(hh);
  // int l = gbz.graph.get_length(hh);
  // std::cerr << name << " " << l << " " << gbz.graph.get_sequence(hh)
  //           << std::endl;

  std::map<gbwt::size_type, std::vector<gbwt::size_type>> result;

  for (int strand = 0; strand < 2; ++strand) {
    int vv = gbwt::Node::encode(v, strand);
    if (!gbz.index.contains(vv))
      continue;
    std::vector<gbwt::size_type> xxx = fl.decompressSA(vv);
    for (const gbwt::size_type &x : xxx) {
      gbwt::size_type pidx = fl.seqId(x);
      gbwt::size_type pp = gbwt::Path::id(pidx);
      int is_reverse = gbwt::Path::is_reverse(pidx);
      if (!is_reverse)
        // here we check on - strand so that distance from end is distance from
        // start on + strand
        continue;
      int o = fl.seqOffset(x);
      // int o = is_reverse ? fl.seqOffset(x) : pl - fl.seqOffset(x) - 1;
      // std::cerr << "  " << o << std::endl;
      result[pp].push_back(o);
    }
  }

  for (auto &[_, val] : result) {
    std::sort(val.begin(), val.end());
  }
  return result;
}

// size_t Graph::get_nsamples() const { return gbz.index.da_samples.size(); }

std::string Graph::get_gfa_idx(gbwtgraph::nid_t v) const {
  return gbz.graph.get_segment_name(gbz.graph.get_handle(v));
}

std::map<gbwt::size_type, gbwt::vector_type>
Graph::get_subpaths(gbwtgraph::nid_t v1, gbwtgraph::nid_t v2) const {
  // XXX: what about cycles?
  // XXX: what if we want all paths?

  // std::cerr << v1 << ">" << v2 << std::endl;
  std::map<gbwt::size_type, gbwt::vector_type> paths;

  for (int strand1 = 0; strand1 < 2; ++strand1) {
    for (int strand2 = 0; strand2 < 2; ++strand2) {
      int vv1 = gbwt::Node::encode(v1, strand1);
      int vv2 = gbwt::Node::encode(v2, strand2);
      std::vector<gbwt::size_type> pvisits1 = fl.decompressSA(vv1);
      std::vector<gbwt::size_type> pvisits2 = fl.decompressSA(vv2);

      for (int i = 0; i < pvisits1.size(); ++i) {
        for (int j = 0; j < pvisits2.size(); ++j) {
          gbwt::size_type pidx1 = fl.seqId(pvisits1[i]);
          gbwt::size_type pidx2 = fl.seqId(pvisits2[j]);
          if (pidx1 != pidx2)
            continue;
          gbwt::size_type pp = gbwt::Path::id(pidx1);
          std::string sample_name = gbz.index.metadata.fullPath(pp).sample_name;
          if (sample_name.compare("_gbwt_ref") != 0)
            continue;
          if (fl.seqOffset(pvisits1[i]) < fl.seqOffset(pvisits2[j]))
            continue;

          gbwt::edge_type position = std::make_pair(vv1, i);
          //   if (!gbz.index.contains(position)) {
          //     vv1 = gbwt::Node::encode(v1, 1);
          //     position = std::make_pair(vv1, k);
          //     assert(gbz.index.contains(position));
          //   }

          gbwt::vector_type path;
          while (position.first != vv2) {
            // std::cout << (position.first >> 1) << " ";
            path.push_back(position.first);
            position = gbz.index.LF(position);
          }
          // std::cout << (position.first >> 1) << std::endl;
          path.push_back(position.first);
          paths[pp] = path;
        }
      }
    }
  }
  // std::map<gbwt::size_type, std::vector<gbwt::size_type>> occs1 =
  // locate(v1); std::map<gbwt::size_type, std::vector<gbwt::size_type>> occs2
  // = locate(v2);

  // for (auto &[k, val1] : occs1) {
  //   // k is the path identifier without strand
  //   std::string sample_name = gbz.index.metadata.fullPath(k).sample_name;
  //   if (sample_name.compare("_gbwt_ref") != 0)
  //     continue;
  //   if (occs2.find(k) == occs2.end())
  //     continue;
  //   auto val2 = occs2[k];
  //   if (val1.size() > 1 && val2.size() > 1) {
  //     std::cerr << "Skipping " << v1 << ":" << v2 << " (cycle)" <<
  //     std::endl; continue;
  //   }

  //   // element of path are node identifier + strand information
  //   // gbwt::vector_type path0 = gbz.index.extract(gbwt::Path::encode(k,
  //   0));
  //   // paths[k] = {path0.begin() + val1[0], path0.begin() + val2[0] + 1};

  //   // TODO: LF from v2 to v1 along path k (how to relate k to offset from
  //   // searchstate?)

  //   gbwt::vector_type path;
  //   gbwt::node_type vv1 = gbwt::Node::encode(v1, 0);
  //   gbwt::SearchState searchstate = gbz.index.find(vv1);
  //   for (gbwt::size_type offset = searchstate.range.first;
  //        offset <= searchstate.range.second; ++offset) {
  //     std::cerr << "--- " << offset << " " << k << std::endl;
  //   }

  //   gbwt::node_type vv1m = gbwt::Node::encode(v1, 1);
  //   searchstate = gbz.index.find(vv1m);
  //   for (gbwt::size_type offset = searchstate.range.first;
  //        offset <= searchstate.range.second; ++offset) {
  //     std::cerr << "--- " << offset << " " << k << std::endl;
  //   }

  //   exit(1);
  //   gbwt::edge_type position = std::make_pair(vv1, k);
  //   if (!gbz.index.contains(position)) {
  //     vv1 = gbwt::Node::encode(v1, 1);
  //     position = std::make_pair(vv1, k);
  //     assert(gbz.index.contains(position));
  //   }

  //   while ((position.first >> 1) != v2) {
  //     // std::cout << (position.first >> 1) << " ";
  //     path.push_back(position.first);
  //     position = gbz.index.LF(position);
  //   }
  //   paths[k] = path;
  //   // std::cout << (position.first >> 1) << " ";
  // }
  return paths;
}

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

std::string Graph::get_sequence(gbwtgraph::nid_t v) const {
  return gbz.graph.get_sequence(gbz.graph.get_handle(v));
}
