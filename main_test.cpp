#include <filesystem>
#include <fstream>
#include <getopt.h>
#include <string>

#include "gbwt/fast_locate.h"
#include "gbwt/gbwt.h"
#include "gbwtgraph/gbz.h"
#include "kmc_api/kmc_file.h"
#include "sdsl/simple_sds.hpp"

#include "kmer.hpp"
#include "misc.hpp"
#include "sketch.hpp"
#include "usage.hpp"

typedef struct {
  uint64_t kmer;
  uint32_t vertex; // gbwt identifier, without strand information
} solid_anchor_t;

// std::string print_kmer(uint64_t kmer_d, int klen) {
//   char *kmer = (char *)malloc(sizeof(char) * (klen + 1));
//   d23(kmer_d, klen, kmer);
//   for (int i = 0; i < klen; ++i)
//     kmer[i] = "0ACGT0"[kmer[i]];
//   return kmer;
// }

int main_test(int argc, char *argv[]) {
  std::string gbz_fn = argv[1];

  gbwtgraph::GBZ gbz;
  sdsl::simple_sds::load_from(gbz, gbz_fn);
  gbwt::printStatistics(gbz.index, gbz_fn, std::cerr);

  // R-index
  gbwt::FastLocate fl;
  fl = gbwt::FastLocate(gbz.index);
  std::ofstream out;
  out.open(gbz_fn + ".ri", std::ofstream::out);
  fl.serialize(out);
  out.close();

  const gbwt::Metadata &metadata = gbz.index.metadata;
  int ret;
  for (gbwt::Metadata::size_type sample_id = 0; sample_id < metadata.samples();
       ++sample_id) {

    std::string sample_name = metadata.sample(sample_id);
    fprintf(stderr, "[M::%s] %s (%ld)\n", __func__, sample_name.c_str(),
            sample_id);

    std::vector<gbwt::size_type> paths = metadata.pathsForSample(sample_id);
    for (size_t pp = 0; pp < paths.size(); ++pp) {
      gbwt::size_type path_id = paths[pp];
      gbwt::size_type seq_id = gbwt::Path::encode(path_id, false);

      gbwt::edge_type curr = gbz.index.start(seq_id);
      while (curr.first != gbwt::ENDMARKER) {
        // curr.first
        gbwtgraph::handle_t handle =
            gbwtgraph::GBWTGraph::node_to_handle(curr.first);
        std::string gfa_idx = gbz.graph.get_segment_name(handle);
        gbwtgraph::view_type view = gbz.graph.get_sequence_view(handle);
        std::string seq(view.first, view.second);
        std::cerr << curr.first << " " << ((curr.first & 1) == 0 ? '+' : '-')
                  << " " << gfa_idx << " " << seq << std::endl;

        uint32_t v = (curr.first >> 1);

        std::cerr << gbz.graph.get_segment_name(
                         gbwtgraph::GBWTGraph::node_to_handle(v << 1))
                  << std::endl;

        curr = gbz.index.LF(curr);
      }
    }
  }
  return 0;
}