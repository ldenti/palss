#include <string>

#include "gbwt/gbwt.h"
#include "gbwtgraph/gbz.h"
#include "sdsl/simple_sds.hpp"

int main_extract(int argc, char *argv[]) {
  std::string gbz_fn = argv[1];
  std::string sample_name = argv[2];

  gbwtgraph::GBZ gbz;
  sdsl::simple_sds::load_from(gbz, gbz_fn);
  // gbwt::printStatistics(gbz.index, gbz_fn, std::cerr);

  const gbwt::Metadata &metadata = gbz.index.metadata;
  gbwt::size_type sample_id = metadata.sample(sample_name);
  std::vector<gbwt::size_type> paths = metadata.pathsForSample(sample_id);
  for (size_t pp = 0; pp < paths.size(); ++pp) {
    gbwt::size_type path_id = paths[pp];
    gbwt::size_type seq_id = gbwt::Path::encode(path_id, false);
    // std::string contig_name = gbz.index.metadata.fullPath(path_id).contig_name;
    // int haplotype = gbz.index.metadata.fullPath(path_id).haplotype;
    // int offset = gbz.index.metadata.fullPath(path_id).offset;

    std::string sequence = ""; // full path sequence
    gbwt::edge_type curr = gbz.index.start(seq_id);
    while (curr.first != gbwt::ENDMARKER) {
      gbwtgraph::handle_t handle =
          gbwtgraph::GBWTGraph::node_to_handle(curr.first);
      // view_type: in-place view of the sequence: (start, length)
      gbwtgraph::view_type view = gbz.graph.get_sequence_view(handle);
      sequence.append(view.first, view.second);
      curr = gbz.index.LF(curr);
    }
    // CHM13#chr20
    // HG00733#1#JAHEPQ010000055.1#0[32408859-32512270]
    // printf(">{}#{}#{}#{}\n");
    printf(">%ld\n", pp);
    printf("%s\n", sequence.c_str());
  }

  return 0;
}
