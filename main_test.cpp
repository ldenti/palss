// #include <filesystem>
// #include <fstream>
// #include <getopt.h>
#include <string>

// #include "graph.hpp"
// #include "kmer.hpp"
// #include "misc.hpp"
// #include "sketch.hpp"
// #include "usage.hpp"

#include "reads.hpp"

// void print_path_visits(const gbwt::FastLocate &fl, uint32_t v) {
//   const gbwt::Metadata &metadata = fl.index->metadata;

//   // std::vector<gbwt::size_type> offsets = fl.locate(fl.index->find(v));
//   std::vector<gbwt::size_type> offsets = fl.decompressSA(v);
//   for (const auto &off : offsets) {
//     gbwt::size_type path_id = fl.seqId(off) >> 1; // remove strand
//     information std::string sample_name =
//     metadata.fullPath(path_id).sample_name; std::string contig_name =
//     metadata.fullPath(path_id).contig_name;
//     // int haplotype = metadata.fullPath(path_id).haplotype;
//     // int offset = metadata.fullPath(path_id).offset;
//     std::cerr << "--- " << path_id << " " << sample_name << " " <<
//     contig_name
//               << std::endl;
//   }
//   std::cerr << std::endl;
// }

int main_test(int argc, char *argv[]) {
  // std::string gbz_fn = argv[1];
  // Graph graph(gbz_fn);
  // graph.load();
  // graph.load_fl();
  // uint32_t x = 274900;
  // uint32_t y = 562222;
  // std::vector<path_t> paths = graph.get_paths(x, y);
  // std::cerr << paths.size() << std::endl;

  char *fn = argv[1];

  rbatch_t *rb = rbx_init(fn, 10000);

  int x;
  while ((x = rbx_load(rb)) > 0) {
    printf("%s\n", rb->reads[0]->name);
  }

  rbx_destroy(rb);

  return 0;
}