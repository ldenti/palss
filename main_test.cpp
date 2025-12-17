// #include <filesystem>
// #include <fstream>
// #include <getopt.h>
#include <string>

#include "graph.hpp"
// #include "kmer.hpp"
// #include "misc.hpp"
// #include "sketch.hpp"
// #include "usage.hpp"

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

// Root/leaf path in the tree rooted at a given vertex
typedef struct {
  std::vector<gbwt::node_type> vertices;
  size_t l; // length in bp after first vertex (root)
} rlpath2_t;

// TODO : we could improve this using some sort of caching
std::vector<rlpath2_t> kdfs2(const gbwtgraph::GBZ &gbz,
                             const gbwt::FastLocate &fl,
                             const gbwt::node_type &source,
                             const size_t &klen) {
  std::vector<rlpath2_t> paths; // fully extended path
  std::queue<rlpath2_t> queue;  // paths we still need to "extend"
  queue.push({{source}, 0});
  while (!queue.empty()) {
    rlpath2_t path = queue.front();
    queue.pop();
    gbwt::node_type vertex = path.vertices.back();
    gbwtgraph::handle_t handle = gbz.graph.node_to_handle(vertex);

    // get outgoings edges
    std::vector<gbwt::node_type> outs;
    gbwt::SearchState state = gbz.graph.get_state(handle);
    gbz.graph.follow_paths(
        state, [&outs](const gbwt::SearchState &next_state) -> bool {
          if (!next_state.empty())
            outs.push_back(next_state.node);
          return true;
        });

    bool extended = false;
    for (const gbwt::node_type &v : outs) {
      const gbwtgraph::handle_t &h = gbz.graph.node_to_handle(v);
      size_t l = gbz.graph.get_length(h);
      rlpath2_t new_path = path;
      new_path.vertices.push_back(v);

      // check if path is consistent with haplotypes
      gbwt::size_type first;
      gbwt::SearchState state =
          fl.find(new_path.vertices.begin(), new_path.vertices.end(), first);
      if (state.empty()) {
        continue;
      }
      extended = true;
      // add path to solutions or queue depending on sequence length
      new_path.l += l;
      if (new_path.l >= (size_t)klen - 1)
        paths.push_back(new_path);
      else
        queue.push(new_path);
    }
    if (!extended) {
      // add path to solutions since we cannot extend it, only if we added at
      // least one vertex to it
      if (path.vertices.size() > 1)
        paths.push_back(path);
    }
  }
  return paths;
}

int main_test(int argc, char *argv[]) {
  std::string gbz_fn = argv[1];
  size_t klen = std::stoi(argv[2]);

  Graph graph(gbz_fn, "CHM13");
  graph.load();
  graph.build_fl();

  const gbwtgraph::GBZ &gbz = graph.get_gbz();
  const gbwt::FastLocate &fl = graph.get_fl();
  const gbwt::Metadata &metadata = graph.get_metadata();

  // Get all anchors by visiting each vertex (kDFS)
  size_t total_vertices = gbz.graph.max_node_id() - gbz.graph.min_node_id();
  std::cerr << "---" << gbz.graph.min_node_id() << " "
            << gbz.graph.max_node_id() << std::endl;

  for (gbwtgraph::nid_t source_id = gbz.graph.min_node_id();
       source_id <= gbz.graph.max_node_id(); ++source_id) {

    /* note 1: segments longer than 1024 are split, so max_node_id >= #S
     * note 2: not all IDs in [min_node_id, max_node_id] are actually vertices
     * note 3: we need to check both strand for a node since outgoing edges
     *         depends on node strand
     * note 4: we will limit to haplotype paths on the + strand
     **/

    for (int strand = 0; strand < 2; ++strand) {
      gbwtgraph::handle_t source_handle =
          gbz.graph.get_handle(source_id, strand);

      std::cout << source_id << " " << (strand ? "-" : "+")
                << gbz.graph.get_segment_name(source_handle) << " "
                << gbz.graph.get_length(source_handle) << std::endl;
      continue;

      gbwt::node_type source = gbz.graph.handle_to_node(source_handle);
      // iterate over all paths starting from this vertex (that are
      // consistent with indexed haplotypes)

      std::vector<rlpath2_t> paths = kdfs2(gbz, fl, source, klen);

      for (const rlpath2_t &p : paths) {
        for (const auto &v : p.vertices) {
          std::cout << graph.get_gfa_name(v >> 1) << " ";
        }
        std::cout << std::endl;
      }

      for (const rlpath2_t &path : paths) {
        // keep only paths on + strand
        gbwt::size_type first;
        gbwt::SearchState state =
            fl.find(path.vertices.begin(), path.vertices.end(), first);
        assert(!state.empty());

        std::vector<gbwt::size_type> p_offsets = fl.locate(state, first);
        bool on_plus = false;
        bool is_ref = false;
        //     // for (const auto &v : path.vertices)
        //     //   std::cerr << graph.get_gfa_name(v >> 1) << ((v & 1) ? "-" :
        //     "+")
        //     //             << " ";
        //     // std::cerr << std::endl;
        for (const gbwt::size_type &p : p_offsets) {
          if (gbwt::Path::is_reverse(p))
            continue;
          std::cerr << metadata.fullPath(gbwt::Path::id(p)).contig_name << " "
                    << metadata.fullPath(gbwt::Path::id(p)).offset << std::endl;

          auto f = metadata.findFragment(metadata.path(p));

          std::string path_sequence;
          // build path sequence and path info (vertex for each position)
          for (const gbwt::node_type &v : path.vertices) {
            gbwtgraph::handle_t h = gbwtgraph::GBWTGraph::node_to_handle(v);
            gbwtgraph::view_type view = gbz.graph.get_sequence_view(h);
            path_sequence.append(view.first, view.second);
          }
          std::cerr << path_sequence << std::endl;
        }
      }
    }
  }

  return 0;
}
