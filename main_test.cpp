// #include <filesystem>
// #include <fstream>
// #include <getopt.h>
#include <string>

#include "graph.hpp"
#include "reads.hpp"
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

typedef struct {
  std::deque<gbwt::node_type> vertices;
  std::string sequence;              // sequence after first vertex
  std::vector<gbwt::node_type> info; // vertex for each character in sequence
} subpath_t;

std::vector<subpath_t> kdfs_hap(const gbwtgraph::GBZ &gbz, subpath_t path,
                                int klen) {
  std::vector<subpath_t> paths; // fully extended path
  std::queue<subpath_t> queue;  // paths we still need to "extend"
  queue.push(path);
  while (!queue.empty()) {
    subpath_t sp = queue.front();
    queue.pop();
    gbwt::node_type vertex = sp.vertices.back();
    gbwtgraph::handle_t handle = gbz.graph.node_to_handle(vertex);

    // get all outgoings states
    std::vector<gbwt::SearchState> outs;
    gbwt::SearchState state = gbz.graph.get_state(handle);
    gbz.graph.follow_paths(
        state, [&outs](const gbwt::SearchState &next_state) -> bool {
          if (!next_state.empty())
            outs.push_back(next_state);
          return true;
        });

    for (const gbwt::SearchState &s : outs) {
      std::cerr << "### " << s.node << std::endl;
      std::vector<gbwt::size_type> pvisits = gbz.index.locate(s);
      for (const auto &p : pvisits)
        std::cerr << p << " ";
      std::cerr << std::endl;

      const gbwtgraph::handle_t &h = gbz.graph.node_to_handle(s.node);
      subpath_t sp2;
      sp2.vertices = sp.vertices;
      sp2.vertices.push_back(s.node);
      sp2.sequence = sp.sequence + gbz.graph.get_sequence(h);
      sp2.info = sp.info;
      for (size_t i = 0; i < gbz.graph.get_length(h); ++i)
        sp2.info.push_back(s.node);
      if (sp2.sequence.size() >= (size_t)klen - 1)
        paths.push_back(sp2);
      else
        queue.push(sp2);
    }
  }
  return paths;
}

std::string sp_print_2(const subpath_t &sp) {
  std::string s = "[";
  for (const auto &v : sp.vertices)
    s += " " + std::to_string(v >> 1) + (v & 1 ? "-" : "+");
  s += " ] : " + sp.sequence;
  return s;
}





// ###################################################################

// Root/leaf path in the tree rooted at a given vertex
typedef struct {
  std::vector<gbwt::node_type> vertices;
  size_t l; // length in bp after first vertex (root)
} rlpath_t;

std::vector<rlpath_t> kdfs_test(const gbwtgraph::GBZ &gbz,
                           const gbwt::FastLocate &fl,
                           const gbwt::node_type &source, const size_t &klen) {
  std::vector<rlpath_t> paths; // fully extended path
  std::queue<rlpath_t> queue;  // paths we still need to "extend"
  queue.push({{source}, 0});
  while (!queue.empty()) {
    rlpath_t path = queue.front();
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

    for (const gbwt::node_type &v : outs) {
      const gbwtgraph::handle_t &h = gbz.graph.node_to_handle(v);
      size_t l = gbz.graph.get_length(h);
      rlpath_t new_path = path;
      new_path.vertices.push_back(v);

      // check if path is consistent with haplotypes
      gbwt::size_type first;
      gbwt::SearchState state =
          fl.find(new_path.vertices.begin(), new_path.vertices.end(), first);
      if (state.empty())
        continue;

      // add path to solutions or queue depending on sequence length
      new_path.l += l;
      if (new_path.l >= (size_t)klen - 1)
        paths.push_back(new_path);
      else
        queue.push(new_path);
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
  for (gbwtgraph::nid_t source_id = gbz.graph.min_node_id();
       source_id <= gbz.graph.max_node_id(); ++source_id) {

    /* note 1: segments longer than 1024 are split, so max_node_id >= #S
     * note 2: not all IDs in [min_node_id, max_node_id] are actually vertices
     * note 3: we need to check both strand for a node since outgoing edges
     *         depends on node strand
     **/

    for (int strand = 1; strand >= 0; --strand) {
      gbwtgraph::handle_t source_handle =
	gbz.graph.get_handle(source_id, strand);

      gbwt::node_type source = gbz.graph.handle_to_node(source_handle);
      // gbwtgraph::view_type source_view =
      //     gbz.graph.get_sequence_view(source_handle);
      // size_t source_length = source_view.second;
      // std::string source_sequence(source_view.first, source_view.second);
      // std::vector<gbwt::node_type> source_info(source_sequence.size(), source);

      std::vector<rlpath_t> paths = kdfs_test(gbz, fl, source, klen);
      for (const rlpath_t &path : paths) {
	for (const auto &v : path.vertices) {
	  gbwtgraph::handle_t h = gbwtgraph::GBWTGraph::node_to_handle(v);
	  std::cout << gbz.graph.get_segment_name(h) << (v & 1 ? "-" : "+") << " ";
	}
	std::cout << std::endl;
      }
    }
  }
  return 0;
}

int main_test_cacheddfs(int argc, char *argv[]) {
  std::string gbz_fn = argv[1];
  int klen = std::stoi(argv[2]);

  Graph graph(gbz_fn, "CHM13");
  graph.load();
  graph.build_fl();
  const gbwtgraph::GBZ &gbz = graph.get_gbz();
  const gbwt::Metadata &metadata = graph.get_metadata();

  std::vector<bool> visited(gbz.graph.get_node_count(), false);

  for (size_t path_id = 0; path_id < metadata.paths(); ++path_id) {
    gbwt::size_type seq_id = gbwt::Path::encode(path_id, false);
    gbwt::edge_type curr = gbz.index.start(seq_id);
    while (curr.first != gbwt::ENDMARKER) {
      gbwtgraph::handle_t handle =
          gbwtgraph::GBWTGraph::node_to_handle(curr.first);
      gbwtgraph::nid_t handle_id = gbz.graph.get_id(handle);
      if (visited[handle_id - 1]) {
      } else {
        std::cerr << "Path " << path_id << " "
                  << gbz.graph.get_segment_name(handle)
                  << ((curr.first & 1) ? "-" : "+") << " - visiting"
                  << std::endl;
        subpath_t subpath = {{curr.first}, "", {}};
        std::vector<subpath_t> paths = kdfs_hap(gbz, subpath, klen);
        for (const auto &p : paths) {
          std::cout << sp_print_2(p) << std::endl;
        }
      }
      curr = gbz.index.LF(curr);
      visited[handle_id - 1] = true;
    }
  }

  return 0;
}

int main_test_rbatch(int argc, char *argv[]) {
  char *fn = argv[1];
  rbatch_t *rb = rbx_init(fn, 10000);
  int x;
  while ((x = rbx_load(rb)) > 0) {
    printf("%s\n", rb->reads[0]->name);
  }
  rbx_destroy(rb);

  return 0;
}
