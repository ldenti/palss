#include <algorithm>
#include <filesystem>
#include <fstream>
#include <getopt.h>
#include <string>

#include "graph.hpp"
#include "kmer.hpp"
#include "misc.hpp"
#include "sketch.hpp"
#include "usage.hpp"

// Each anchor is a pair <kmer, vertices>
// Vertices are encoded as 32bits for starting vertex, 32bits for ending vertex,
// using gbwt identifier w/ strand
typedef std::vector<std::pair<uint64_t, uint64_t>> anchors_t;

// Root/leaf path in the tree rooted at a given vertex
typedef struct {
  std::vector<gbwt::node_type> vertices;
  size_t l; // length in bp after first vertex (root)
} rlpath_t;

// Store anchors starting in positions [start, end) of sequence. Tag using vinfo
void get_anchors(anchors_t &anchors, const std::string &sequence,
                 const std::vector<gbwt::node_type> &vinfo, const size_t &klen,
                 const size_t &start, const size_t &end) {
  assert(sequence.size() == vinfo.size());
  char kmer[klen + 1];   // first kmer on sequence (plain)
  uint64_t kmer_d = 0;   // kmer
  uint64_t rckmer_d = 0; // reverse and complemented kmer
  uint64_t ckmer_d = 0;  // canonical kmer
  uint8_t c;             // new character to append
  size_t pos;            // current position
  uint64_t value;        // value (start/end vertices)

  // first kmer
  pos = start;
  strncpy(kmer, sequence.c_str() + pos, klen);
  kmer_d = k2d(kmer, klen);
  rckmer_d = rc(kmer_d, klen);
  ckmer_d = std::min(kmer_d, rckmer_d);
  value = (vinfo[pos] << 32) | vinfo[pos + klen - 1];
  anchors.push_back(std::make_pair(ckmer_d, value));

  // all other kmers
  ++pos;
  for (; pos < end; ++pos) {
    c = to_int[(uint8_t)sequence[pos + klen - 1]] - 1; // A is 1 but we want 0
    kmer_d = lsappend(kmer_d, c, klen);
    rckmer_d = rsprepend(rckmer_d, reverse_char(c), klen);
    ckmer_d = std::min(kmer_d, rckmer_d);
    value = (vinfo[pos] << 32) | vinfo[pos + klen - 1];
    anchors.push_back(std::make_pair(ckmer_d, value));
  }
}

// TODO : we could improve this using some sort of caching
std::vector<rlpath_t> kdfs(const gbwtgraph::GBZ &gbz,
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

// Flag (set value to -1UL) for all anchors with repeated kmers (assuming
// anchors are sorted by kmer)
// XXX: could we set to 0 instead of -1UL?
size_t flag_repetitions(anchors_t &anchors) {
  size_t total = anchors.size();
  for (uint64_t i = 0; i < anchors.size();) {
    uint64_t kmer_d = anchors[i].first;
    ++i;
    if (anchors[i].first == kmer_d)
      --total;
    while (i < anchors.size() && anchors[i].first == kmer_d) {
      anchors[i - 1].second = -1UL;
      anchors[i].second = -1UL;
      ++i;
      --total;
    }
  }
  return total;
}

// flag repetitions checking also value of anchor, by comparing GFA identifiers
size_t flag_repetitions_wvalue(anchors_t &anchors, const Graph &graph) {
  size_t total = anchors.size();
  for (size_t i = 0; i < anchors.size();) {
    uint64_t kmer_d = anchors[i].first;
    size_t j = i + 1;
    while (j < anchors.size() && anchors[j].first == kmer_d)
      ++j;

    uint64_t value = anchors[i].second;
    uint32_t s1 = value >> 33, s2 = (uint32_t)value >> 1;
    // XXX: can we avoid checking GFA but directly using internal identifiers?
    std::string s1_gfa = graph.get_gfa_name(s1);
    std::string s2_gfa = graph.get_gfa_name(s2);

    size_t i2 = i + 1;
    for (; i2 < j; ++i2) {
      uint32_t n1 = anchors[i2].second >> 33,
               n2 = (uint32_t)anchors[i2].second >> 1;
      std::string n1_gfa = graph.get_gfa_name(n1);
      std::string n2_gfa = graph.get_gfa_name(n2);
      if (anchors[i2].second != value &&
          (s1_gfa.compare(n2_gfa) != 0 && s2_gfa.compare(n1_gfa) != 0))
        break;
    }
    if (i2 < j) {
      anchors[i].second = -1UL;
      --total;
    }
    for (i2 = i + 1; i2 < j; ++i2) {
      anchors[i2].second = -1UL;
      --total;
    }
    i = j;
  }
  return total;
}

int main_sketch(int argc, char *argv[]) {
  double rt = realtime();
  size_t klen = 31;               // kmer size
  std::string ref_path = "CHM13"; // reference paths
  // int nThreads = 4;               // number of threads

  int _c;
  while ((_c = getopt(argc, argv, "k:r:h")) != -1) {
    switch (_c) {
    case 'k':
      klen = std::stoi(optarg);
      break;
    // case '@':
    //   nThreads = std::stoi(optarg);
    //   break;
    case 'r':
      ref_path = optarg;
      break;
    case 'h':
      fprintf(stderr, "%s", SKETCH_USAGE_MESSAGE);
      return 0;
    }
  }

  if (argc - optind != 1) {
    fprintf(stderr, "%s", SKETCH_USAGE_MESSAGE);
    return 1;
  }
  std::string gbz_fn = argv[optind++];

  Graph graph(gbz_fn, ref_path);
  graph.load();
  graph.build_fl();
  // graph.load_fl();

  const gbwtgraph::GBZ &gbz = graph.get_gbz();
  const gbwt::FastLocate &fl = graph.get_fl();
  const gbwt::Metadata &metadata = graph.get_metadata();

  // STEP 1: get anchors from reference paths
  anchors_t anchors;
  gbwt::size_type sample_id = metadata.sample(ref_path);
  std::vector<gbwt::size_type> paths_identifiers =
      metadata.pathsForSample(sample_id);
  for (size_t pp = 0; pp < paths_identifiers.size(); ++pp) {
    gbwt::size_type path_id = paths_identifiers[pp];
    gbwt::size_type seq_id = gbwt::Path::encode(path_id, false);
    std::string sequence = ""; // full path sequence
    std::vector<gbwt::size_type>
        vinfo; // vertex identifier (w/ strand) for each position
    gbwt::edge_type curr = gbz.index.start(seq_id);
    while (curr.first != gbwt::ENDMARKER) {
      gbwtgraph::handle_t handle =
          gbwtgraph::GBWTGraph::node_to_handle(curr.first);
      // view_type: in-place view of the sequence: (start, length)
      gbwtgraph::view_type view = gbz.graph.get_sequence_view(handle);
      sequence.append(view.first, view.second);
      for (size_t x = 0; x < view.second; ++x)
        vinfo.push_back(curr.first);
      curr = gbz.index.LF(curr);
    }
    get_anchors(anchors, sequence, vinfo, klen, 0, sequence.size() - klen + 1);
  }
  fprintf(stderr,
          "[M::%s] Extracted %ld anchors from reference paths in %.3f secs\n",
          __func__, anchors.size(), realtime() - rt);

  rt = realtime();
  std::sort(anchors.begin(), anchors.end());
  fprintf(stderr, "[M::%s] Sorted anchors in %.3f secs\n", __func__,
          realtime() - rt);

  rt = realtime();
  size_t total_anchors = flag_repetitions(anchors);
  fprintf(stderr,
          "[M::%s] Flagged anchors in %.3f sec (resulting anchors: %ld)\n",
          __func__, realtime() - rt, total_anchors);

  rt = realtime();
  int shift = std::ceil(log2(total_anchors));
  sketch_t *sketch = sk_init((uint64_t)1 << shift, klen, 9); // XXX: hardcoded
  for (const auto &[kmer, value] : anchors) {
    if (value != -1UL)
      sk_insert(sketch, kmer, value);
  }
  fprintf(stderr,
          "[M::%s] Build reference sketch from %ld anchors in %.3f sec\n",
          __func__, sketch->n, realtime() - rt);

  // STEP 2: get all anchors by visiting each vertex (kDFS)
  int n_visited = 0;
  anchors.clear();
  rt = realtime();
  for (gbwtgraph::nid_t source_id = gbz.graph.min_node_id();
       source_id <= gbz.graph.max_node_id(); ++source_id) {

    /* note 1: segments longer than 1024 are split, so max_node_id >= #S
     * note 2: not all IDs in [min_node_id, max_node_id] are actually vertices
     * note 3: we need to check both strand for a node since outgoing edges
     *         depends on node strand
     **/

    for (int strand = 0; strand < 2; ++strand) {
      gbwtgraph::handle_t source_handle =
          gbz.graph.get_handle(source_id, strand);
      gbwt::node_type source = gbz.graph.handle_to_node(source_handle);

      gbwtgraph::view_type source_view =
          gbz.graph.get_sequence_view(source_handle);
      size_t source_length = source_view.second;
      std::string source_sequence(source_view.first, source_view.second);
      std::vector<gbwt::node_type> source_info(source_sequence.size(), source);
      anchors_t local_anchors;
      if (source_length >= klen)
        get_anchors(local_anchors, source_sequence, source_info, klen, 0,
                    source_length - klen + 1);

      // iterate over all paths starting from this vertex (that are
      // consistent with indexed haplotypes)
      std::vector<rlpath_t> paths = kdfs(gbz, fl, source, klen);
      for (const rlpath_t &path : paths) {
        std::string path_sequence;
        std::vector<gbwt::node_type> vinfo;

        // build path sequence and path info (vertex for each position)
        for (const gbwt::node_type &v : path.vertices) {
          gbwtgraph::handle_t h = gbwtgraph::GBWTGraph::node_to_handle(v);
          gbwtgraph::view_type view = gbz.graph.get_sequence_view(h);
          path_sequence.append(view.first, view.second);
          for (size_t i = 0; i < view.second; ++i)
            vinfo.push_back(v);
        }

        get_anchors(local_anchors, path_sequence, vinfo, klen,
                    source_length < klen ? 0 : source_length - klen + 1,
                    source_length);

        // std::cout << gbz.graph.get_segment_name(source_handle)
        //           << ((curr.first & 1) ? "-" : "+") << ","
        //           << path.vertices.size() << " ";
        // std::cout << graph.get_gfa_name(path.vertices.front() >> 1)
        //           << ((path.vertices.front() & 1) ? "-" : "+");
        // for (size_t v = 1; v < path.vertices.size(); ++v) {
        //   std::cout << "," << graph.get_gfa_name(path.vertices[v] >> 1)
        //             << ((path.vertices[v] & 1) ? "-" : "+");
        // }
        // std::cout << std::endl;
        // std::cout << path_sequence << std::endl;
      }

      for (const auto &[kmer, value] : local_anchors) {
        uint64_t hit = sk_get(sketch, kmer, 0);
        if (hit == -1UL) {
          // anchor is not in the reference
          anchors.push_back(std::make_pair(kmer, value));
        } else if (hit == 0) {
          // anchor was in the reference but we already found it repeated
          // somewhere else
        } else {
          if (hit == value) {
            // we just keep the reference anchor
          } else {
            uint32_t h1 = hit >> 33, h2 = (uint32_t)hit >> 1;
            uint32_t v1 = value >> 33, v2 = (uint32_t)value >> 1;
            // compare GFA identifier

            /* If a vertex is split (since longer than 1024), we may get two
             * different internal identifiers ("inverted") if we have same kmer
             * on different paths with different strand */

            // XXX: can we avoid checking GFA but directly using internal
            // identifiers?
            std::string h1_gfa = graph.get_gfa_name(h1);
            std::string h2_gfa = graph.get_gfa_name(h2);
            std::string v1_gfa = graph.get_gfa_name(v1);
            std::string v2_gfa = graph.get_gfa_name(v2);
            if (h1_gfa.compare(v2_gfa) != 0 && h2_gfa.compare(v1_gfa) != 0) {
              // anchor is repeated somewhere else
              // XXX: this may fail in some cases (?)
              sk_set(sketch, kmer, 0); // flag reference anchor as invalid
            }
          }
        }
      }

      ++n_visited;
      if (n_visited % 50000 == 0) {
        fprintf(stderr, "Visited %d vertices in %.3f sec\n", n_visited,
                realtime() - rt);
      }
    }
  }
  fprintf(stderr,
          "[M::%s] Extracted %ld anchors from other paths in %.3f secs\n",
          __func__, anchors.size(), realtime() - rt);

  rt = realtime();
  std::sort(anchors.begin(), anchors.end());
  fprintf(stderr, "[M::%s] Sorted other anchors in %.3f secs\n", __func__,
          realtime() - rt);

  rt = realtime();
  total_anchors = flag_repetitions_wvalue(anchors, graph);
  fprintf(
      stderr,
      "[M::%s] Flagged other anchors in %.3f sec (resulting anchors: %ld)\n",
      __func__, realtime() - rt, total_anchors);

  // Get final sketch by merging good reference and other paths anchors
  rt = realtime();
  size_t final_anchors = sketch->n + total_anchors;
  std::cerr << final_anchors << std::endl;
  shift = std::ceil(log2(final_anchors));
  sketch_t *final_sketch =
      sk_init((uint64_t)1 << shift, klen, 9); // XXX: hardcoded
  for (int64_t i = 0; i < sketch->n; ++i) {
    if (sketch->vls[i] != 0)
      sk_insert(final_sketch, sketch->sxs[i],
                (sketch->vls[i] << 1) | 1); // tag as 1 (last bit)
  }
  for (const auto &[kmer, value] : anchors) {
    if (value != -1UL)
      sk_insert(final_sketch, kmer, value << 1); // tag as 0 (last bit)
  }
  fprintf(stderr, "[M::%s] Build final sketch from %ld anchors in %.3f sec\n",
          __func__, final_sketch->n, realtime() - rt);

  sk_store(final_sketch, "-");
  sk_destroy(sketch);
  sk_destroy(final_sketch);

  return 0;
}
