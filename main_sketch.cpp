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

typedef std::vector<std::pair<uint64_t, uint64_t>> kmers_t;

typedef struct path {
  std::vector<gbwt::node_type> vertices;
  size_t l;
} subpath_t;

// Store anchors starting in positions [start, end) of sequence. Tag using info
void count_kmers(kmers_t &kmers, const std::string &sequence,
                 const std::vector<gbwt::node_type> &info, int klen,
                 size_t start, size_t end) {
  assert(sequence.size() == info.size());
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
  value = (info[pos] << 32) | info[pos + klen - 1];
  kmers.push_back(std::make_pair(ckmer_d, value));

  // all other kmers
  ++pos;
  for (; pos < end; ++pos) {
    c = to_int[(uint8_t)sequence[pos + klen - 1]] -
        1; // A is 1 but it should be 0
    kmer_d = lsappend(kmer_d, c, klen);
    rckmer_d = rsprepend(rckmer_d, reverse_char(c), klen);
    ckmer_d = std::min(kmer_d, rckmer_d);
    value = (info[pos] << 32) | info[pos + klen - 1];
    kmers.push_back(std::make_pair(ckmer_d, value));
  }
}

// TODO : we could improve this using some sort of caching

std::vector<subpath_t> kdfs(const gbwtgraph::GBZ &gbz,
                            const gbwt::FastLocate &fl, gbwt::node_type source,
                            int klen) {
  std::vector<subpath_t> paths; // fully extended path
  std::queue<subpath_t> queue;  // paths we still need to "extend"
  queue.push({{source}, 0});
  while (!queue.empty()) {
    subpath_t path = queue.front();
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
      subpath_t new_path = path;
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

// Flag (set value to -1UL) for all repeated kmers
size_t flag_repetitions(kmers_t &kmers) {
  size_t total_kmers = kmers.size();
  for (uint64_t i = 0; i < kmers.size();) {
    uint64_t kmer_d = kmers[i].first;
    ++i;
    if (kmers[i].first == kmer_d)
      --total_kmers;
    while (i < kmers.size() && kmers[i].first == kmer_d) {
      kmers[i - 1].second = -1UL;
      kmers[i].second = -1UL;
      ++i;
      --total_kmers;
    }
  }
  return total_kmers;
}

int main_sketch(int argc, char *argv[]) {
  double rt = realtime();
  uint klen = 31;                 // kmer size
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

  // std::cerr << gbz.graph.get_node_count() << " [" << gbz.graph.min_node_id()
  //           << ", " << gbz.graph.max_node_id() << "]" << std::endl;
  /* note 1: segments longer than 1024 are split, so max_node_id >= #S
   * note 2: not all IDs in [min_node_id, max_node_id] are actually vertices
   **/

  // STEP 1: get anchors from reference paths
  kmers_t kmers;
  gbwt::size_type sample_id = metadata.sample(ref_path);
  std::vector<gbwt::size_type> paths_identifiers =
      metadata.pathsForSample(sample_id);
  for (size_t pp = 0; pp < paths_identifiers.size(); ++pp) {
    gbwt::size_type path_id = paths_identifiers[pp];
    gbwt::size_type seq_id = gbwt::Path::encode(path_id, false);
    std::string sequence = ""; // full path sequence
    std::vector<gbwt::size_type>
        info; // vertex identifier (w/ strand) for each position
    gbwt::edge_type curr = gbz.index.start(seq_id);
    while (curr.first != gbwt::ENDMARKER) {
      gbwtgraph::handle_t handle =
          gbwtgraph::GBWTGraph::node_to_handle(curr.first);
      // view_type: in-place view of the sequence: (start, length)
      gbwtgraph::view_type view = gbz.graph.get_sequence_view(handle);
      sequence.append(view.first, view.second);
      for (size_t x = 0; x < view.second; ++x)
        info.push_back(curr.first);
      curr = gbz.index.LF(curr);
    }
    count_kmers(kmers, sequence, info, klen, 0, sequence.size() - klen + 1);
  }
  fprintf(stderr,
          "[M::%s] Extracted %ld kmers from reference paths in %.3f secs\n",
          __func__, kmers.size(), realtime() - rt);

  rt = realtime();
  std::sort(kmers.begin(), kmers.end());
  fprintf(stderr, "[M::%s] Sorted kmers in %.3f secs\n", __func__,
          realtime() - rt);

  rt = realtime();
  size_t total_kmers = flag_repetitions(kmers);
  fprintf(stderr, "[M::%s] Flagged kmers in %.3f sec (resulting kmers: %ld)\n",
          __func__, realtime() - rt, total_kmers);

  rt = realtime();
  int shift = std::ceil(log2(total_kmers));
  sketch_t *sketch = sk_init((uint64_t)1 << shift, klen, 9); // XXX: hardcoded
  for (const auto &[kmer, value] : kmers) {
    if (value != -1UL)
      sk_insert(sketch, kmer, value);
  }
  kmers.clear();
  fprintf(stderr, "[M::%s] Build reference sketch from %ld kmers in %.3f sec\n",
          __func__, sketch->n, realtime() - rt);

  // STEP 2: get all kmers by visiting each vertex (kDFS)
  int n_visited = 0;
  kmers_t other_kmers;
  rt = realtime();
  for (gbwtgraph::nid_t source_id = gbz.graph.min_node_id();
       source_id <= gbz.graph.max_node_id(); ++source_id) {
    for (int strand = 0; strand < 2; ++strand) {
      gbwtgraph::handle_t source_handle =
          gbz.graph.get_handle(source_id, strand);
      gbwt::node_type source = gbz.graph.handle_to_node(source_handle);

      gbwtgraph::view_type source_view =
          gbz.graph.get_sequence_view(source_handle);
      // view_type: in-place view of the sequence: (start, length)
      size_t source_length = source_view.second;
      std::string source_sequence(source_view.first, source_view.second);
      std::vector<gbwt::node_type> source_info(source_sequence.size(), source);
      kmers_t local_kmers;
      if (source_length >= klen)
        count_kmers(local_kmers, source_sequence, source_info, klen, 0,
                    source_length - klen + 1);

      // iterate over all paths starting from this vertex (that are
      // consistent with indexed haplotypes)
      std::vector<subpath_t> paths = kdfs(gbz, fl, source, klen);
      for (const subpath_t &path : paths) {
        std::string path_sequence;
        std::vector<gbwt::node_type> info;

        // build path sequence and path info (vertex for each position)
        for (const gbwt::node_type &v : path.vertices) {
          gbwtgraph::handle_t h = gbwtgraph::GBWTGraph::node_to_handle(v);
          gbwtgraph::view_type view = gbz.graph.get_sequence_view(h);
          path_sequence.append(view.first, view.second);
          for (size_t i = 0; i < view.second; ++i)
            info.push_back(v);
        }

        count_kmers(local_kmers, path_sequence, info, klen,
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

      for (const auto &[kmer, value] : local_kmers) {
        uint64_t hit = sk_get(sketch, kmer);
        if (hit == -1UL) {
          // anchor is not in the reference
          other_kmers.push_back(std::make_pair(kmer, value));
          // other_kmers[kmer] = value;
        } else if (hit == 0) {
        } else {
          if (hit != value) {
            uint32_t h1 = hit >> 33, h2 = (uint32_t)hit >> 1;
            uint32_t v1 = value >> 33, v2 = (uint32_t)value >> 1;

            std::string h1_gfa = graph.get_gfa_name(h1);
            std::string h2_gfa = graph.get_gfa_name(h2);
            std::string v1_gfa = graph.get_gfa_name(v1);
            std::string v2_gfa = graph.get_gfa_name(v2);
            if (h1_gfa.compare(v2_gfa) != 0 && h2_gfa.compare(v1_gfa) != 0) {
              // if (h1 != v1 && h2 != v2) {
              // char ks[klen];
              // d2s(kmer, klen, ks);
              // std::cerr << kmer << " " << ks << std::endl;
              // std::cerr << h1 << "/" << h2 << std::endl;
              // std::cerr << h1_gfa << "/" << h2_gfa << std::endl;
              // std::cerr << v1 << "/" << v2 << std::endl;
              // std::cerr << v1_gfa << "/" << v2_gfa << std::endl;
              // std::cerr << "=== === ===" << std::endl;

              // anchor is repeated somewhere else
              // XXX: this may fail in some cases
              sk_set(sketch, kmer, 0);
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
  fprintf(stderr, "[M::%s] Extracted %ld kmers from other paths in %.3f secs\n",
          __func__, other_kmers.size(), realtime() - rt);

  rt = realtime();
  std::sort(other_kmers.begin(), other_kmers.end());
  fprintf(stderr, "[M::%s] Sorted other kmers in %.3f secs\n", __func__,
          realtime() - rt);

  rt = realtime();
  size_t total_other_kmers = flag_repetitions(other_kmers);
  fprintf(stderr,
          "[M::%s] Flagged other kmers in %.3f sec (resulting kmers: %ld)\n",
          __func__, realtime() - rt, total_other_kmers);

  // Get final sketch by merging good reference and other paths anchors
  rt = realtime();
  size_t final_kmers = sketch->n + total_other_kmers;
  std::cerr << final_kmers << std::endl;
  shift = std::ceil(log2(final_kmers));
  sketch_t *final_sketch =
      sk_init((uint64_t)1 << shift, klen, 9); // XXX: hardcoded
  for (int64_t i = 0; i < sketch->n; ++i) {
    if (sketch->vls[i] != 0)
      sk_insert(final_sketch, sketch->sxs[i],
                (sketch->vls[i] << 1) | 1); // tag as 1 (last bit)
  }
  for (const auto &[kmer, value] : other_kmers) {
    if (value != -1UL)
      sk_insert(final_sketch, kmer, value << 1); // tag as 0 (last bit)
  }
  fprintf(stderr, "[M::%s] Build final sketch from %ld kmers in %.3f sec\n",
          __func__, final_sketch->n, realtime() - rt);

  sk_store(final_sketch, "-");
  sk_destroy(sketch);
  sk_destroy(final_sketch);

  return 0;
}
