#include <algorithm>
#include <filesystem>
#include <fstream>
#include <getopt.h>
#include <omp.h>
#include <set>
#include <string>
// #include <execution>

#include "anchor.hpp"
#include "graph.hpp"
#include "kmer.hpp"
#include "misc.hpp"
#include "sketch.hpp"
#include "usage.hpp"

#include "xxhash.h"

// Root/leaf path in the tree rooted at a given vertex
typedef struct {
  std::vector<gbwt::node_type> vertices;
  size_t l; // length in bp after first vertex (root)
} rlpath_t;

// Store anchors starting in positions [start, end) of sequence. Tag using vinfo
// (gbwt vertices w/ strand information)
void get_anchors(anchors_t &anchors, const std::string &sequence,
                 const std::vector<gbwt::node_type> &vinfo,
                 const std::map<gbwt::node_type, size_t> vlengths,
                 const size_t &klen, const size_t &start, const size_t &end,
                 bool is_ref, const uint64_t &density) {
  assert(sequence.size() == vinfo.size());
  char kmer[klen + 1];    // first kmer on sequence (plain)
  uint64_t kmer_d = 0;    // kmer
  uint64_t rckmer_d = 0;  // reverse and complemented kmer
  uint64_t ckmer_d = 0;   // canonical kmer
  uint8_t c;              // new character to append
  uint32_t pos;           // current position
  uint32_t pos2;          // position on "following" vertex
  gbwt::node_type v1, v2; // vertices
  XXH64_hash_t hash;

  // first kmer
  pos = start;
  strncpy(kmer, sequence.c_str() + pos, klen);
  kmer_d = k2d(kmer, klen);
  rckmer_d = rc(kmer_d, klen);
  ckmer_d = std::min(kmer_d, rckmer_d);

  // get positions
  v1 = vinfo[pos];
  v2 = vinfo[pos + klen - 1];
  if (v1 == v2) {
    pos2 = pos + klen - 1;
  } else {
    pos2 = 0;
    while (vinfo[pos + klen - 1 - pos2 - 1] == v2) {
      ++pos2;
    }
  }

  hash = XXH64(&ckmer_d, 8, 0); // xxseed = 0
  if (hash <= density) {
    if (kmer_d != ckmer_d) {
      // we do not have canonical on path sequence, so we have to reverse the
      // positions
      anchors.push_back(
          anchor_t{ckmer_d, (uint32_t)(v2 ^ 1), (uint32_t)(v1 ^ 1),
                   (uint32_t)vlengths.at(v2) - pos2 - 1,
                   (uint32_t)vlengths.at(v1) - pos - 1, is_ref, true});
    } else {
      anchors.push_back(anchor_t{ckmer_d, (uint32_t)v1, (uint32_t)v2, pos, pos2,
                                 is_ref, true});
    }
  }

  // all other kmers
  ++pos;
  for (; pos < end; ++pos) {
    // ACGT -> 0123
    c = to_int[(uint8_t)sequence[pos + klen - 1]] - 1; // A is 1 but we want 0
    kmer_d = lsappend(kmer_d, c, klen);
    rckmer_d = rsprepend(rckmer_d, reverse_char(c), klen);
    ckmer_d = std::min(kmer_d, rckmer_d);

    v1 = vinfo[pos];
    v2 = vinfo[pos + klen - 1];
    if (v1 == v2 || vinfo[pos + klen - 1 - 1] == v2) {
      // same vertex or same vertex as before
      ++pos2;
    } else {
      pos2 = 0;
    }

    hash = XXH64(&ckmer_d, 8, 0); // xxseed = 0
    if (hash <= density) {
      if (kmer_d != ckmer_d) {
        // we do not have canonical on path sequence, so we have to reverse the
        // positions
        anchors.push_back(anchor_t{
            ckmer_d, (uint32_t)(vinfo[pos + klen - 1] ^ 1),
            (uint32_t)(vinfo[pos] ^ 1),
            (uint32_t)vlengths.at(vinfo[pos + klen - 1]) - pos2 - 1,
            (uint32_t)vlengths.at(vinfo[pos]) - pos - 1, is_ref, true});
      } else {
        anchors.push_back(anchor_t{ckmer_d, (uint32_t)vinfo[pos],
                                   (uint32_t)vinfo[pos + klen - 1], pos, pos2,
                                   is_ref, true});
      }
    }
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

    bool extended = false;
    for (const gbwt::node_type &v : outs) {
      const gbwtgraph::handle_t &h = gbz.graph.node_to_handle(v);
      size_t l = gbz.graph.get_length(h);
      rlpath_t new_path = path;
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

// flag repetitions checking vertices/positions of anchor
size_t flag_repetitions(anchors_t &anchors, const Graph &graph) {
  size_t total = anchors.size();
  for (size_t i = 0; i < anchors.size();) {
    uint64_t kmer_d = anchors[i].kmer;
    size_t j = i + 1;
    while (j < anchors.size() && anchors[j].kmer == kmer_d)
      ++j;

    // We have same kmer in [i,j)

    bool is_ref = anchors[i].is_reference;
    bool has_both =
        false; // set to true if we also see v2>v1 with inverted strand

    uint32_t a_v1 = anchors[i].v1, a_v2 = anchors[i].v2;

    bool a_ir1 = a_v1 & 1, a_ir2 = a_v2 & 1;
    // remove strand bit
    a_v1 = a_v1 >> 1;
    a_v2 = a_v2 >> 1;
    size_t a_l1 = graph.get_vertex_len(a_v1), a_l2 = graph.get_vertex_len(a_v2);
    uint32_t a_p1 = anchors[i].pos1, a_p2 = anchors[i].pos2;
    // get offsets on + strand
    a_p1 = a_ir1 ? (a_l1 - a_p1 - 1) : a_p1;
    a_p2 = a_ir2 ? (a_l2 - a_p2 - 1) : a_p2;

    size_t i2 = i + 1;
    for (; i2 < j; ++i2) {
      // just flag as reference if one is reference
      is_ref |= anchors[i2].is_reference;

      uint32_t b_v1 = anchors[i2].v1, b_v2 = anchors[i2].v2;
      bool b_ir1 = b_v1 & 1, b_ir2 = b_v2 & 1;
      b_v1 = b_v1 >> 1;
      b_v2 = b_v2 >> 1;
      size_t b_l1 = graph.get_vertex_len(b_v1),
             b_l2 = graph.get_vertex_len(b_v2);
      uint32_t b_p1 = anchors[i2].pos1, b_p2 = anchors[i2].pos2;
      b_p1 = b_ir1 ? (b_l1 - b_p1 - 1) : b_p1;
      b_p2 = b_ir2 ? (b_l2 - b_p2 - 1) : b_p2;

      if ((a_v1 == b_v1 && a_v2 == b_v2) && (a_p1 == b_p1 && a_p2 == b_p2) &&
          (a_ir1 == b_ir1 && a_ir2 == b_ir2)) {
        // "same strand"
      } else if ((a_v1 == b_v2 && a_v2 == b_v1) &&
                 (a_p1 == b_p2 && a_p2 == b_p1) &&
                 (a_ir1 != b_ir2 && a_ir1 != b_ir2)) {
        // XXX: not sure if we should keep these or not, still they should be a
        // very tiny fraction On chr20 (8 samples), it is just 6/63696792
        // (9.4*10^-8)
        //
        // these are some sort of palindromes
        has_both = true;
      } else {
        // repetition
        break;
      }
    }

    // update reference and has_inverted bits just in case
    anchors[i].is_reference |= is_ref;
    anchors[i].has_both = has_both;

    if (i2 < j) {
      // tag first kmer as invalid, since we stopped earlier
      anchors[i].is_valid = false;
      --total;
    }
    for (i2 = i + 1; i2 < j; ++i2) {
      // tag all other kmers as invalid
      anchors[i2].is_valid = false;
      --total;
    }
    i = j;
  }
  return total;
}

int main_sketch(int argc, char *argv[]) {
  double rt = realtime();
  size_t klen = 31;               // kmer size
  uint64_t density = -1;          // max hash (based on density [0,1])
  std::string ref_path = "CHM13"; // reference paths
  int nThreads = 4;               // number of threads
  bool use_edges =
      true; // store kmers over edges (still used to flag duplicates)
  bool only_ref =
      false; // use only anchors from reference vertices longer than k

  int _c;
  while ((_c = getopt(argc, argv, "k:d:r:e1@:h")) != -1) {
    switch (_c) {
    case 'k':
      klen = std::stoi(optarg);
      break;
    case '1':
      only_ref = true;
      break;
    case 'd':
      density =
          atof(optarg) == 1.0 ? (uint64_t)-1 : (uint64_t)-1 * atof(optarg);
      break;
    case '@':
      nThreads = std::stoi(optarg);
      break;
    case 'r':
      ref_path = optarg;
      break;
    case 'e':
      use_edges = false;
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

  // Get all anchors by visiting each vertex (kDFS)
  std::vector<anchors_t> anchors_perthread(nThreads);
  // int n_visited = 0;
  rt = realtime();

  fprintf(stderr, "[M::%s] Vertex range: [%lld, %lld]\n", __func__,
          gbz.graph.min_node_id(), gbz.graph.max_node_id());

  size_t total_vertices = gbz.graph.max_node_id() - gbz.graph.min_node_id();

#pragma omp parallel for num_threads(nThreads) schedule(static, 1)
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
      gbwt::node_type source = gbz.graph.handle_to_node(source_handle);

      gbwtgraph::view_type source_view =
          gbz.graph.get_sequence_view(source_handle);
      size_t source_length = source_view.second;
      std::string source_sequence(source_view.first, source_view.second);
      std::vector<gbwt::node_type> source_info(source_sequence.size(), source);
      anchors_t local_anchors;
      bool source_isref = false;

      // iterate over all paths starting from this vertex (that are
      // consistent with indexed haplotypes)
      std::vector<rlpath_t> paths = kdfs(gbz, fl, source, klen);

      for (const rlpath_t &path : paths) {
        gbwt::size_type first;
        gbwt::SearchState state =
            fl.find(path.vertices.begin(), path.vertices.end(), first);
        assert(!state.empty());

        std::vector<gbwt::size_type> p_offsets = fl.locate(state, first);
        bool on_plus = false;
        bool is_ref = false;

        for (const gbwt::size_type &p : p_offsets) {
          // keep only paths on + strand
          if (gbwt::Path::is_reverse(p))
            continue;
          on_plus = true;
          if (metadata.fullPath(gbwt::Path::id(p))
                  .sample_name.compare(ref_path) == 0)
            is_ref = true;
        }

        if (!on_plus)
          continue;
        // we have at least one path on + strand

        if (is_ref)
          source_isref = true;

        std::string path_sequence;
        std::vector<gbwt::node_type> vinfo;
        std::map<gbwt::node_type, size_t> vlengths;

        // build path sequence and path info (vertex for each position)
        for (const gbwt::node_type &v : path.vertices) {
          gbwtgraph::handle_t h = gbwtgraph::GBWTGraph::node_to_handle(v);
          vlengths[v] = gbz.graph.get_length(h);
          gbwtgraph::view_type view = gbz.graph.get_sequence_view(h);
          path_sequence.append(view.first, view.second);
          for (size_t i = 0; i < view.second; ++i)
            vinfo.push_back(v);
        }

        if (path_sequence.size() < klen)
          continue;

        if (!only_ref) {
          get_anchors(local_anchors, path_sequence, vinfo, vlengths, klen,
                      source_length < klen ? 0 : source_length - klen + 1,
                      path_sequence.size() - source_length >= klen - 1
                          ? source_length
                          : path_sequence.size() - klen + 1,
                      is_ref, density);
        }
      }

      // add anchors from source vertex (tagging as reference if vertex is on
      // reference path)
      if (source_length >= klen && strand == 0) {
        if (!only_ref || source_isref) {
          std::map<gbwt::node_type, size_t> vlengths;
          vlengths[source] = source_length;
          get_anchors(local_anchors, source_sequence, source_info, vlengths,
                      klen, 0, source_length - klen + 1, source_isref, density);
        }
      }

      for (const anchor_t &a : local_anchors) {
        anchors_perthread[omp_get_thread_num()].push_back(a);
      }
    }

    // TODO: removed due to multithreading
    // ++n_visited;
    // if (n_visited % 500000 == 0) {
    //   fprintf(stderr, "Visited %d/%ld vertices in %.3f sec\n", n_visited,
    //           total_vertices, realtime() - rt);
    // }
  }

  size_t total_anchors = 0;
  std::vector<size_t> cumulative_counts(anchors_perthread.size());
  for (size_t aa = 0; aa < anchors_perthread.size(); ++aa) {
    cumulative_counts[aa] = total_anchors;
    total_anchors += anchors_perthread[aa].size();
  }

  fprintf(stderr,
          "[M::%s] Extracted %ld anchors from %ld vertices in %.3f secs\n",
          __func__, total_anchors, total_vertices, realtime() - rt);

  // TODO: we might sort, flag and shift anchors in parallel here so that we
  // (theoretically) need to sort less anchors later. This will require to get
  // new cumulative counts and new total_anchors before "inserting"

  anchors_t anchors(total_anchors);
  anchors.reserve(total_anchors);

  rt = realtime();
#pragma omp parallel for num_threads(nThreads) schedule(static, 1)
  for (size_t aa = 0; aa < anchors_perthread.size(); ++aa) {
    // std::sort(
    //     anchors_perthread[aa].begin(), anchors_perthread[aa].end(),
    //     [](const anchor_t &x, const anchor_t &y) { return x.kmer < y.kmer;
    //     });
    std::move(anchors_perthread[aa].begin(), anchors_perthread[aa].end(),
              anchors.begin() + cumulative_counts[aa]);
    anchors_perthread[aa].clear();
  }
  fprintf(stderr, "[M::%s] Sorted and linearized anchors in %.3f secs\n",
          __func__, realtime() - rt);

  rt = realtime();
  std::sort(
      // std::execution::par,
      anchors.begin(), anchors.end(),
      [](const anchor_t &x, const anchor_t &y) { return x.kmer < y.kmer; });
  fprintf(stderr, "[M::%s] Sorted all anchors in %.3f secs\n", __func__,
          realtime() - rt);

  rt = realtime();
  total_anchors = flag_repetitions(anchors, graph);
  fprintf(
      stderr,
      "[M::%s] Flagged other anchors in %.3f sec (resulting anchors: %ld)\n",
      __func__, realtime() - rt, total_anchors);

  // Build sketch
  rt = realtime();
  int shift = std::ceil(log2(total_anchors));
  sketch_t *sketch =
      sk_init((uint64_t)1 << shift, klen, 9); // XXX: hardcoded prefix size to 9
  if (sketch == nullptr) {
    fprintf(stderr, "[E::%s] Cannot allocate sketch. Halting..\n", __func__);
    return -1;
  }

  for (const anchor_t &a : anchors) {
    if (a.is_valid) {
      if (use_edges || a.v1 == a.v2)
        // if (a.v1 != a.v2)
        sk_insert(sketch, a.kmer, a.v1, a.v2, a.pos1, a.pos2, a.has_both,
                  a.is_reference);
    }
  }
  fprintf(stderr, "[M::%s] Built final sketch from %ld anchors in %.3f sec\n",
          __func__, sketch->n, realtime() - rt);

  sk_store(sketch, "-");
  sk_destroy(sketch);

  return 0;
}
