// #include <algorithm>
#include <assert.h>
#include <getopt.h>
#include <omp.h>
#include <vector>

extern "C" {
#include "fm-index.h"
}

#include "anchor.hpp"
#include "graph.hpp"
#include "misc.hpp"
#include "reads.hpp"
#include "sfs.hpp"
#include "sketch.hpp"
#include "usage.hpp"

// Compute SFS strings from P and store them into solutions
std::vector<sfs_t> ping_pong_search(const rb3_fmi_t *index, uint8_t *P, int l,
                                    char *name) {
  std::vector<sfs_t> out;
  rb3_sai_t ik;
  int begin = l - 1;
  while (begin >= 0) {
    // Backward search. Stop at first mismatch.
    // int bmatches = 0;
    rb3_fmd_set_intv(index, P[begin], &ik);
    while (ik.size != 0 && begin > 0) {
      --begin;
      // ++bmatches;
      rb3_sai_t ok[RB3_ASIZE]; // output SA intervals (one for each symbol
                               // between 0 and 5)
      rb3_fmd_extend(index, &ik, ok, 1);
      ik = ok[P[begin]];
    }
    // last checked char (i.e., first of the query) was a match
    if (begin == 0 && ik.size != 0) {
      break;
    }

    // Forward search
    int end = begin;
    // int fmatches = 0;
    rb3_fmd_set_intv(index, P[end], &ik);
    while (ik.size != 0) {
      ++end;
      // ++fmatches;
      rb3_sai_t ok[RB3_ASIZE];
      rb3_fmd_extend(index, &ik, ok, 0);
      ik = ok[P[end] >= 1 && P[end] <= 4 ? 5 - P[end] : P[end]];
    }

    out.push_back({name, begin, end - begin + 1});

    if (begin == 0)
      break;
    // if (config->overlap == 0) // Relaxed
    //   begin -= 1;
    // else
    begin = end - 1;
  }
  return out;
}

// Merge specifics strings that are too close (d-bp apart) on read
void assemble(std::vector<sfs_t> &S, int d = 0) {
  // Reverse the vector
  std::reverse(S.begin(), S.end());

  size_t i = 0;
  while (i < S.size()) {
    size_t j;
    for (j = i + 1; j <= S.size(); ++j) {
      if (j == S.size() || S[j - 1].s + S[j - 1].l <= S[j].s - d) {
        // non-overlapping: update first, clean others
        S[i].l = S[j - 1].s + S[j - 1].l - S[i].s;
        for (size_t j2 = i + 1; j2 < j; ++j2)
          S[j2].l = 0;
        break;
      }
    }
    i = j;
  }

  // Remove gaps by shifting left and resize
  int new_n = 0;
  i = 0;
  while (i < S.size()) {
    if (S[i].l > 0) {
      S[new_n].s = S[i].s;
      S[new_n].l = S[i].l;
      ++new_n;
    }
    ++i;
  }
  S.resize(new_n);
}

// 1234-read to 0123-kmer counts
std::map<uint64_t, int> count_kmers(uint8_t *read, int readl, int klen) {
  uint8_t c; // current char
  std::map<uint64_t, int> kcounts;

  uint64_t kmer_d = 0;
  for (uint8_t i = 0; i < klen; ++i) {
    c = read[i] - 1;
    kmer_d = (kmer_d << 2) | (c < 4 ? c : rand() % 4);
  }

  uint64_t rckmer_d = rc(kmer_d, klen);
  uint64_t ckmer_d = std::min(kmer_d, rckmer_d);
  ++kcounts[ckmer_d];
  for (int p = klen; p < readl; ++p) {
    c = read[p] - 1;
    kmer_d = lsappend(kmer_d, c, klen);
    rckmer_d = rsprepend(rckmer_d, reverse_char(c), klen);
    ckmer_d = std::min(kmer_d, rckmer_d);
    ++kcounts[ckmer_d];
  }

  return kcounts;
}

// Get paths containing an anchor (since we store only start/end vertices, we
// need to subset paths over bubbles)
void set_paths(const Graph &graph, anchor_t &anchor, int klen) {
  uint32_t v1 = anchor.v1;
  uint32_t v2 = anchor.v2;

  // std::cerr << anchor.kmer << " " << anchor.is_canonical << std::endl;
  // std::cerr << graph.get_gfa_name(v1 >> 1) << "/" << (v1 & 1) << ">"
  //           << graph.get_gfa_name(v2 >> 1) << "/" << (v2 & 1) << std::endl;
  // std::cerr << std::endl;

  uint8_t strand_flags = 3; // take paths from both strands
  if (v1 != v2 && !anchor.has_both)
    strand_flags = 1; // take paths on + strand only
  // strand_flags = 2; // take paths on - strand only

  std::vector<path_t> paths = graph.get_paths(v1, v2, strand_flags, false);

  // get only paths containing the anchor kmer
  // since paths can have vertices on both strand, we need to check for both
  // "versions" of the kmer
  char skmer[klen];
  char skmer_rc[klen];
  if (anchor.is_canonical) {
    d2s(anchor.kmer, klen, skmer);
    d2s(rc(anchor.kmer, klen), klen, skmer_rc);
  } else {
    d2s(rc(anchor.kmer, klen), klen, skmer);
    d2s(anchor.kmer, klen, skmer_rc);
  }

  // std::cerr << graph.get_gfa_name(anchor.v1 >> 1) << ">"
  //           << graph.get_gfa_name(anchor.v2 >> 1) << "/" << anchor.is_reverse
  //           << std::endl;

  for (size_t i = 0; i < paths.size(); ++i) {
    const path_t &path = paths[i];
    // std::cerr << graph.get_path_contig(path.id >> 1) << " " << (path.id & 1)
    //           << std::endl;
    size_t p1 = path.sequence.find(skmer);
    size_t p0 = path.sequence.find(skmer_rc);
    if (p1 != std::string::npos || p0 != std::string::npos) {
      assert(p1 != p0);
      // std::cerr << "O " << graph.get_path_sample(path.id >> 1) << " "
      //           << graph.get_path_contig(path.id >> 1) << std::endl;
      anchor.paths[path.id] = {(uint32_t)path.offset1, (uint32_t)path.offset2,
                               path.is_reference};
    } else {
      // std::cerr << "X " << graph.get_path_sample(path.id >> 1) << " "
      //           << graph.get_path_contig(path.id >> 1) << std::endl;
    }
  }
  // std::cerr << std::endl << std::endl;
}

// Compute distance in bp between v1 and v2 on path with identifier pid
// XXX: can we get the "offset" directly when building the path so we do not
// need to decompressSA here again?
size_t compute_distance_bp(const Graph &graph, const gbwt::node_type &v1,
                           const gbwt::node_type &v2,
                           const gbwt::size_type &pid) {
  const gbwt::FastLocate &fl = graph.get_fl();

  std::vector<gbwt::size_type> intervals1 = fl.decompressSA(v1);

  for (size_t i = 0; i < intervals1.size(); ++i) {
    gbwt::FastLocate::size_type int1 = intervals1[i];
    if (fl.seqId(int1) != pid)
      continue;
    gbwt::edge_type position = std::make_pair(v1, i);

    size_t l = 0;
    // skip v1
    position = graph.get_gbz().index.LF(position);
    while (position.first != v2) {
      gbwtgraph::handle_t handle =
          gbwtgraph::GBWTGraph::node_to_handle(position.first);
      l += graph.get_gbz().graph.get_length(handle);
      position = graph.get_gbz().index.LF(position);
    }
    return l;
  }
  return -1UL;
}

// Chaining heuristic
std::pair<anchors_t, bool> chain(const anchors_t &anchors, const Graph &graph,
                                 const gbwt::size_type &pid, int klen) {
  anchors_t best_chain;
  bool best_strand;
  anchors_t chain;
  bool strand;

  for (size_t x = 0; x < anchors.size(); ++x) {
    chain.push_back(anchors.at(x));
    for (size_t y = x + 1; y < anchors.size(); ++y) {
      const anchor_t &yy = anchors.at(y);

      int curr_strand;

      int qd = yy.qp - chain.back().qp;
      assert(qd > 0);

      int rd;
      /*
        chain.back() > yy
        v1---v2        v1---v2

        or

        yy           > chain.back()
        v2---v1        v2---v1
       */
      if (chain.back().v2 == yy.v1) {
        // "internal" ends are on the same vertex
        curr_strand = chain.back().pos1 < yy.pos1;
        if (curr_strand)
          rd = yy.pos1 - chain.back().pos2 + klen - 1;
        else
          rd = chain.back().pos2 - yy.pos1 + klen - 1;
        // qd = std::abs(qd);
      } else {
        curr_strand = chain.back().paths.at(pid).offset1 >
                      chain.back().paths.at(pid).offset2;
        if (chain.back().v1 == yy.v1 && chain.back().v2 == yy.v2) {
          // overlapping anchors over same edge (same pair of vertices)
          if (curr_strand)
            rd = yy.pos1 - chain.back().pos1;
          else
            rd = chain.back().pos1 - yy.pos1;
        } else {
          // anchors on two different vertices
          if (curr_strand) {
            int vl = graph.get_vertex_len(chain.back().v2 >> 1);
            rd = compute_distance_bp(graph, chain.back().v2, yy.v1, pid) +
                 (vl - chain.back().pos1) + yy.pos1;
          } else {
            int vl = graph.get_vertex_len(yy.v2 >> 1);
            rd = compute_distance_bp(graph, yy.v2, chain.back().v1, pid) +
                 (vl - yy.pos1) + chain.back().pos1;
          }
        }
      }
      assert(rd > 0);

      if (std::max(rd, qd) == 0 || std::min(rd, qd) / (float)std::max(rd, qd) >
                                       0.9) { // FIXME: hardcoded
        if (chain.size() == 1) {
          strand = curr_strand;
        }
        if (strand == curr_strand)
          chain.push_back(yy);
      }
    }
    if (chain.size() / (float)anchors.size() > 0.5) {
      best_chain = chain;
      best_strand = strand;
      break;
    } else {
      if (chain.size() > best_chain.size()) {
        best_chain = chain;
        best_strand = strand;
      }
    }
    chain.clear();
  }
  return std::make_pair(best_chain, best_strand);
}

typedef struct {
  anchor_t anchor1;
  anchor_t anchor2;
  int chain_size;
  int total_anchors;
  bool strand;
} anchoring_t;

std::map<gbwt::size_type, anchoring_t>
chaining(const Graph &graph, const anchors_t &anchors, int klen) {
  std::map<gbwt::size_type, anchoring_t> result;

  // Get how many anchors we have on "each" path
  std::map<uint32_t, std::vector<uint32_t>> pcounts;
  for (size_t a = 0; a < anchors.size(); ++a) {
    const anchor_t &aa = anchors[a];

    // std::cerr << a << ": " << aa.kmer << " " << aa.qp << " "
    //           << graph.get_gfa_name(aa.v1 >> 1) << ">"
    //           << graph.get_gfa_name(aa.v2 >> 1) << " " << aa.v1 << ">" <<
    //           aa.v2
    //           << "/" << aa.has_both << ":" << aa.pos1 << "-" << aa.pos2
    //           << " >> ";
    for (const auto &[pid, p] : aa.paths) {
      // std::cerr << (pid >> 1) << "/" << (pid & 1) << ":"
      //           << graph.get_path_contig(pid >> 1) << " ";
      pcounts[pid].push_back(a);
    }
    // std::cerr << std::endl;
  }
  //
  // We may want to merge paths with same anchors, but then we need to process
  // them "all" anyway since two anchors can have different distance on
  // different path
  //
  // Sort paths by decreasing number of anchors
  std::vector<std::pair<uint32_t, uint32_t>> sorted_pcounts;
  for (const auto &[k, v] : pcounts)
    sorted_pcounts.push_back(std::make_pair(v.size(), k));
  std::sort(sorted_pcounts.begin(), sorted_pcounts.end(),
            [](const auto &a, const auto &b) { return a > b; });
  //
  // Iterate over paths and chain
  for (const auto &[na, p] : sorted_pcounts) {
    // p is path id w/ strand bit
    anchors_t path_anchors;
    // std::cerr << p << " " << graph.get_path_contig(p >> 1) << " " << (p & 1)
    //           << " " << na << " :\n";
    for (const size_t aa : pcounts[p]) {
      anchor_t a = anchors[aa];
      if (p & 1) {
        //
        // we do not need to reverse the anchor since path strand can be - if
        // vertex was + but on path is < in the case of single vertex anchors,
        // offset is always wrt + strand
        //
        // std::swap(a.pos1, a.pos2);
        // std::swap(a.v1, a.v2);
        // a.pos1 = graph.get_vertex_len(a.v1 >> 1) - a.pos1 - 1;
        // a.pos2 = graph.get_vertex_len(a.v2 >> 1) - a.pos2 - 1;
      }
      // assert(a.pos1 < a.pos2);
      // std::cerr << "\t" << a.qp << "\t" << a.v1 << ":" << a.pos1 << " / "
      //           << a.v2 << ":" << a.pos2 << "\t"
      //           << graph.get_gfa_name(a.v1 >> 1) << ">"
      //           << graph.get_gfa_name(a.v2 >> 1) << "\t" << a.kmer << "\t";

      // std::cerr << p << ":" << a.paths[p].offset1 << "-" <<
      // a.paths[p].offset2
      //           << " ";
      // std::cerr << std::endl;

      path_anchors.push_back(a);
    }
    // std::cerr << std::endl;

    // chain
    std::pair<anchors_t, bool> bc = chain(path_anchors, graph, p, klen);
    result[p] = {bc.first.front(), bc.first.back(), (int)bc.first.size(),
                 (int)path_anchors.size(), bc.second};
    // std::cerr << bc.first.size() << " " << (bc.second ? "+" : "-") <<
    // std::endl;
  }
  return result;
}

// Anchor specific strings on graph using graph sketch
void anchor(const Graph &graph, sketch_t *sketch, std::vector<sfs_t> &sfs,
            uint8_t *read, int readl, std::map<uint64_t, int> kcounts, int klen,
            size_t NA, bool reference_only) {
  int beg, end;
  uint8_t *kmer = (uint8_t *)malloc(klen);
  uint64_t kmer_d;       // kmer
  uint64_t rckmer_d = 0; // reverse and complemented kmer
  uint64_t ckmer_d = 0;  // canonical kmer
  int c;                 // current char
  hit_t hit;             // hit from sketch
  int count;

  //   for (size_t i = 0; i < readl; ++i)
  //     std::cerr << (int)read[i];
  //   std::cerr << std::endl;

  for (uint sidx = 0; sidx < sfs.size(); ++sidx) {
    sfs_t &s = sfs[sidx];

    // Finding anchors on the left
    beg = s.s - klen + 1;
    beg = beg < 0 ? 0 : beg;
    memcpy(kmer, read + beg, klen);
    kmer[klen] = '\0';

    kmer_d = k2d((char *)kmer, klen);
    rckmer_d = rc(kmer_d, klen);
    ckmer_d = std::min(kmer_d, rckmer_d);

    anchors_t sanchors;
    while (sanchors.size() < NA) {
      hit = sk_get(sketch, ckmer_d, reference_only);
      count = kcounts[ckmer_d];
      if (hit.value != -1UL && count == 1) {
        anchor_t a;
        a.kmer = ckmer_d;
        a.v1 = hit.value >> 32;
        a.v2 = (uint32_t)hit.value;
        a.pos1 = (hit.info >> 17) & 0x7FFF;
        a.pos2 = (hit.info >> 2) & 0x7FFF;
        a.qp = beg;
        a.has_both = (hit.info >> 1) & 1;
        a.is_reference = hit.info & 1;
        a.is_canonical = ckmer_d == kmer_d;

        sanchors.push_back(a);
        set_paths(graph, sanchors.back(), klen);
      }
      --beg;
      if (beg < 0)
        break;

      c = read[beg] - 1;
      kmer_d = rsprepend(kmer_d, c, klen);
      rckmer_d = lsappend(rckmer_d, reverse_char(c), klen);
      ckmer_d = std::min(kmer_d, rckmer_d);
    }

    // Finding anchors on the right
    end = s.s + s.l - 1;
    end = end > readl - klen ? readl - klen : end;
    memcpy(kmer, read + end, klen);

    kmer_d = k2d((char *)kmer, klen);
    rckmer_d = rc(kmer_d, klen);
    ckmer_d = std::min(kmer_d, rckmer_d);
    anchors_t eanchors;
    while (eanchors.size() < NA) {
      hit = sk_get(sketch, ckmer_d, reference_only);
      count = kcounts[ckmer_d];
      memcpy(kmer, read + end, klen);
      if (hit.value != -1UL && count == 1) {
        anchor_t a;
        a.kmer = ckmer_d;
        a.v1 = hit.value >> 32;
        a.v2 = (uint32_t)hit.value;
        a.pos1 = (hit.info >> 17) & 0x7FFF;
        a.pos2 = (hit.info >> 2) & 0x7FFF;
        a.qp = end;
        a.has_both = (hit.info >> 1) & 1;
        a.is_reference = hit.info & 1;
        a.is_canonical = ckmer_d == kmer_d;

        eanchors.push_back(a);
        set_paths(graph, eanchors.back(), klen);
      }
      ++end;
      if (end == readl - klen + 1)
        break;

      c = read[end + klen - 1] - 1;
      kmer_d = lsappend(kmer_d, c, klen);
      rckmer_d = rsprepend(rckmer_d, reverse_char(c), klen);
      ckmer_d = std::min(kmer_d, rckmer_d);
    }

    // std::cerr << sanchors.size() << " " << eanchors.size() << std::endl;
    if (sanchors.size() == 0 || eanchors.size() == 0) {
      s.flag |= 1; // flag as invalid
      continue;
    }

    std::reverse(sanchors.begin(), sanchors.end());
    std::map<gbwt::size_type, anchoring_t> schains =
        chaining(graph, sanchors, klen);
    std::map<gbwt::size_type, anchoring_t> echains =
        chaining(graph, eanchors, klen);

    std::map<std::pair<anchor_t, anchor_t>, std::vector<gbwt::size_type>>
        anchoring;
    for (const auto &[p, sanchoring] : schains) {
      auto x = echains.find(p);
      if (x == echains.end())
        continue;
      if (sanchoring.strand != x->second.strand)
        continue;
      // std::cerr << p << " " << graph.get_path_contig(p >> 1) << std::endl;
      // std::cerr << sanchoring.anchor2.v1 << " " << x->second.anchor1.v2
      //           << std::endl;
      anchoring[std::make_pair(sanchoring.anchor2, x->second.anchor1)]
          .push_back(((p >> 1) << 1) | sanchoring.strand);
    }
    if (anchoring.empty()) {
      s.flag |= 2; // flag as invalid
      continue;
    }

    std::pair<anchor_t, anchor_t> best_anchors;
    std::vector<uint64_t> best_paths;
    for (const auto &[as, ps] : anchoring) {
      if (ps.size() > best_paths.size()) {
        best_paths = ps;
        best_anchors = as;
      }
    }

    // std::cerr << "======= ";
    // for (const auto &p : best_paths)
    //   std::cerr << graph.get_path_contig(p >> 1) << "/" << (p & 1 ? "+" :
    //   "-")
    //             << " ";
    // std::cerr << "=======" << std::endl;
    // std::cerr << best_anchors.first.qp << " " << best_anchors.first.v1 << ":"
    //           << best_anchors.first.pos1 << " > " << best_anchors.first.v2
    //           << ":" << best_anchors.first.pos2 << std::endl;
    // std::cerr << best_anchors.second.qp << " " << best_anchors.second.v1 <<
    // ":"
    //           << best_anchors.second.pos1 << " > " << best_anchors.second.v2
    //           << ":" << best_anchors.second.pos2 << std::endl;

    assert(best_anchors.first.qp < best_anchors.second.qp);

    // Assigning the anchors
    s.s = best_anchors.first.qp;
    s.l = best_anchors.second.qp + klen - s.s;
    s.sv = best_anchors.first.v1;
    s.ev = best_anchors.second.v2;
    // XXX: not sure wrt which strand I have this positions
    s.soff = best_anchors.first.pos1;
    s.eoff = best_anchors.second.pos2;
    s.paths = best_paths;
    s.skmer = best_anchors.first.kmer;
    s.ekmer = best_anchors.second.kmer;
  }
  free(kmer);
}

// void remove_duplicates(std::vector<sfs_t> &S) {
//   // XXX: what about strings that are prefix/suffix or overlapping?
//   for (size_t i = 0; i < S.size(); ++i) {
//     if (S[i].flag != 0)
//       continue;
//     for (size_t j = i + 1; j < S.size(); ++j) {
//       if (S[j].flag != 0)
//         continue;
//       if (S[i].s == S[j].s && S[i].l == S[j].l)
//         S[j].flag = 3;
//     }
//   }
// }

void fill_sequence(sfs_t &s, uint8_t *read) {
  s.plain_seq = std::string((char *)read + s.s, s.l);
  for (char &c : s.plain_seq)
    c = "0ACGT0"[(uint8_t)c];
}

int main_sfs(int argc, char *argv[]) {
  double rt;
  int bsize = 10000; // batch size
  size_t NA = 20;    // number of kmers to check for anchoring
  int nth = 1;       // number of threads
  bool reference_only = false;
  bool anchoring = true;

  int _c;
  while ((_c = getopt(argc, argv, "a:b:nr@:h")) != -1) {
    switch (_c) {
    case 'a':
      NA = std::stoi(optarg);
      break;
    case 'b':
      bsize = std::stoi(optarg);
      break;
    case 'n':
      anchoring = false;
      break;
    case 'r':
      reference_only = true;
      break;
    case '@':
      nth = std::stoi(optarg);
      break;
    case 'h':
      fprintf(stderr, "%s", SEARCH_USAGE_MESSAGE);
      return 0;
    }
  }

  if (argc - optind != 4) {
    fprintf(stderr, "%s", SEARCH_USAGE_MESSAGE);
    return 1;
  }
  std::string gbz_fn = argv[optind++];
  std::string skt_fn = argv[optind++];
  std::string fmd_fn = argv[optind++];
  std::string fx_fn = argv[optind++];

  // FMD-index loading
  rt = realtime();
  rb3_fmi_t fmd;
  rb3_fmi_restore(&fmd, fmd_fn.c_str(), 0);
  if (fmd.e == 0 && fmd.r == 0) {
    fprintf(stderr, "[E::%s] FMD-index cannot be restored", __func__);
    exit(EXIT_FAILURE);
  }
  fprintf(stderr, "[M::%s] restored FMD index in %.3f sec\n", __func__,
          realtime() - rt);

  // Graph
  Graph graph(gbz_fn, "CHM13");
  graph.load();
  graph.load_fl();
  fprintf(stderr, "[M::%s] Load graph in in %.3f sec\n", __func__,
          realtime() - rt);

  // Sketch
  rt = realtime();
  sketch_t *sketch = sk_load(skt_fn);
  size_t klen = sketch->k;
  fprintf(stderr, "[M::%s] Restored sketch in %.3f sec\n", __func__,
          realtime() - rt);

  rbatch_t *rb = rbx_init(fx_fn.c_str(), bsize);
  std::vector<std::vector<sfs_t>> output(bsize);

  int total = 0, nreads = 0;
  int x;
  rt = realtime();
  while ((x = rbx_load(rb)) > 0) {
    nreads += x;
    // #pragma omp parallel for num_threads(nth) schedule(static, 1)
    for (int qq = 0; qq < x; ++qq) {
      uint8_t *seq = (uint8_t *)rb->reads[qq]->seq;
      int seql = rb->reads[qq]->seq_l;
      rb3_char2nt6(seql, seq);
      output[qq] = ping_pong_search(&fmd, seq, seql, rb->reads[qq]->name);
      assemble(output[qq], 0);

      if (anchoring) {
        std::map<uint64_t, int> kmers = count_kmers(seq, seql, sketch->k);
        anchor(graph, sketch, output[qq], seq, seql, kmers, klen, NA,
               reference_only);
        // remove_duplicates(sfs);
      }
      for (sfs_t &s : output[qq])
        fill_sequence(s, seq);
    }

    for (int qq = 0; qq < x; ++qq) {
      for (uint j = 0; j < output[qq].size(); ++j) {
        const sfs_t &s = output[qq][j];
        ++total;
        if (anchoring) {
          if (s.flag == 0) {
            std::cout << (int)s.flag << "\t" << s.rname << "\t" << s.s << "\t"
                      << s.l << "\t" << s.s + s.l << "\t" << s.sv << "\t"
                      << s.ev << "\t" << s.soff << "\t" << s.eoff << "\t"
                      << graph.get_gfa_name(s.sv >> 1) << "\t"
                      << graph.get_gfa_name(s.ev >> 1) << "\t" << s.skmer
                      << "\t" << s.ekmer << "\t";
            std::cout << s.paths[0];
            for (size_t p = 1; p < s.paths.size(); ++p)
              std::cout << "," << s.paths[p];
            std::cout << "\t" << s.plain_seq << std::endl;
          } else {
            std::cout << (int)s.flag << "\t" << s.rname << "\t" << s.s << "\t"
                      << s.l << "\t" << s.s + s.l << "\t"
                      << "."
                      << "\t"
                      << "."
                      << "\t"
                      << "."
                      << "\t"
                      << "."
                      << "\t"
                      << "."
                      << "\t"
                      << "."
                      << "\t"
                      << "."
                      << "\t"
                      << "."
                      << "\t";
            std::cout << ".";
            std::cout << "\t"
                      << "." << std::endl;
          }
        } else {
          std::cout << (int)s.flag << "\t" << s.rname << "\t" << s.s << "\t"
                    << s.l << "\t" << s.s + s.l << "\t" << s.plain_seq
                    << std::endl;
        }
      }
    }
    fprintf(stderr,
            "[M::%s] computed %d specific strings from %d reads in %.3f sec\n",
            __func__, total, nreads, realtime() - rt);
  }

  rbx_destroy(rb);
  sk_destroy(sketch);
  rb3_fmi_free(&fmd);

  return 0;
}
