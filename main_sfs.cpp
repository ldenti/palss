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

    out.push_back({name, begin, /*l - end - 1,*/ end - begin + 1});

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

// flip vertex strand bit
uint32_t flip(uint32_t v) { return v ^ 1; }

// // Get paths containing an anchor (since we store only start/end vertices, we
// // need to subset paths over bubbles)
// void set_paths(const Graph &graph, anchor_t &anchor, int klen) {
//   uint32_t v1 = anchor.v1;
//   uint32_t v2 = anchor.v2;

//   // Get paths on + strand only (1) - 3 is both, 2 is - only
//   std::vector<path_t> paths = graph.get_paths(v1, v2, 1, false);
//   std::vector<path_t> paths2 = graph.get_paths(flip(v2), flip(v1), 1, false);
//   for (path_t &p : paths) {
//     p.reversed = false;
//   }
//   for (path_t &p : paths2) {
//     p.reversed = true;
//     paths.push_back(p);
//   }

//   // get only paths containing the anchor kmer
//   char skmer[klen];
//   char skmer_rc[klen];
//   if (anchor.is_canonical) {
//     // the kmer stored in anchor is the same as in the read
//     d2s(anchor.kmer, klen, skmer);
//     d2s(rc(anchor.kmer, klen), klen, skmer_rc);
//   } else {
//     d2s(rc(anchor.kmer, klen), klen, skmer);
//     d2s(anchor.kmer, klen, skmer_rc);
//   }

//   // iterate over paths and find kmer and assign path to anchor if good
//   for (size_t i = 0; i < paths.size(); ++i) {
//     const path_t &path = paths[i];
//     size_t p1 = path.sequence.find(skmer);
//     size_t p0 = path.sequence.find(skmer_rc);
//     if (p1 != std::string::npos || p0 != std::string::npos) {
//       assert(p1 != p0);

//       // XXX: this assert may fail if we have cycles (but it depends if
//       // sketching filter out kmers occurring on same vertex/vertices two or
//       // more times along the same path)
//       assert(anchor.paths.find(path.id << 1) == anchor.paths.end() &&
//              anchor.paths.find((path.id << 1) | 1) == anchor.paths.end());

//       anchor.paths[(path.id << 1) | (p1 != std::string::npos)] = {
//           (uint32_t)path.offset1, (uint32_t)path.offset2,
//           p1 != std::string::npos, path.reversed, path.is_reference};
//       // key is path identifier w/ strand bit + consistency bit: "is the kmer
//       // along the path consistent with the kmer extracted from the read?"
//       (to
//       // model strand). They are consistent (bit to 1) if both canonical or
//       both
//       // non canonical
//     }
//   }
// }

// Get paths containing an anchor (since we store only start/end vertices, we
// need to subset paths over bubbles)
void set_paths(const Graph &graph, anchor_t &anchor, int klen) {
  uint32_t v1 = anchor.v1;
  uint32_t v2 = anchor.v2;

  char akmer_can[klen];
  char akmer_rc[klen];
  d2s(anchor.kmer, klen, akmer_can);
  d2s(rc(anchor.kmer, klen), klen, akmer_rc);

  // Get paths on + strand only (1) - 3 is both, 2 is - only

  // From v1 to v2 I must check for canonical kmer (since sketch stores
  // canonical direction)
  std::vector<path_t> paths = graph.get_paths(v1, v2, 1, false);
  for (path_t &p : paths) {
    if (p.sequence.find(akmer_can) != std::string::npos) {
      // if the anchor from read was already canonical, then strand is + (we
      // have consistency)
      anchor.paths[(p.id << 1) | anchor.is_canonical] = {
          (uint32_t)p.offset1, (uint32_t)p.offset2, anchor.is_canonical, false,
          p.is_reference};
    }
  }

  std::vector<path_t> paths2 = graph.get_paths(flip(v2), flip(v1), 1, false);
  for (path_t &p : paths2) {
    if (p.sequence.find(akmer_rc) != std::string::npos) {
      // if the anchor from read was canonical, then strand is - (we do not have
      // consistency since here we don't have the canonical kmer in the path for
      // sure, due to how sketching works)
      anchor.paths[(p.id << 1) | !anchor.is_canonical] = {
          (uint32_t)p.offset1, (uint32_t)p.offset2, !anchor.is_canonical, true,
          p.is_reference};
    }
  }

  // XXX: what about cycles?

  // in anchor.paths, the key is path identifier w/ strand bit + consistency
  // bit: "is the kmer along the path consistent with the kmer extracted from
  // the read?" (to model strand). They are consistent (bit to 1) if both
  // canonical or both non canonical
}

// Compute distance in bp between v1 and v2 on path with identifier pid
// XXX: can we get the "offset" directly when building the path so we do not
// need to decompressSA here again?
size_t compute_distance_bp(const Graph &graph, const gbwt::node_type &v1,
                           const gbwt::node_type &v2,
                           const gbwt::size_type &pid) {

  // std::cerr << "Computing distance: " << graph.get_gfa_name(v1 >> 1)
  //           << ((v1 & 1) ? "-" : "+") << " > " << graph.get_gfa_name(v2 >> 1)
  //           << ((v2 & 1) ? "-" : "+") << std::endl;
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
      if (position.first != gbwt::ENDMARKER)
        return -1UL;
      gbwtgraph::handle_t handle =
          gbwtgraph::GBWTGraph::node_to_handle(position.first);
      l += graph.get_gbz().graph.get_length(handle);
      position = graph.get_gbz().index.LF(position);
    }
    // skip v2
    return l;
  }
  assert(false);
  return -1UL;
}

// Chaining heuristic (both strands)
anchors_t chain(const anchors_t &anchors, const Graph &graph,
                const gbwt::size_type &pid, int klen, bool use_reverse) {
  anchors_t best_chain;
  anchors_t chain;

  for (size_t x = 0; x < anchors.size(); ++x) {
    chain.push_back(anchors.at(x));
    for (size_t y = x + 1; y < anchors.size(); ++y) {
      const anchor_t &yy = anchors.at(y);

      int qd = yy.qp - chain.back().qp;
      if (use_reverse)
        qd = yy.qp_rev - chain.back().qp_rev;
      assert(qd > 0);
      int rd = 0;
      if (chain.back().v1 == yy.v1) {
        // same vertex, so check positions
        rd = yy.pos1 - chain.back().pos1;
      } else {
        // different vertex, so check offset along path
        int vl = graph.get_vertex_len(chain.back().v1 >> 1);
        size_t d =
            compute_distance_bp(graph, chain.back().v1, yy.v1,
                                pid >> 1); // we need to remove consistency bit
        if (d == -1UL)
          rd = 0;
        else
          rd = d + vl - chain.back().pos1 + yy.pos1;
        assert(rd >= 0);
      }

      if (rd <= 0)
        // XXX1: if rd < 0, for sure we are on the same vertex and we might have
        // some strange event
        // XXX2: since anchors are unique on the graph AND the read, we cannot
        // have distance 0 between them. Even if we have a cycle, we filtered
        // out anchors that are repeated on the read
        continue;

      if (std::max(rd, qd) == 0 || std::min(rd, qd) / (float)std::max(rd, qd) >
                                       0.9) { // FIXME: hardcoded
        chain.push_back(yy);
      }
    }
    if (chain.size() / (float)anchors.size() > 0.5) {
      best_chain = chain;
      break;
    } else {
      if (chain.size() > best_chain.size()) {
        best_chain = chain;
      }
    }
    chain.clear();
  }
  return best_chain;
}

typedef struct {
  anchor_t anchor1;
  anchor_t anchor;
  anchor_t anchor2;
  int chain_size;
  int total_anchors;
  // bool strand;
} anchoring_t;

std::map<gbwt::size_type, anchoring_t>
chaining(const Graph &graph, const anchors_t &anchors, int klen) {
  std::map<gbwt::size_type, anchoring_t> result;

  // Get how many anchors we have on "each" path
  std::map<uint32_t, std::vector<uint32_t>> pcounts;
  for (size_t a = 0; a < anchors.size(); ++a) {
    const anchor_t &aa = anchors[a];
    for (const auto &[pid, p] : aa.paths)
      pcounts[pid].push_back(a);
  }

  // XXX: we can't merge paths with same anchors since they might have different
  // vertices in between

  // Sort paths by decreasing number of anchors
  std::vector<std::pair<uint32_t, uint32_t>> sorted_pcounts;
  for (const auto &[k, v] : pcounts)
    sorted_pcounts.push_back(std::make_pair(v.size(), k));
  std::sort(sorted_pcounts.begin(), sorted_pcounts.end(),
            [](const auto &a, const auto &b) { return a > b; });

  // Iterate over paths and chain
  for (const auto &[na, p] : sorted_pcounts) {
    // p is path id (w/ strand bit) plus "consistency" bit
    anchors_t path_anchors;
    for (const size_t aa : pcounts[p]) {
      path_anchors.push_back(anchors[aa]);
      if (path_anchors.back().paths[p].reversed) {
        std::swap(path_anchors.back().v1, path_anchors.back().v2);
        path_anchors.back().v1 = flip(path_anchors.back().v1);
        path_anchors.back().v2 = flip(path_anchors.back().v2);

        std::swap(path_anchors.back().pos1, path_anchors.back().pos2);
        size_t vl1 = graph.get_vertex_len(path_anchors.back().v1 >> 1);
        size_t vl2 = graph.get_vertex_len(path_anchors.back().v2 >> 1);
        path_anchors.back().pos1 = vl1 - path_anchors.back().pos1 - 1;
        path_anchors.back().pos2 = vl2 - path_anchors.back().pos2 - 1;
      }
    }

    if (!(p & 1))
      // read on reverse strand of path
      std::reverse(path_anchors.begin(), path_anchors.end());

    // std::pair<anchors_t, bool> bc =
    anchors_t bc = chain(path_anchors, graph, p, klen, !(p & 1));
    if (!(p & 1))
      // read on reverse strand, so chain is "inverted"
      result[p] = {bc.back(), bc[bc.size() / 2], bc.front(), (int)bc.size(),
                   (int)path_anchors.size()};
    else
      result[p] = {bc.front(), bc[bc.size() / 2], bc.back(), (int)bc.size(),
                   (int)path_anchors.size()};
  }
  return result;
}

// Anchor specific strings on graph using graph sketch
void anchor(const Graph &graph, sketch_t *sketch, std::vector<sfs_t> &sfs,
            uint8_t *read, int readl, std::map<uint64_t, int> kcounts, int klen,
            size_t NA, bool reference_only) {
  int beg, end;
  uint8_t kmer[klen];
  uint64_t kmer_d;       // kmer
  uint64_t rckmer_d = 0; // reverse and complemented kmer
  uint64_t ckmer_d = 0;  // canonical kmer
  int c;                 // current char
  hit_t hit;             // hit from sketch
  int count;

  for (uint sidx = 0; sidx < sfs.size(); ++sidx) {
    sfs_t &s = sfs[sidx];

    // Finding anchors on the left
    anchors_t sanchors;
    beg = s.s -
          klen; // XXX: we do not want anchors overlapping the specific string
    if (beg >= 0) {
      memcpy(kmer, read + beg, klen);
      kmer[klen] = '\0';
      kmer_d = k2d((char *)kmer, klen);
      rckmer_d = rc(kmer_d, klen);
      ckmer_d = std::min(kmer_d, rckmer_d);

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
          a.qp_rev = readl - (beg + klen);
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
    }

    // Finding anchors on the right
    anchors_t eanchors;
    end = s.s +
          s.l; // XXX: we do not want anchors overlapping the specific string
    if (end <= readl - klen) {
      memcpy(kmer, read + end, klen);
      kmer[klen] = '\0';
      kmer_d = k2d((char *)kmer, klen);
      rckmer_d = rc(kmer_d, klen);
      ckmer_d = std::min(kmer_d, rckmer_d);

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
          a.qp_rev = readl - (end + klen);
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
    }

    if (sanchors.size() == 0 || eanchors.size() == 0) {
      s.flag |= 1; // flag as invalid
      continue;
    }

    std::reverse(sanchors.begin(), sanchors.end());
    std::map<gbwt::size_type, anchoring_t> schains =
        chaining(graph, sanchors, klen);
    std::map<gbwt::size_type, anchoring_t> echains =
        chaining(graph, eanchors, klen);

    std::map<std::pair<uint32_t, uint32_t>, sfs_t> anchoring;
    for (const auto &[p, sanchoring] : schains) {
      // p is path id (w/ strand bit) plus "consistency" bit
      auto x = echains.find(p);
      if (x == echains.end())
        continue;

      // XXX: if we want anchor in between the chain, use .anchor
      anchor_t a1 = sanchoring.anchor2;
      anchor_t a2 = x->second.anchor1;

      uint32_t qp = a1.qp;
      uint32_t l = a2.qp + klen - a1.qp;
      std::pair<uint32_t, uint32_t> key = std::make_pair(qp, l);
      if (anchoring.find(key) == anchoring.end()) {
        anchoring[key] = {};
        sfs_t &s = anchoring[key];
        s.s = qp;
        s.l = l;
        // s.strand = sanchoring.strand;
        if (p & 1) {
          s.sv = a1.v1;
          s.ev = a2.v2;
          s.soff = a1.pos1;
          s.eoff = a2.pos2;
        } else {
          s.sv = a2.v1;
          s.soff = a2.pos1;
          s.ev = a1.v2;
          s.eoff = a1.pos2;
          // s.sv = ((s.sv >> 1) << 1) | (~(s.sv & 1) & 1);
          // s.ev = ((s.ev >> 1) << 1) | (~(s.ev & 1) & 1);
          // s.skmer = a2.kmer;
          // s.ekmer = a1.kmer;
        }
        s.skmer = a1.kmer;
        s.ekmer = a2.kmer;
      }
      anchoring[key].paths.push_back(p);
    }

    if (anchoring.empty()) {
      s.flag |= 2; // flag as invalid
      continue;
    }

    uint32_t best = 0;
    std::pair<uint32_t, uint32_t> best_key;
    for (const auto &[k, s] : anchoring) {
      if (s.paths.size() > best) {
        best = s.paths.size();
        best_key = k;
      }
    }

    // assert(best_anchors.first.qp < best_anchors.second.qp);

    // Assigning the anchors
    s.s = anchoring[best_key].s;
    s.l = anchoring[best_key].l;
    // s.strand = anchoring[best_key].strand;
    s.sv = anchoring[best_key].sv;
    s.ev = anchoring[best_key].ev;
    s.soff = anchoring[best_key].soff;
    s.eoff = anchoring[best_key].eoff;
    s.paths = anchoring[best_key].paths;
    s.skmer = anchoring[best_key].skmer;
    s.ekmer = anchoring[best_key].ekmer;
  }
}

void remove_duplicates(std::vector<sfs_t> &S) {
  // sfs are sorted by position on read
  size_t last_i = 0;
  while (last_i < S.size() && S[last_i].flag != 0)
    ++last_i;
  for (size_t i = last_i + 1; i < S.size(); ++i) {
    if (S[i].flag != 0)
      continue;
    if (S[i].s == S[last_i].s) {
      if (S[i].l < S[last_i].l) {
        S[i].flag = 3;
      } else {
        S[last_i].flag = 3;
        last_i = i;
      }
    } else if (S[i].s + S[i].l == S[last_i].s + S[last_i].l) {
      S[i].flag = 3;
    } else {
      last_i = i;
    }
  }
}

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

  // Reads
  rbatch_t *rb = rbx_init(fx_fn.c_str(), bsize);
  std::vector<std::vector<sfs_t>> output(bsize);

  int total = 0, nreads = 0;
  int x;
  rt = realtime();
  while ((x = rbx_load(rb)) > 0) {
    nreads += x;
    // #pragma omp parallel for num_threads(nth) schedule(static, 1)
    for (int qq = 0; qq < x; ++qq) {
      std::cerr << "=== " << rb->reads[qq]->name << " ===" << std::endl;
      // convert to 1234
      uint8_t *seq = (uint8_t *)rb->reads[qq]->seq;
      int seql = rb->reads[qq]->seq_l;
      rb3_char2nt6(seql, seq);

      // search specific strings
      output[qq] = ping_pong_search(&fmd, seq, seql, rb->reads[qq]->name);

      // assemble overlapping specific strings
      assemble(output[qq], 0);

      // anchor specific strings to graph
      if (anchoring) {
        std::map<uint64_t, int> kmers = count_kmers(seq, seql, sketch->k);
        anchor(graph, sketch, output[qq], seq, seql, kmers, klen, NA,
               reference_only);
        remove_duplicates(output[qq]);
      }
      // fill sequence
      for (sfs_t &s : output[qq])
        fill_sequence(s, seq);
    }

    // output
    // TODO: use a thread to do this (as in SVDSS)
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
            std::cout << s.paths[0] << "|"
                      << graph.get_path_sample(s.paths[0] >> 2) << "|"
                      << graph.get_path_contig(s.paths[0] >> 2) << "|"
                      << (s.paths[0] & 1);
            for (size_t p = 1; p < s.paths.size(); ++p)
              std::cout << "," << s.paths[p] << "|"
                        << graph.get_path_sample(s.paths[p] >> 2) << "|"
                        << graph.get_path_contig(s.paths[p] >> 2) << "|"
                        << (s.paths[p] & 1);
            std::cout << "\t" << s.plain_seq << std::endl;
          } else {
            std::cout << (int)s.flag << "\t" << s.rname << "\t" << s.s << "\t"
                      << s.l << "\t" << s.s + s.l << "\t" << "."
                      << "\t" << "."
                      << "\t" << "." << "\t" << "." << "\t" << "." << "\t"
                      << "." << "\t" << "." << "\t" << "." << "\t";
            std::cout << ".";
            std::cout << "\t" << "." << std::endl;
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
