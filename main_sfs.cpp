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

// Get paths containing an anchor (since we store only start/end vertices, we
// need to subset paths over bubbles)
void set_paths(const Graph &graph, anchor_t &anchor, int klen, size_t NP) {
  // XXX: what about cycles?

  anchor.paths.clear();

  uint32_t v1 = anchor.v1;
  uint32_t v2 = anchor.v2;

  char akmer_can[klen];
  char akmer[klen];
  d2s(anchor.kmer, klen, akmer_can);
  d2s(rc(anchor.kmer, klen), klen, akmer);
  if (anchor.inverted) {
    // we can enter here here only after chaining, since we set the inverted bit
    // only after chaining. This means that the anchor was inverted along the
    // path, so we need to flip the vertices to get the correct orientation of
    // the anchor along the path
    v1 = flip(anchor.v2);
    v2 = flip(anchor.v1);
  }

  // Get paths on + strand only (1) - 3 is both, 2 is - only

  // From v1 to v2 I must check for canonical kmer (since sketch stores
  // canonical direction)
  std::vector<path_t> paths = graph.get_paths(v1, v2, 1, NP, false);
  for (path_t &p : paths) {
    if (p.sequence.find(akmer_can) != std::string::npos) {
      // if the anchor from read was already canonical, then strand is + (we
      // have consistency)
      anchor.paths[(p.id << 3) | ((p.is_reference & 1) << 2) |
                   anchor.is_canonical << 1 | 0] =
          ((uint64_t)p.offset1 << 32) | (uint32_t)p.offset2;
    }
  }

  std::vector<path_t> paths2 =
      graph.get_paths(flip(v2), flip(v1), 1, NP, false);
  for (path_t &p : paths2) {
    if (p.sequence.find(akmer) != std::string::npos) {
      // if the anchor from read was canonical, then strand is - (we do not have
      // consistency since here we don't have the canonical kmer in the path for
      // sure, due to how sketching works)
      anchor.paths[(p.id << 3) | ((p.is_reference & 1) << 2) |
                   !anchor.is_canonical << 1 | 1] =
          ((uint64_t)p.offset1 << 32) | (uint32_t)p.offset2;
    }
  }
}

// Compute distance in bp between v1 and v2 on path with identifier pid
// XXX: can we get the "offset" directly when building the path so we do not
// need to decompressSA here again?
size_t compute_distance_bp(const Graph &graph, const gbwt::node_type &v1,
                           const gbwt::node_type &v2,
                           const gbwt::size_type &pid) {

  // std::cerr << "Computing distance: " << graph.get_gfa_name(v1 >> 1)
  //           << ((v1 & 1) ? "-" : "+") << " > " << graph.get_gfa_name(v2 >> 1)
  //           << ((v2 & 1) ? "-" : "+") << " on "
  //           << graph.get_path_contig(pid >> 1) << std::endl;
  if (v1 == v2)
    return 0;

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
      if (position.first == gbwt::ENDMARKER)
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
std::vector<size_t>
chain(const anchors_t &anchors, const std::vector<size_t> &anchor_paths,
      const Graph &graph, const gbwt::size_type &pid, int klen,
      bool use_reverse,
      std::map<gbwt::size_type, std::map<uint64_t, size_t>> &dmemo) {
  std::vector<size_t> best_chain;
  std::vector<size_t> chain;

  uint64_t xy;
  for (size_t x = 0; x < anchor_paths.size(); ++x) {
    chain.push_back(anchor_paths[x]);
    for (size_t y = x + 1; y < anchor_paths.size(); ++y) {
      const anchor_t &yy = anchors.at(anchor_paths[y]);

      const anchor_t &back = anchors[chain.back()];
      int qd = yy.qp - back.qp;
      if (use_reverse)
        qd = yy.qp_rev - back.qp_rev;
      assert(qd > 0);
      int rd = 0;
      if (back.v1 == yy.v1) {
        // same vertex, so check positions
        rd = yy.pos1 - back.pos1;
      } else {
        // different vertex, so check offset along path
        int vl = graph.get_vertex_len(back.v1 >> 1);
        if (dmemo.find(pid) == dmemo.end()) {
          dmemo[pid] = std::map<uint64_t, size_t>();
        }
        xy = ((uint64_t)back.v1 << 32) | (uint32_t)yy.v1;
        size_t dist;
        if (dmemo[pid].find(xy) != dmemo[pid].end()) {
          dist = dmemo[pid][xy];
        } else {
          dist = compute_distance_bp(graph, back.v1, yy.v1, pid);
          dmemo[pid][xy] = dist;
        }
        if (dist == -1UL)
          rd = 0;
        else
          rd = dist + vl - back.pos1 + yy.pos1;
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
        chain.push_back(anchor_paths[y]);
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

anchoring_t chaining(const Graph &graph, anchors_t &anchors,
                     size_t selected_path, int klen) {

  anchoring_t result;

  std::map<gbwt::size_type, std::map<uint64_t, size_t>> dmemo;

  // Get how many anchors we have on "each" path
  std::map<uint32_t, std::vector<uint32_t>> pcounts;
  for (size_t a = 0; a < anchors.size(); ++a) {
    const anchor_t &aa = anchors[a];
    for (const auto &[pid, p] : aa.paths)
      pcounts[pid >> 1].push_back(
          (a << 1) | (pid & 1)); // We need to remove the inverted bit before
                                 // hashing. We store it in the value
  }

  // NOTE: we can't merge paths with same anchors since they might have
  // different vertices in between

  // Sort paths by decreasing number of anchors
  std::vector<std::pair<uint32_t, uint32_t>> sorted_pcounts;
  for (const auto &[k, v] : pcounts)
    sorted_pcounts.push_back(std::make_pair(v.size(), k));
  std::sort(sorted_pcounts.begin(), sorted_pcounts.end(),
            [](const auto &a, const auto &b) { return a > b; });

  // Iterate over paths and chain
  for (const auto &[na, p] : sorted_pcounts) {
    if (p != selected_path)
      continue;

    // p is path id (w/ strand bit) plus two additional bits (we moved the
    // inverted bit to values)
    std::vector<size_t> path_anchors; // indices along anchors

    for (const size_t aa : pcounts[p]) {
      anchor_t &a = anchors[aa >> 1];
      path_anchors.push_back(aa >> 1); // remove inverted bit

      if (aa & 1) {
        // invert if inverted bit is set
        a.inverted = true;

        std::swap(a.v1, a.v2);
        a.v1 = flip(a.v1);
        a.v2 = flip(a.v2);

        std::swap(a.pos1, a.pos2);
        size_t vl1 = graph.get_vertex_len(a.v1 >> 1);
        size_t vl2 = graph.get_vertex_len(a.v2 >> 1);
        a.pos1 = vl1 - a.pos1 - 1;
        a.pos2 = vl2 - a.pos2 - 1;
      }
    }

    bool reverse_strand = (p & 1) == 0;
    if (reverse_strand)
      // read on reverse strand of path
      std::reverse(path_anchors.begin(), path_anchors.end());

    // regardless of strand, anchors are now sorted wrt path. If strand is
    // reverse, we will use the reverse read positions
    std::vector<size_t> bc = chain(anchors, path_anchors, graph, p >> 2, klen,
                                   reverse_strand, dmemo);

    // we want to report first and last anchors along the chain, sorted by
    // original read position
    if (reverse_strand)
      // read on reverse strand, so we need to reverse the chain
      result = {anchors[bc.back()], anchors[bc[bc.size() / 2]],
                anchors[bc.front()], (int)bc.size(), (int)path_anchors.size()};
    else
      result = {anchors[bc.front()], anchors[bc[bc.size() / 2]],
                anchors[bc.back()], (int)bc.size(), (int)path_anchors.size()};
  }
  return result;
}

// Anchor specific strings on graph using graph sketch
void anchor(const Graph &graph, sketch_t *sketch, std::vector<sfs_t> &sfs,
            uint8_t *read, int readl, std::map<uint64_t, int> kcounts, int klen,
            size_t NA, size_t NP, bool reference_only, bool overlapping) {
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

    // char kmer_s[klen + 1];
    // kmer_s[klen] = '\0';
    // d2s(kmer_d, klen, kmer_s);
    // std::cerr << kmer_s << std::endl;

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
          a.inverted = false;

          sanchors.push_back(a);
          set_paths(graph, sanchors.back(), klen, NP);
        }

        --beg;

        if (overlapping || hit.value == -1UL || count != 1) {
          if (beg < 0)
            break;
          c = read[beg] - 1;
          kmer_d = rsprepend(kmer_d, c, klen);
          rckmer_d = lsappend(rckmer_d, reverse_char(c), klen);
          ckmer_d = std::min(kmer_d, rckmer_d);
        } else {
          beg = beg - klen + 1;
          if (beg < 0)
            break;
          memcpy(kmer, read + beg, klen);
          kmer[klen] = '\0';
          kmer_d = k2d((char *)kmer, klen);
          rckmer_d = rc(kmer_d, klen);
          ckmer_d = std::min(kmer_d, rckmer_d);
        }
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
          a.inverted = false;

          eanchors.push_back(a);
          set_paths(graph, eanchors.back(), klen, NP);
        }

        ++end;

        if (overlapping || hit.value == -1UL || count != 1) {
          if (end >= readl - klen + 1)
            break;
          c = read[end + klen - 1] - 1;
          kmer_d = lsappend(kmer_d, c, klen);
          rckmer_d = rsprepend(rckmer_d, reverse_char(c), klen);
          ckmer_d = std::min(kmer_d, rckmer_d);
        } else {
          end = end + klen - 1;
          if (end >= readl - klen + 1)
            break;
          memcpy(kmer, read + end, klen);
          kmer[klen] = '\0';
          kmer_d = k2d((char *)kmer, klen);
          rckmer_d = rc(kmer_d, klen);
          ckmer_d = std::min(kmer_d, rckmer_d);
        }
      }
    }

    /* Here, a set of paths is assigned to each anchor. Each path encodes:
        - path identifier (w/ strand bit)
        - reference bit (1 if path is reference, 0 otherwise)
        - consistency bit (1 if anchor is consistent with path, 0 otherwise)
        - inverted bit (1 if anchor is inverted along path, 0 otherwise)
    */

    if (sanchors.size() == 0 || eanchors.size() == 0) {
      s.flag = 1; // flag as invalid
      continue;
    }

    // Get how many anchors we have on "each" path
    std::map<uint32_t, size_t> spcounts;
    for (size_t a = 0; a < sanchors.size(); ++a) {
      const anchor_t &aa = sanchors[a];
      for (const auto &[pid, p] : aa.paths)
        ++spcounts[pid >> 1]; // We need to remove the inverted bit before
                              // hashing, since it is "path-specific".
                              // Consistency bit gives us information on the
                              // strand of the anchor along the path
    }

    std::map<uint32_t, size_t> epcounts;
    for (size_t a = 0; a < eanchors.size(); ++a) {
      const anchor_t &aa = eanchors[a];
      for (const auto &[pid, p] : aa.paths)
        ++epcounts[pid >> 1]; // We need to remove the inverted bit before
                              // hashing, since it is "path-specific".
                              // Consistency bit gives us information on the
                              // strand of the anchor along the path
    }

    // Find paths that have anchors on both sides of the specific string
    std::vector<uint32_t> paths;
    for (const auto &[k, v] : spcounts) {
      if (epcounts.find(k) != epcounts.end()) {
        paths.push_back(k);
      }
    }

    if (paths.empty()) {
      s.flag = 2; // flag as invalid
      continue;
    }

    // select "best" path
    std::sort(
        paths.begin(), paths.end(),
        [&spcounts, &epcounts](auto const &a, auto const &b) {
          double ha =
              2.0 * spcounts[a] * epcounts[a] / (spcounts[a] + epcounts[a]);
          double hb =
              2.0 * spcounts[b] * epcounts[b] / (spcounts[b] + epcounts[b]);
          if (ha != hb)
            return ha > hb; // prefer higher harmonic mean (balanced+large)
          else if ((spcounts[a] + epcounts[a]) != (spcounts[b] + epcounts[b]))
            return (spcounts[a] + epcounts[a]) >
                   (spcounts[b] +
                    epcounts[b]); // tie-breaker #1: prefer larger sum
          else
            return a < b; // tie-breaker #2: prefer lower ID (reference
                          // path is the lowest)
        });

    uint32_t selected_path = paths.front(); // "best" path

    // reverse start anchors so that they follow read order
    // XXX: this might create problems with our greedy strategy since we start
    // from farthest anchor
    std::reverse(sanchors.begin(), sanchors.end());
    anchoring_t schains = chaining(graph, sanchors, selected_path, klen);
    anchoring_t echains = chaining(graph, eanchors, selected_path, klen);

    // chains are reported following original read positions

    // We could use anchor in between the chain (.anchor)
    anchor_t a1 = schains.anchor2;
    anchor_t a2 = echains.anchor1;

    if (NP != (size_t)-1) {
      set_paths(graph, a1, klen, -1);
      set_paths(graph, a2, klen, -1);
    }

    s.flag = 0;
    s.s = a1.qp;
    s.l = a2.qp + klen - a1.qp;

    // XXX: vvvvvvv check this vvvvvvv
    if (selected_path & 1) {
      // we are on forward strand - there is no inverted bit
      s.sv = a1.v1;
      s.soff = a1.pos1;
      s.ev = a2.v2;
      s.eoff = a2.pos2;
      s.skmer = a1.kmer;
      s.ekmer = a2.kmer;
      s.reverse = false;
    } else {
      // we are on reverse strand - there is no inverted bit
      s.sv = a1.v2;
      s.soff = a1.pos2;
      s.ev = a2.v1;
      s.eoff = a2.pos1;
      s.skmer = a2.kmer;
      s.ekmer = a1.kmer;
      s.reverse = true;
    }
    if (a1.inverted) {
      s.sv = flip(s.sv);
      s.soff = graph.get_vertex_len(s.sv >> 1) - s.soff - 1;
    }
    if (a2.inverted) {
      s.ev = flip(s.ev);
      s.eoff = graph.get_vertex_len(s.ev >> 1) - s.eoff - 1;
    }
    // XXX: ^^^^^^^ check this ^^^^^^^

    // set paths
    for (const auto &[pid, p] : a1.paths) {
      uint8_t a1_is_inverted = pid & 1;
      uint8_t a2_is_inverted;
      if (a2.paths.find(pid) != a2.paths.end())
        a2_is_inverted = a1_is_inverted;
      else if (a2.paths.find(pid ^ 1) != a2.paths.end())
        a2_is_inverted = !a1_is_inverted;
      else
        continue;

      // the new path id will be: path identifier with strand bit + reference
      // bit + consistency/strand bit (1 is +, 0 is -) + inverted anchor1
      // + inverted anchor2

      int pp = (pid << 1) | (a1_is_inverted << 1) | a2_is_inverted;
      s.paths.push_back(pp);
    }

    // we have to check the paths since we might have linked bad chains
    std::vector<uint64_t> path_ids_to_keep;
    for (const uint64_t &pt : s.paths) {
      uint32_t sv = s.sv;
      uint32_t ev = s.ev;
      uint32_t soff = s.soff;
      uint32_t eoff = s.eoff;
      if ((pt >> 1) & 1) {
        sv = s.sv ^ 1;
        soff = graph.get_vertex_len(sv >> 1) - soff - 1;
      }
      if (pt & 1) {
        ev = s.ev ^ 1;
        eoff = graph.get_vertex_len(ev >> 1) - eoff - 1;
      }
      if (((pt >> 2) & 1) == 0) {
        std::swap(sv, ev);
        std::swap(soff, eoff);
      }

      // if two vertices are too distant or are "inverted" due to bad chain
      // linking over the specific string (we have -1UL in this case)
      size_t dist;
      if (sv == ev) {
        assert(soff != eoff);
        if (soff <= eoff)
          dist = eoff - soff;
        else
          dist = -1UL;
      } else {
        dist = compute_distance_bp(graph, sv, ev, pt >> 4);
      }
      if (dist <= 20000) // XXX: hardcoded
        path_ids_to_keep.push_back(pt);
    }

    if (path_ids_to_keep.empty()) {
      s.flag = 3; // flag as invalid
      continue;
    }
    s.paths = path_ids_to_keep;
  }
}

// Flag specific strings that are contained in other specific strings on the
// same read
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
        S[i].flag = 4;
      } else {
        S[last_i].flag = 4;
        last_i = i;
      }
    } else if (S[i].s + S[i].l == S[last_i].s + S[last_i].l) {
      S[i].flag = 4;
    } else {
      last_i = i;
    }
  }
}

void fill_sequence(sfs_t &s, uint8_t *read) {
  s.plain_seq = std::string((char *)read + s.s, s.l);
  std::string map = "0ACGT0";
  if (s.reverse) {
    std::reverse(s.plain_seq.begin(), s.plain_seq.end());
    std::reverse(map.begin(), map.end());
  }
  for (char &c : s.plain_seq)
    c = map[(uint8_t)c];
}

int main_sfs(int argc, char *argv[]) {
  double rt;
  int bsize = 10000; // batch size
  size_t NA = 20;    // number of kmers to check for anchoring
  size_t NP = -1;    // number of paths to check per anchor (default: all paths)
  int nth = 4;       // number of threads
  bool reference_only = false;
  bool search_only = false;
  bool overlapping = true;

  int _c;
  while ((_c = getopt(argc, argv, "a:p:b:sor@:h")) != -1) {
    switch (_c) {
    case 'a':
      NA = std::stoi(optarg);
      break;
    case 'p':
      NP = std::stoi(optarg);
      break;
    case 'b':
      bsize = std::stoi(optarg);
      break;
    case 's':
      search_only = true;
      break;
    case 'o':
      overlapping = false;
      break;
    case 'r':
      reference_only = true;
      break;
    case '@':
      nth = std::stoi(optarg);
      break;
    case 'h':
      fprintf(stderr, "%s", SFS_USAGE_MESSAGE);
      return 0;
    }
  }

  if (argc - optind != 4) {
    fprintf(stderr, "%s", SFS_USAGE_MESSAGE);
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

  // int total = 0;
  int nreads = 0;
  int x;
  rt = realtime();
  while ((x = rbx_load(rb)) > 0) {
    nreads += x;
#pragma omp parallel for num_threads(nth) schedule(static, 1)
    for (int qq = 0; qq < x; ++qq) {
      // std::cerr << "=== " << rb->reads[qq]->name << " ===" << std::endl;
      // convert to 1234
      uint8_t *seq = (uint8_t *)rb->reads[qq]->seq;
      int seql = rb->reads[qq]->seq_l;
      rb3_char2nt6(seql, seq);

      // search specific strings
      output[qq] = ping_pong_search(&fmd, seq, seql, rb->reads[qq]->name);
      // std::cerr << rb->reads[qq]->name << std::endl;

      // assemble overlapping specific strings
      assemble(output[qq], 0);

      // anchor specific strings to graph
      if (!search_only) {
        std::map<uint64_t, int> kmers = count_kmers(seq, seql, sketch->k);
        anchor(graph, sketch, output[qq], seq, seql, kmers, klen, NA, NP,
               reference_only, overlapping);
        remove_duplicates(output[qq]);
      }

      // fill sequence
      for (sfs_t &s : output[qq]) {
        fill_sequence(s, seq);

        if (!search_only & (s.flag == 0)) {
          std::cout << sfs_to_string(s, graph.get_gfa_name(s.sv >> 1),
                                     graph.get_gfa_name(s.ev >> 1));
        } else {
          std::cout << sfs_to_string(s, "", "");
        }
      }
    }
    fprintf(stderr,
            "[M::%s] Processed %d reads in %.3f sec (current rss: %lldGB)\n",
            __func__, nreads, realtime() - rt, current_rss_kb() / 1024 / 1024);
  }

  rbx_destroy(rb);
  sk_destroy(sketch);
  rb3_fmi_free(&fmd);

  return 0;
}
