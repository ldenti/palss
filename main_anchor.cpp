#include <fstream>
#include <getopt.h>
#include <omp.h>
#include <string>
#include <vector>
#include <zlib.h>

extern "C" {
#include "kseq.h"
}

#include "graph.hpp"
#include "kmer.hpp"
#include "misc.hpp"
#include "reads.hpp"
#include "sfs.hpp"
#include "sketch.hpp"
#include "usage.hpp"

typedef struct {
  uint64_t v;   // vertex on graph
  int p;        // position on query
  uint64_t seq; // kmer
} anchor_t;

// ACGT-read to kmer counts
std::map<uint64_t, int> count_kmers(char *read, int readl, int klen) {
  int c; // current char
  char *kmer = (char *)malloc(klen);
  std::map<uint64_t, int> kcounts;
  memcpy(kmer, read, klen);
  uint64_t kmer_d = k2d(kmer, klen);
  uint64_t rckmer_d = rc(kmer_d, klen);
  uint64_t ckmer_d = std::min(kmer_d, rckmer_d);
  ++kcounts[ckmer_d];
  for (int p = klen; p < readl; ++p) {
    c = to_int[(int)read[p]] - 1;
    kmer_d = lsappend(kmer_d, c, klen);
    rckmer_d = rsprepend(rckmer_d, reverse_char(c), klen);
    ckmer_d = std::min(kmer_d, rckmer_d);
    ++kcounts[ckmer_d];
  }
  free(kmer);
  return kcounts;
}

// Anchor specific strings on graph using graph sketch
void anchor(const Graph &graph, sketch_t *sketch, std::vector<sfs_t> &sfs,
            char *read, int readl, std::map<uint64_t, int> kcounts, int klen,
            uint NA, bool reference_only, bool overlapping) {

  int beg, end;
  uint8_t *kmer = (uint8_t *)malloc(klen);
  uint64_t kmer_d;       // kmer
  uint64_t rckmer_d = 0; // reverse and complemented kmer
  uint64_t ckmer_d = 0;  // canonical kmer
  int c;                 // current char
  uint64_t hit;          // hit from sketch

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
    std::vector<anchor_t> sanchors;
    while (beg > 0 && sanchors.size() < NA) {
      if (overlapping) {
        if ((hit = sk_get(sketch, ckmer_d)) != -1UL && kcounts[ckmer_d] == 1)
          sanchors.push_back({hit, beg, ckmer_d});
        --beg;
        c = to_int[(int)read[beg]] - 1;
        kmer_d = rsprepend(kmer_d, c, klen);
        rckmer_d = lsappend(rckmer_d, reverse_char(c), klen);
        ckmer_d = std::min(kmer_d, rckmer_d);
      } else {
        if ((hit = sk_get(sketch, ckmer_d)) != -1UL && kcounts[ckmer_d] == 1) {
          sanchors.push_back({hit, beg, ckmer_d});
          beg -= klen;
          if (beg > 0) {
            memcpy(kmer, read + beg, klen);
            kmer_d = k2d((char *)kmer, klen);
            rckmer_d = rc(kmer_d, klen);
            ckmer_d = std::min(kmer_d, rckmer_d);
          }
        } else {
          --beg;
          c = to_int[(int)read[beg]] - 1;
          kmer_d = rsprepend(kmer_d, c, klen);
          rckmer_d = lsappend(rckmer_d, reverse_char(c), klen);
          ckmer_d = std::min(kmer_d, rckmer_d);
        }
      }
    }

    // Finding anchors on the right
    end = s.s + s.l - 1;
    end = end > readl - klen ? readl - klen : end;
    memcpy(kmer, read + end, klen);
    kmer_d = k2d((char *)kmer, klen);
    rckmer_d = rc(kmer_d, klen);
    ckmer_d = std::min(kmer_d, rckmer_d);
    std::vector<anchor_t> eanchors;
    while (end < readl - klen + 1 && eanchors.size() < NA) {
      if (overlapping) {
        if ((hit = sk_get(sketch, ckmer_d)) != -1UL && kcounts[ckmer_d] == 1)
          eanchors.push_back({hit, end, ckmer_d});
        ++end;
        c = to_int[(int)read[end + klen - 1]] - 1;
        kmer_d = lsappend(kmer_d, c, klen);
        rckmer_d = rsprepend(rckmer_d, reverse_char(c), klen);
        ckmer_d = std::min(kmer_d, rckmer_d);
      } else {
        if ((hit = sk_get(sketch, ckmer_d)) != -1UL && kcounts[ckmer_d] == 1) {
          eanchors.push_back({hit, end, ckmer_d});
          end += klen;
          if (end < readl - klen + 1) {
            memcpy(kmer, read + end, klen);
            kmer_d = k2d((char *)kmer, klen);
            rckmer_d = rc(kmer_d, klen);
            ckmer_d = std::min(kmer_d, rckmer_d);
          }
        } else {
          ++end;
          c = to_int[(int)read[end + klen - 1]] - 1;
          kmer_d = lsappend(kmer_d, c, klen);
          rckmer_d = rsprepend(rckmer_d, reverse_char(c), klen);
          ckmer_d = std::min(kmer_d, rckmer_d);
        }
      }
    }

    if (sanchors.size() == 0 || eanchors.size() == 0) {
      s.flag |= 1; // flag as invalid
      continue;
    }

    // Finding best pair of anchors
    size_t sa_i, ea_i;
    for (sa_i = 0; sa_i < sanchors.size(); ++sa_i) {
      anchor_t &sa = sanchors[sa_i];
      gbwt::size_type sv1 = sa.v >> 32;
      gbwt::size_type sv2 = (uint32_t)sa.v;
      std::vector<path_t> spaths = graph.get_paths(sv1, sv2, reference_only);

      // get only paths containing the anchor kmer
      char skmer1[klen];
      char skmer1_rc[klen];
      memcpy(skmer1, read + sa.p, klen);
      kmer_d = k2d((char *)skmer1, klen);
      rckmer_d = rc(kmer_d, klen);
      d2s(rckmer_d, klen, (char *)skmer1_rc);

      std::vector<std::pair<size_t, bool>> good_spaths;
      for (size_t i = 0; i < spaths.size(); ++i) {
        const path_t &path = spaths[i];
        if (path.sequence.find((char *)skmer1) != std::string::npos)
          good_spaths.push_back(std::make_pair(i, 0));
        else if (path.sequence.find((char *)skmer1_rc) != std::string::npos)
          good_spaths.push_back(std::make_pair(i, 1));
      }

      for (ea_i = 0; ea_i < eanchors.size(); ++ea_i) {
        anchor_t &ea = eanchors[ea_i];
        gbwt::size_type ev1 = ea.v >> 32;
        gbwt::size_type ev2 = (uint32_t)ea.v;
        std::vector<path_t> epaths = graph.get_paths(ev1, ev2, reference_only);
        // get only paths containing the anchor kmer
        char ekmer1[klen];
        char ekmer1_rc[klen];
        memcpy(ekmer1, read + ea.p, klen);
        kmer_d = k2d((char *)ekmer1, klen);
        rckmer_d = rc(kmer_d, klen);
        d2s(rckmer_d, klen, (char *)ekmer1_rc);

        std::vector<std::pair<size_t, bool>> good_epaths;
        for (size_t i = 0; i < epaths.size(); ++i) {
          const path_t &path = epaths[i];
          if (path.sequence.find((char *)ekmer1) != std::string::npos)
            good_epaths.push_back(std::make_pair(i, 0));
          else if (path.sequence.find((char *)ekmer1_rc) != std::string::npos)
            good_epaths.push_back(std::make_pair(i, 1));
        }

        for (const auto &[si, sstrand] : good_spaths) {
          for (const auto &[ei, estrand] : good_epaths) {
            if (spaths[si].id == epaths[ei].id && sstrand == estrand) {
              std::vector<std::pair<gbwt::size_type, gbwt::size_type>> offsets;
              offsets.push_back(std::make_pair(sv1, spaths[si].offset1));
              offsets.push_back(std::make_pair(sv2, spaths[si].offset2));
              offsets.push_back(std::make_pair(ev1, epaths[ei].offset1));
              offsets.push_back(std::make_pair(ev2, epaths[ei].offset2));

              // XXX: filter paths here if anchors/strands do not agree. how?

              //   std::cerr << sstrand << " "
              //             << graph.get_path_contig(spaths[si].id >> 1)
              //             << std::endl;
              //   for (const auto &x : offsets)
              //     std::cerr << graph.get_gfa_name(x.first) << " " << x.second
              //               << std::endl;
              //   std::cerr << std::endl;

              std::sort(
                  offsets.begin(), offsets.end(),
                  [](const std::pair<gbwt::size_type, gbwt::size_type> &a,
                     const std::pair<gbwt::size_type, gbwt::size_type> &b) {
                    return a.second > b.second;
                  });
              // offsets measure the distance to the end of the sequence in
              // nodes
              sa.v = ((uint64_t)offsets[0].first << 32) | offsets[1].first;
              ea.v = ((uint64_t)offsets[2].first << 32) | offsets[3].first;

              goto endloops;
            }
          }
        }
      }
    }

  endloops:
    if (sa_i == sanchors.size()) {
      s.flag |= 2; // flag as invalid
      continue;
    }
    assert(sanchors[sa_i].p < eanchors[ea_i].p);

    // Assigning the anchors
    s.s = sanchors[sa_i].p;
    s.l = eanchors[ea_i].p + klen - sanchors[sa_i].p;
    s.sv1 = sanchors[sa_i].v >> 32;
    s.sv2 = (uint32_t)sanchors[sa_i].v;
    s.ev1 = eanchors[ea_i].v >> 32;
    s.ev2 = (uint32_t)eanchors[ea_i].v;
    // XXX: kmers do not follow "strand" induced by selected path
    s.skmer = std::min(sanchors[sa_i].seq, eanchors[ea_i].seq);
    s.ekmer = std::max(sanchors[sa_i].seq, eanchors[ea_i].seq);
    //
    s.plain_seq = std::string(read + s.s, s.l);
  }
  free(kmer);
}

void remove_duplicates(std::vector<sfs_t> &S) {
  // XXX: what about strings that are prefix/suffix or overlapping?
  for (size_t i = 0; i < S.size(); ++i) {
    if (S[i].flag != 0)
      continue;
    for (size_t j = i + 1; j < S.size(); ++j) {
      if (S[j].flag != 0)
        continue;
      if (S[i].s == S[j].s && S[i].l == S[j].l)
        S[j].flag = 3;
    }
  }
}

int main_anchor(int argc, char *argv[]) {
  double rt;
  int bsize = 10000; // batch size
  int NA = 20;       // number of kmers to check for anchoring
  int nth = 4;       // number of threads
  bool overlapping = true;
  bool reference_only = true;

  int _c;
  while ((_c = getopt(argc, argv, "n:b:@:oah")) != -1) {
    switch (_c) {
    case 'n':
      NA = std::stoi(optarg);
      break;
    case 'b':
      bsize = std::stoi(optarg);
      break;
    case 'o':
      overlapping = false;
      break;
    case 'a':
      reference_only = false;
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
  std::string fx_fn = argv[optind++];
  std::string sfs_fn = argv[optind++];

  // Graph
  Graph graph(gbz_fn, "CHM13");
  graph.load();
  graph.load_fl();

  // Sketch
  rt = realtime();
  sketch_t *sketch = sk_load(skt_fn);
  uint klen = sketch->k;
  fprintf(stderr, "[M::%s] Restored sketch in %.3f sec\n", __func__,
          realtime() - rt);

  std::map<std::string, std::vector<sfs_t>> specific_strings = load_sfs(sfs_fn);
  fprintf(stderr, "[M::%s] Loaded specific strings in %.3f sec\n", __func__,
          realtime() - rt);

  rbatch_t *rb = rbx_init(fx_fn.c_str(), bsize);

  std::vector<int> info(4, 0);
  int total = 0;
  int nreads = 0;
  int x;
  rt = realtime();
  while ((x = rbx_load(rb)) > 0) {
    nreads += x;
#pragma omp parallel for num_threads(nth) schedule(static, 1)
    for (int qq = 0; qq < x; ++qq) {
      // fprintf(stderr, "%s\n", rb->reads[qq]->name);
      char *read_seq = rb->reads[qq]->seq;
      int read_l = rb->reads[qq]->seq_l;
      std::map<uint64_t, int> kmers = count_kmers(read_seq, read_l, klen);
      std::vector<sfs_t> &sfs = specific_strings[rb->reads[qq]->name];
      anchor(graph, sketch, sfs, read_seq, read_l, kmers, klen, NA,
             reference_only, overlapping);
      remove_duplicates(sfs);
    }

    for (int qq = 0; qq < x; ++qq) {
      std::vector<sfs_t> &sfs = specific_strings[rb->reads[qq]->name];
      for (const sfs_t &s : sfs) {
        ++total;
        ++info[s.flag];

        if (s.flag == 3)
          continue;

        printf("%d\t%s\t%d\t%d\t%d\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%s\t%s\t%s\t%"
               "s\t%"
               "s\n",
               s.flag, s.rname.c_str(), s.s, s.l, s.s + s.l, s.sv1, s.sv2,
               s.ev1, s.ev2, s.skmer, s.ekmer,
               graph.get_gfa_name(s.sv1).c_str(),
               graph.get_gfa_name(s.sv2).c_str(),
               graph.get_gfa_name(s.ev1).c_str(),
               graph.get_gfa_name(s.ev2).c_str(),
               s.flag == 0 ? s.plain_seq.c_str() : "*");
      }
    }
    fprintf(stderr,
            "[M::%s] anchored %d specific strings (%d/%d/%d/%d) from %d reads "
            "in %.3f sec\n",
            __func__, total, info[0], info[1], info[2], info[3], nreads,
            realtime() - rt);
  }

  rbx_destroy(rb);
  sk_destroy(sketch);

  return 0;
}