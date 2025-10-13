#include <fstream>
#include <getopt.h>
#include <omp.h>
#include <string>
#include <vector>
#include <zlib.h>

#include "gbwt/fast_locate.h"
#include "gbwt/gbwt.h"
#include "gbwtgraph/gbz.h"
#include "sdsl/simple_sds.hpp"

extern "C" {
#include "fm-index.h"
#include "kseq.h"
}

#include "kmer.hpp"
#include "misc.hpp"
#include "sketch.hpp"
#include "usage.hpp"

KSEQ_INIT(gzFile, gzread)

typedef struct {
  std::string idx;
  std::string seq;
} read_t;

typedef struct {
  uint64_t v;   // vertex on graph
  int p;        // position on query
  uint64_t seq; // kmer
} anchor_t;

typedef struct {
  int s; // start on query
  int l; // length
  uint8_t flag;
  //
  uint64_t sv;
  uint64_t ev;
  uint64_t skmer;
  uint64_t ekmer;
  //
  std::string rname; // plain read name
  uint8_t *seq;      // 1-4encoded sequence
} sfs_t;

// Compute SFS strings from P and store them into solutions
std::vector<sfs_t> ping_pong_search(const rb3_fmi_t *index, uint8_t *P, int l) {
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

    out.push_back({begin, end - begin + 1, 0});

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

std::map<uint64_t, int> count_kmers(uint8_t *read, int readl, int klen) {
  int c; // current char
  uint8_t *kmer = (uint8_t *)malloc(klen);
  std::map<uint64_t, int> kcounts;
  memcpy(kmer, read, klen);
  uint64_t kmer_d = k2d((char *)kmer, klen);
  uint64_t rckmer_d = rc(kmer_d, klen);
  uint64_t ckmer_d = std::min(kmer_d, rckmer_d);
  ++kcounts[ckmer_d];
  for (int p = klen; p < readl; ++p) {
    c = read[p] < 5 ? read[p] - 1 : rand() % 4;
    kmer_d = lsappend(kmer_d, c, klen);
    rckmer_d = rsprepend(rckmer_d, reverse_char(c), klen);
    ckmer_d = std::min(kmer_d, rckmer_d);
    ++kcounts[ckmer_d];
  }
  return kcounts;
}

// Anchor specific strings on graph using graph sketch
void anchor(const gbwtgraph::GBZ &gbz, const gbwt::FastLocate &fl,
            sketch_t *sketch, std::vector<sfs_t> &sfs, uint8_t *read, int readl,
            std::map<uint64_t, int> kcounts, int klen, uint NA,
            bool overlapping) {

  int beg, end;
  uint8_t *kmer = (uint8_t *)malloc(klen);
  uint64_t kmer_d;       // kmer
  uint64_t rckmer_d = 0; // reverse and complemented kmer
  uint64_t ckmer_d = 0;  // canonical kmer
  int c;                 // current char
  uint32_t hit;          // hit from sketch

  for (uint sidx = 0; sidx < sfs.size(); ++sidx) {
    sfs_t &s = sfs[sidx];

    // Finding anchors on the left
    beg = s.s - klen + 1;
    beg = beg < 0 ? 0 : beg;
    memcpy(kmer, read + beg, klen);
    kmer_d = k2d((char *)kmer, klen);
    rckmer_d = rc(kmer_d, klen);
    ckmer_d = std::min(kmer_d, rckmer_d);
    std::vector<anchor_t> sanchors;
    while (beg > 0 && sanchors.size() < NA) {
      if (overlapping) {
        if ((hit = sk_get(sketch, ckmer_d)) != -1U && kcounts[ckmer_d] == 1)
          sanchors.push_back({hit, beg, ckmer_d});
        --beg;
        c = read[beg] < 5 ? read[beg] - 1 : rand() % 4;
        kmer_d = rsprepend(kmer_d, c, klen);
        rckmer_d = lsappend(rckmer_d, reverse_char(c), klen);
        ckmer_d = std::min(kmer_d, rckmer_d);
      } else {
        if ((hit = sk_get(sketch, ckmer_d)) != -1U && kcounts[ckmer_d] == 1) {
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
          c = read[beg] < 5 ? read[beg] - 1 : rand() % 4;
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
        if ((hit = sk_get(sketch, ckmer_d)) != -1U && kcounts[ckmer_d] == 1)
          eanchors.push_back({hit, end, ckmer_d});
        ++end;
        c = read[end + klen - 1] < 5 ? read[end + klen - 1] - 1 : rand() % 4;
        kmer_d = lsappend(kmer_d, c, klen);
        rckmer_d = rsprepend(rckmer_d, reverse_char(c), klen);
        ckmer_d = std::min(kmer_d, rckmer_d);
      } else {
        if ((hit = sk_get(sketch, ckmer_d)) != -1U && kcounts[ckmer_d] == 1) {
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
          c = read[end + klen - 1] < 5 ? read[end + klen - 1] - 1 : rand() % 4;
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

    if (false) {
      // std::cout << read_name << ",";
      for (size_t x = 0; x < NA; ++x) {
        if (x >= sanchors.size()) {
          std::cout << ".,";
          continue;
        }

        const auto &sa = sanchors[x];
        assert(sa.v != -1U);
        if (sa.v == -1U) {
          std::cout << ".,";
          continue;
        }

        gbwtgraph::handle_t handle =
            gbwtgraph::GBWTGraph::node_to_handle(sa.v << 1);
        std::string gfa_idx = gbz.graph.get_segment_name(handle);

        std::vector<gbwt::size_type> paths = fl.decompressSA(sa.v << 1);

        for (const auto &p : paths) {
          gbwt::size_type x = fl.seqId(p);
          // gbwt::size_type xx = fl.seqOffset(p);
          // std::string sample_name =
          //         gbz.index.metadata.fullPath(x).sample_name;
          // std::cerr << p << " " << x << " " << gbwt::Path::id(x) << " " << xx
          // << " " << sample_name << std::endl;
          std::cout << gbwt::Path::id(x)
                    << (gbwt::Path::is_reverse(x) ? "-" : "+") << "/";
        }
        std::cout << ",";
      }

      for (size_t x = 0; x < NA; ++x) {
        if (x >= eanchors.size()) {
          std::cout << ".,";
          continue;
        }

        const auto &ea = eanchors[x];
        assert(ea.v != -1U);
        if (ea.v == -1U) {
          std::cout << ".,";
          continue;
        }
        std::vector<gbwt::size_type> paths = fl.decompressSA(ea.v << 1);
        for (const auto &p : paths) {
          gbwt::size_type x = fl.seqId(p);
          std::cout << gbwt::Path::id(x)
                    << (gbwt::Path::is_reverse(x) ? "-" : "+") << "/";
        }
        std::cout << ",";
      }
      std::cout << std::endl;
    }

    // Finding best pair of anchors
    size_t sa, ea;
    for (sa = 0; sa < sanchors.size(); ++sa) {
      std::vector<gbwt::size_type> soffsets =
          fl.decompressSA(sanchors[sa].v << 1);
      std::set<gbwt::size_type> spaths;
      for (const auto &off : soffsets)
        spaths.insert(fl.seqId(off));

      for (ea = 0; ea < eanchors.size(); ++ea) {
        std::vector<gbwt::size_type> eoffsets =
            fl.decompressSA(eanchors[ea].v << 1);
        std::set<gbwt::size_type> epaths;
        for (const auto &off : eoffsets)
          epaths.insert(fl.seqId(off));

        std::set<gbwt::size_type> res;
        set_intersection(spaths.begin(), spaths.end(), epaths.begin(),
                         epaths.end(), std::inserter(res, res.begin()));
        if (!res.empty())
          break;
      }
      if (ea < eanchors.size())
        break;
    }

    if (sa == sanchors.size()) {
      s.flag |= 2; // flag as invalid
      continue;
    }

    assert(sanchors[sa].p < eanchors[ea].p);

    // Assigning the anchors
    s.s = sanchors[sa].p;
    s.l = eanchors[ea].p + klen - sanchors[sa].p;
    s.sv = sanchors[sa].v;
    s.ev = eanchors[ea].v;
    s.skmer = sanchors[sa].seq;
    s.ekmer = eanchors[ea].seq;
  }
  free(kmer);
}

int load_batch(kseq_t *seq, std::vector<read_t> &entries, int nb) {
  int l = 0;
  int n = 0;
  while (n < nb && (l = kseq_read(seq)) >= 0) {
    entries[n++] = {seq->name.s, seq->seq.s};
  }
  return n;
}

int main_search(int argc, char *argv[]) {
  double rt;
  int bsize = 10000; // batch size
  int NA = 20;       // number of kmers to check for anchoring
  int nth = 4;       // number of threads
  bool overlapping = true;

  int _c;
  while ((_c = getopt(argc, argv, "b:@:oh")) != -1) {
    switch (_c) {
    case 'b':
      bsize = std::stoi(optarg);
      break;
    case 'o':
      overlapping = false;
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

  rt = realtime();
  gbwtgraph::GBZ gbz;

  sdsl::simple_sds::load_from(gbz, gbz_fn);
  fprintf(stderr, "[M::%s] Loaded GBZ in %.3fs\n", __func__, realtime() - rt);

  // R-index
  rt = realtime();
  gbwt::FastLocate fl;
  std::ifstream in;
  in.open(gbz_fn + ".ri", std::ifstream::in);
  fl.load(in);
  fl.setGBWT(gbz.index);
  in.close();
  fprintf(stderr, "[M::%s] Restored R-index in %.3f sec\n", __func__,
          realtime() - rt);

  rt = realtime();
  sketch_t *sketch = sk_load(skt_fn);
  uint klen = sketch->k;
  fprintf(stderr, "[M::%s] Restored sketch in %.3f sec\n", __func__,
          realtime() - rt);

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

  std::vector<read_t> entries(bsize);
  std::vector<std::vector<sfs_t>> output(bsize);

  gzFile fp = gzopen(fx_fn.c_str(), "r");
  kseq_t *seq = kseq_init(fp);

  int x;
  while ((x = load_batch(seq, entries, bsize)) > 0) {
    rt = realtime();
    fprintf(stderr, "[M::%s] loaded %d reads in %.3f sec\n", __func__, x,
            realtime() - rt);
    rt = realtime();
#pragma omp parallel for num_threads(nth) schedule(static, 1)
    for (int qq = 0; qq < x; ++qq) {
      uint8_t *seq = (uint8_t *)entries[qq].seq.data();
      int l = entries[qq].seq.size();
      rb3_char2nt6(l, seq);
      output[qq] = ping_pong_search(&fmd, seq, l);
      assemble(output[qq]);
      std::map<uint64_t, int> kmers = count_kmers(seq, l, klen);
      anchor(gbz, fl, sketch, output[qq], seq, l, kmers, klen, NA, overlapping);

      for (sfs_t &s : output[qq]) {
        s.rname = entries[qq].idx;

        s.seq = (uint8_t *)malloc(s.l + 1);
        for (int i = 0; i < s.l; ++i)
          s.seq[i] = seq[s.s + i];
        s.seq[s.l] = '\0';
      }
    }

    int total = 0;
    std::vector<int> info(3, 0);
    for (int qq = 0; qq < x; ++qq) {
      for (uint j = 0; j < output[qq].size(); ++j) {
        const sfs_t &s = output[qq][j];
        ++total;
        ++info[s.flag];
        gbwtgraph::handle_t sh =
            gbwtgraph::GBWTGraph::node_to_handle(s.sv << 1);
        gbwtgraph::handle_t eh =
            gbwtgraph::GBWTGraph::node_to_handle(s.ev << 1);
        std::string s_gfa_idx = gbz.graph.get_segment_name(sh);
        std::string e_gfa_idx = gbz.graph.get_segment_name(eh);

        char sseq[s.l];
        for (int i = 0; i < s.l; ++i)
          sseq[i] = "NACGTN"[(int)s.seq[i]];
        sseq[s.l] = '\0';

        std::cout << (int)s.flag << "\t" << s.rname << "\t" << s.s << "\t"
                  << s.l << "\t" << s.s + s.l << "\t" << s.sv << "\t" << s.ev
                  << "\t" << s.skmer << "\t" << s.ekmer << "\t" << s_gfa_idx
                  << ">" << e_gfa_idx << "\t" << sseq << std::endl;
      }
    }
    fprintf(stderr,
            "[M::%s] computed %d (%d/%d/%d) specific strings in %.3f sec\n",
            __func__, total, info[0], info[1], info[2], realtime() - rt);
  }

  kseq_destroy(seq);
  gzclose(fp);
  sk_destroy(sketch);
  rb3_fmi_free(&fmd);

  return 0;
}
