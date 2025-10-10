#include <fstream>
#include <getopt.h>
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
  uint64_t v;   // vertex on graph
                //   int offset;   // offset on vertex
  int p;        // position on query
  uint64_t seq; // kmer
  //   int closest;  // is this anchor the closest one to the specific string?
} anchor_t;

typedef struct {
  int qidx; // read index
  int s;    // start on query
  int l;    // length
  //   anchor_t a;   // left anchor
  //   anchor_t b;   // right anchor
  //   int strand;   // inferred strand
  //   uint64_t esk; // expected starting kmer (from cluster)
  //   uint64_t eek; // expected ending kmer (from cluster)
  //   int good;     // is it good for calling step?
  //   char *rname;  // plain read name
  //   uint8_t *seq; // 1-4encoded sequence
} sfs_t;

// Compute SFS strings from P and store them into solutions
std::vector<sfs_t> ping_pong_search(const rb3_fmi_t *index, uint8_t *P, int l,
                                    int qidx) {
  std::vector<sfs_t> out;
  rb3_sai_t ik;
  int begin = l - 1;
  while (begin >= 0) {
    // Backward search. Stop at first mismatch.
    int bmatches = 0;
    rb3_fmd_set_intv(index, P[begin], &ik);
    while (ik.size != 0 && begin > 0) {
      --begin;
      ++bmatches;
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
    int fmatches = 0;
    rb3_fmd_set_intv(index, P[end], &ik);
    while (ik.size != 0) {
      ++end;
      ++fmatches;
      rb3_sai_t ok[RB3_ASIZE];
      rb3_fmd_extend(index, &ik, ok, 0);
      ik = ok[P[end] >= 1 && P[end] <= 4 ? 5 - P[end] : P[end]];
    }

    out.push_back({qidx, begin, end - begin + 1});

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
      S[new_n].qidx = S[i].qidx;
      S[new_n].s = S[i].s;
      S[new_n].l = S[i].l;
      ++new_n;
    }
    ++i;
  }
  S.resize(new_n);
}

// Anchor specific strings on graph using graph sketch
void anchor(const gbwtgraph::GBZ &gbz, gbwt::FastLocate fl, sketch_t *sketch,
            std::vector<sfs_t> &sfs, uint8_t *read, int readl, int klen,
            uint NA) {

  int beg, end;
  uint8_t *kmer = (uint8_t *)malloc(klen);
  uint64_t kmer_d;       // kmer
  uint64_t rckmer_d = 0; // reverse and complemented kmer
  uint64_t ckmer_d = 0;  // canonical kmer
  int c;                 // current char
  uint64_t hit;          // hit from sketch

  for (uint sidx = 0; sidx < sfs.size(); ++sidx) {
    sfs_t &s = sfs[sidx];

    // Finding anchors in flanking regions
    beg = s.s - klen + 1;
    beg = beg < 0 ? 0 : beg;
    end = s.s + s.l - 1;
    end = end > readl - klen ? readl - klen : end;

    memcpy(kmer, read + beg, klen);
    kmer_d = k2d((char *)kmer, klen);
    rckmer_d = rc(kmer_d, klen);
    ckmer_d = std::min(kmer_d, rckmer_d);

    std::vector<anchor_t> sanchors;
    // int nsa = 0;
    while (beg > 0 && sanchors.size() < NA) {
      if ((hit = sk_get(sketch, ckmer_d)) != -1)
        sanchors.push_back({hit, beg, ckmer_d});
      --beg;
      c = read[beg] < 5 ? read[beg] - 1 : rand() % 4;
      kmer_d = rsprepend(kmer_d, c, klen);
      rckmer_d = lsappend(rckmer_d, reverse_char(c), klen);
      ckmer_d = std::min(kmer_d, rckmer_d);
    }

    memcpy(kmer, read + end, klen);
    kmer_d = k2d((char *)kmer, klen);
    rckmer_d = rc(kmer_d, klen);
    ckmer_d = std::min(kmer_d, rckmer_d);
    std::vector<anchor_t> eanchors;
    while (end < readl - klen + 1 && eanchors.size() < NA) {
      if ((hit = sk_get(sketch, ckmer_d)) != -1) {
        eanchors.push_back({hit, end, ckmer_d});
      }
      ++end;
      c = read[end + klen - 1] < 5 ? read[end + klen - 1] - 1 : rand() % 4;
      kmer_d = lsappend(kmer_d, c, klen);
      rckmer_d = rsprepend(rckmer_d, reverse_char(c), klen);
      ckmer_d = std::min(kmer_d, rckmer_d);
    }

    if (sanchors.size() == 0 || eanchors.size() == 0) {
      // s.qidx = -1; // tag as invalid
      continue;
    }

    // keep only unique anchors on read
    for (size_t i1 = 0; i1 < sanchors.size(); ++i1) {
      // Starting anchors vs starting anchors
      if (sanchors[i1].v == -1)
        continue;
      for (size_t i2 = i1 + 1; i2 < sanchors.size(); ++i2) {
        if (sanchors[i2].v == -1)
          continue;
        if (sanchors[i1].seq == sanchors[i2].seq) {
          sanchors[i1].v = -1;
          sanchors[i2].v = -1;
        }
      }
      // Starting anchors vs ending anchors
      for (size_t i2 = 0; i2 < eanchors.size(); ++i2) {
        if (eanchors[i2].v == -1)
          continue;
        if (sanchors[i1].seq == eanchors[i2].seq) {
          sanchors[i1].v = -1;
          eanchors[i2].v = -1;
        }
      }
    }
    // Ending anchors vs ending anchors
    for (size_t i1 = 0; i1 < eanchors.size(); ++i1) {
      if (eanchors[i1].v == -1)
        continue;
      for (size_t i2 = i1 + 1; i2 < eanchors.size(); ++i2) {
        if (eanchors[i2].v == -1)
          continue;
        if (eanchors[i1].seq == eanchors[i2].seq) {
          eanchors[i1].v = -1;
          eanchors[i2].v = -1;
        }
      }
    }

    // Finding best pair of anchors
    //     int mind = 100, dist;
    //     int sax = -1, eax = -1; // index for selected anchors

    //     std::map<std::pair<int, int>, int> memo;
    //     std::map<std::pair<int, int>, int>::iterator hhit;
    //     int x = 0, y = 0;
    //     int xoff = 0, yoff = 0;
    //     std::pair<int, int> xy = {x, y};
    for (size_t i = 0; i < sanchors.size(); ++i) {
      // int sv = sanchors[i].v;
      gbwt::node_type sv = sanchors[i].v;

      std::cerr << sanchors[i].seq << " " << sv << " "
                << gbz.graph.get_segment_name(
                       gbwtgraph::GBWTGraph::node_to_handle(sv << 1))
                << std::endl;

      if (!gbz.graph.has_node(sv)) {
        std::cerr << "not" << std::endl;
      }
      // std::set<gbwt::size_type> spaths;

      //   std::vector<gbwt::size_type> xxx = fl.decompressSA(sv);
      //   std::cerr << xxx.size() << std::endl;

      //   for (int strand = 0; strand < 2; ++strand) {
      //     int vv = sv; // gbwt::Node::encode(sv, strand);
      //     if (!gbz.index.contains(vv)) {
      //       std::cerr << "wrong" << std::endl;
      //       continue;
      //     }

      //     std::vector<gbwt::size_type> xxx = fl.decompressSA(vv);
      //     std::cerr << xxx.size() << std::endl;
      //     for (const gbwt::size_type &x : xxx) {
      //       std::cerr << "ok" << std::endl;
      //       gbwt::size_type pidx = fl.seqId(x);
      //       gbwt::size_type pp = gbwt::Path::id(pidx);
      //       std::cerr << pidx << " " << pp << std::endl;
      //     }
      // }

      //       if (x == -1)
      //         // anchor has been filtered out since it was repeated in the
      //         read continue;
      //       xoff = sanchors[i].offset;
      //       xy.first = x;
      //       for (int j = 0; j < nea; ++j) {
      //         y = eanchors[j].v;
      //         if (y == -1)
      //           // anchor has been filtered out since it was repeated in the
      //           read continue;
      //         yoff = eanchors[j].offset;
      //         xy.second = y;
      //         if ((hhit = memo.find(xy)) == memo.end()) {
      //           memo[xy] = graph.distance(sanchors[i].v, eanchors[j].v);
      //           memo[std::make_pair(y, x)] = memo[xy];
      //         }
      //         dist = memo[xy];
      //         if (dist < 0)
      //           continue;
      //         if (x == y && (xoff == yoff || (xoff < yoff && xoff + klen >=
      //         yoff)
      //         ||
      //                        (xoff > yoff && yoff + klen >= xoff)))
      //           continue;
      //         if (dist < mind) {
      //           sax = i;
      //           eax = j;
      //           mind = dist;
      //         }
      //       }
    }
    //     if (sax == -1 || eax == -1) {
    //       s.qidx = -1; // tag as invalid
    //       continue;
    //     }

    // Assigning the anchors
    //     anchor_t sa = sanchors[sax];
    //     anchor_t ea = eanchors[eax];
    //     int b = sa.p;
    //     int l = ea.p + klen - sa.p;
    //     int strand = 1;
    //     if (sa.v > ea.v || (sa.v == ea.v && sa.offset > ea.offset)) {
    //       anchor_t tmp = sa;
    //       sa = ea;
    //       ea = tmp;
    //       strand = 0;
    //     }

    //     s.s = b;
    //     s.l = l;
    //     s.a = sa;
    //     s.a.closest = sax;
    //     s.b = ea;
    //     s.b.closest = eax;
    //     s.strand = strand;

    //     assert(sa.v <= ea.v);
  }
  free(kmer);
}

int main_search(int argc, char *argv[]) {
  double rt;
  int NA = 20;
  // int nth = 1; // number of threads

  int _c;
  while ((_c = getopt(argc, argv, "@:h")) != -1) {
    switch (_c) {
    // case '@':
    //   nth = std::stoi(optarg);
    //   break;
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

  rt = realtime();
  gbwt::FastLocate fl;
  std::ifstream in;
  in.open(gbz_fn + ".ri", std::ifstream::in | std::ifstream::app);
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

  gzFile fp = gzopen(fx_fn.c_str(), "r");
  kseq_t *seq = kseq_init(fp);
  int l;
  int qidx = 0;
  while ((l = kseq_read(seq)) >= 0) {
    uint8_t *useq = (uint8_t *)seq->seq.s;
    rb3_char2nt6(l, useq);
    std::vector<sfs_t> sfs = ping_pong_search(&fmd, useq, l, qidx);
    // std::cerr << sfs.size() << std::endl;
    assemble(sfs);
    // std::cerr << sfs.size() << std::endl;

    anchor(gbz, fl, sketch, sfs, useq, seq->seq.l, klen, NA);

    // // strand stuff
    // int strands[2] = {0, 0};
    // for (uint ss = 0; ss < output[qq].size(); ++ss) {
    //   if (output[qq][ss].qidx != -1)
    //     ++strands[output[qq][ss].strand];
    // }
    // // + strand if tie
    // int strand = strands[0] > strands[1] ? 0 : 1;
    // for (uint ss = 0; ss < output[qq].size(); ++ss) {
    //   if (output[qq][ss].qidx == -1)
    //     continue;
    //   if (output[qq][ss].strand == 0)
    //     // reverse
    //     output[qq][ss].s = l - (output[qq][ss].s + output[qq][ss].l);
    //   if (output[qq][ss].strand != strand)
    //     output[qq][ss].strand = 2;
    // }

    // // Finally reassemble good ones
    // assemble2(output[qq], strand);

    ++qidx;
  }

  kseq_destroy(seq);
  gzclose(fp);
  sk_destroy(sketch);
  rb3_fmi_free(&fmd);

  return 0;
}
