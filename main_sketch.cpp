#include <algorithm>
#include <assert.h>
#include <getopt.h>
#include <omp.h>
#include <set>
#include <stdio.h>
#include <zlib.h>

extern "C" {
#include "fm-index.h"
#include "kseq.h"
}

#include "graph.hpp"
#include "kmer.hpp"
#include "misc.hpp"
#include "sketch.hpp"
#include "usage.h"

KSTREAM_INIT(gzFile, gzread, 65536)

// Backward search until qinterval size is > than min_size
int search(const rb3_fmi_t *fmd, rb3_sai_t ik, uint8_t *kmer, int k,
           int min_size) {
  int begin = k;
  while (ik.size != 0 && begin > 0) {
    --begin;
    rb3_sai_t ok[RB3_ASIZE]; // output SA intervals (one for each symbol
                             // between 0 and 5)
    rb3_fmd_extend(fmd, &ik, ok, 1);
    ik = ok[kmer[begin]];
    if (ik.size <= min_size)
      return ik.size;
  }
  return ik.size;
}

// Backward search a given kmer and set qinterval ik
void set_qint(const rb3_fmi_t *fmd, rb3_sai_t *ik, uint8_t *kmer, int k) {
  int begin;
  begin = k - 1;
  rb3_fmd_set_intv(fmd, kmer[begin], ik);

  while (ik->size != 0 && begin > 0) {
    --begin;
    rb3_sai_t ok[RB3_ASIZE]; // output SA intervals (one for each symbol
                             // between 0 and 5)
    rb3_fmd_extend(fmd, ik, ok, 1);
    memcpy(ik, &ok[kmer[begin]], sizeof(rb3_sai_t));
  }
}

// Q-intervals memoization (all mklen-mers)
rb3_sai_t *memoize(const rb3_fmi_t *fmd, const int mklen) {
  int n = 1 << (mklen * 2);
  rb3_sai_t *qints = (rb3_sai_t *)malloc(n * sizeof(rb3_sai_t));

  char *kmer = (char *)malloc(mklen + 1);
  uint8_t *qkmer;
  kmer[mklen] = '\0';
  for (int i = 0; i < n; ++i) {
    d23(i, mklen, kmer);
    qkmer = (uint8_t *)kmer;
    rb3_char2nt6(mklen, qkmer);
    set_qint(fmd, &qints[i], qkmer, mklen);
  }
  free(kmer);
  return qints;
}

int main_sketch(int argc, char *argv[]) {
  double rt0 = realtime();
  double rt = rt0, rt1;

  int klen = 27;      // kmer size
  int mklen = 9;      // memoization kmer size
  int nh = INT32_MAX; // expected number of haplotypes
  bool txt_flag = 0;  // dump sketch in txt
  bool big_flag = 0;  // do we need big sketch?
  int nth = 4;        // number of threads
  std::string path_prefix = "CHM13";
  std::string out_fn = "-";

  int _c;
  while ((_c = getopt(argc, argv, "k:m:g:@:p:tbh")) != -1) {
    switch (_c) {
    case 'k':
      klen = std::stoi(optarg);
      break;
    case 'm':
      mklen = std::stoi(optarg);
      break;
    case 'g':
      nh = std::stoi(optarg);
      break;
    case 't':
      txt_flag = true;
      break;
    case 'b':
      big_flag = true;
      break;
    case 'p':
      path_prefix = optarg;
      break;
    case '@':
      nth = std::stoi(optarg);
      break;
    case 'h':
      fprintf(stderr, "%s", SKETCH_USAGE_MESSAGE);
      return 0;
    }
  }

  if (argc - optind != 2) {
    fprintf(stderr, "%s", SKETCH_USAGE_MESSAGE);
    return 1;
  }
  char *gfa_fn = argv[optind++];
  char *fmd_fn = argv[optind++];

  // FMD-index loading
  rb3_fmi_t fmd;
  rb3_fmi_restore(&fmd, fmd_fn, 0);
  if (fmd.e == 0 && fmd.r == 0) {
    fprintf(stderr, "[E::%s] FMD-index cannot be restored", __func__);
    exit(EXIT_FAILURE);
  }
  fprintf(stderr, "[M::%s] restored FMD index in %.3f sec\n", __func__,
          realtime() - rt);
  rt = realtime();

  // Q-intervals memoization
  assert(mklen < 16);
  int memo_bm = (1 << (mklen * 2)) - 1;
  rb3_sai_t *qints = memoize(&fmd, mklen);
  fprintf(stderr, "[M::%s] memoization (all %d-mers) in %.3f sec\n", __func__,
          mklen, realtime() - rt);
  rt = realtime();

  // First pass over graph. Store all kmers in an array
  sketch_t *sketch = sk_init((uint64_t)1 << (big_flag ? 32 : 30), klen,
                             mklen); // XXX: hardcoded

  uint64_t *kmers = sketch->vls; // we are going to "reuse" values array
  uint64_t totkmers = 0;         // index for insertion in kmers
  seg_t seg;

  char *kmer = (char *)malloc(sizeof(char) *
                              (klen + 1)); // first kmer on sequence (plain)
  uint64_t kmer_d = 0;                     // kmer
  uint64_t rckmer_d = 0;                   // reverse and complemented kmer
  uint64_t ckmer_d = 0;                    // canonical kmer
  uint8_t c;                               // new character to append
  int p = 0;                               // current position on segment
  int nvertices = 0;                       // number of vertices processed

  kstring_t s = {0, 0, 0};
  int dret;
  gzFile fp = gzopen(gfa_fn, "r");
  if (fp == 0)
    return 0;
  kstream_t *ks = ks_init(fp);
  std::set<uint64_t> ref_vs; // vertices on reference paths (GFA space)
  rt1 = realtime();
  while (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) >= 0) {
    if (s.s[0] == 'S') {
      ++nvertices;
      if (nvertices % 1000000 == 0) {
        fprintf(stderr, "[M::%s] parsed 1M vertices in %.3f sec\n", __func__,
                realtime() - rt1);
        rt1 = realtime();
      }

      gfa_parse_S(s.s, seg);
      if (seg.l < klen)
        continue;

      strncpy(kmer, seg.seq.c_str(), klen);
      // first kmer
      kmer_d = k2d(kmer, klen);
      rckmer_d = rc(kmer_d, klen);
      ckmer_d = std::min(kmer_d, rckmer_d);
      kmers[totkmers++] = ckmer_d;
      // all other kmers
      for (p = klen; p < seg.l; ++p) {
        c = to_int[seg.seq[p]] - 1; // A is 1 but it should be 0
        kmer_d = lsappend(kmer_d, c, klen);
        rckmer_d = rsprepend(rckmer_d, reverse_char(c), klen);
        ckmer_d = std::min(kmer_d, rckmer_d);
        kmers[totkmers++] = ckmer_d;
      }
    } else if (s.s[0] == 'P' || s.s[0] == 'W') {
      path_t path;
      if (s.s[0] == 'P')
        gfa_parse_P(s.s, path, path_prefix);
      else if (s.s[0] == 'W')
        gfa_parse_W(s.s, path, path_prefix);
      for (uint v = 0; v < path.vertices.size(); ++v)
        ref_vs.insert(path.vertices[v]);
    }
  }
  gzclose(fp);
  ks_destroy(ks);
  fprintf(stderr,
          "[M::%s] loaded %ld kmers (from %d vertices, reference vertices: "
          "%ld, %f) "
          "in %.3f sec\n",
          __func__, totkmers, nvertices, ref_vs.size(),
          ref_vs.size() / (float)nvertices, realtime() - rt);

  // Sort kmers
  rt = realtime();
  std::sort(kmers, kmers + totkmers);
  fprintf(stderr, "[M::%s] sorted %ld kmers in %.3f sec\n", __func__, totkmers,
          realtime() - rt);

  // Flag all repeated anchors
  rt = realtime();
  for (uint64_t i = 0; i < totkmers;) {
    kmer_d = kmers[i];
    ++i;
    while (i < totkmers && kmers[i] == kmer_d) {
      kmers[i - 1] = 0;
      kmers[i] = 0;
      ++i;
    }
  }
  fprintf(stderr, "[M::%s] flagged kmers in %.3f sec\n", __func__,
          realtime() - rt);

  fprintf(
      stderr,
      "[M::%s] searching kmers in the FMD-index using %d threads (nh : %d)\n",
      __func__, nth, nh);
  rt = realtime();

  // Flag all non-solid anchors by querying the FMD-index
  char **kmers_t = (char **)malloc(nth * sizeof(char *));
  for (int i = 0; i < nth; ++i)
    kmers_t[i] = (char *)malloc(sizeof(char) * (klen + 1));
  uint8_t **qkmers_t = (uint8_t **)malloc(nth * sizeof(uint8_t *));
#pragma omp parallel for num_threads(nth)
  for (uint64_t i = 0; i < totkmers; ++i) {
    if (kmers[i] == 0)
      continue;
    int tidx = omp_get_thread_num();

    d23(kmers[i], klen, kmers_t[tidx]);
    qkmers_t[tidx] = (uint8_t *)kmers_t[tidx];
    rb3_char2nt6(klen, qkmers_t[tidx]);

    rb3_sai_t qint = qints[kmers[i] & memo_bm]; // memoization interval
    int hits = 1;                               // hits in the FMD index
    // XXX: assuming nh>0
    if (qint.size > nh) {
      hits = search(&fmd, qint, qkmers_t[tidx], klen - mklen, nh);
      assert(hits > 0);
    }
    if (hits > nh)
      kmers[i] = 0;
  }
  free(qkmers_t);
  for (int i = 0; i < nth; ++i)
    free(kmers_t[i]);
  free(kmers_t);
  fprintf(stderr, "[M::%s] searched kmers in %.3f sec\n", __func__,
          realtime() - rt);
  rt = realtime();

  // Add solid anchors to sketch
  for (uint64_t i = 0; i < totkmers; ++i) {
    kmer_d = kmers[i];
    if (kmer_d == 0)
      continue;
    sk_add(sketch, kmer_d);
  }
  fprintf(stderr, "[M::%s] added %ld kmers to sketch in %.3f sec\n", __func__,
          sketch->n, realtime() - rt);

  // Reiterate over graph to assign values to solid anchors
  fp = gzopen(gfa_fn, "r");
  ks = ks_init(fp);

  rt1 = realtime();
  nvertices = 0;
  while (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) >= 0) {
    if (s.s[0] == 'S') {
      ++nvertices;
      if (nvertices % 1000000 == 0) {
        fprintf(stderr, "[M::%s] parsed 1M vertices in %.3f sec\n", __func__,
                realtime() - rt1);
        rt1 = realtime();
      }

      gfa_parse_S(s.s, seg);
      if (seg.l < klen)
        continue;

      // first kmer
      strncpy(kmer, seg.seq.c_str(), klen);
      kmer_d = k2d(kmer, klen);
      rckmer_d = rc(kmer_d, klen);
      ckmer_d = std::min(kmer_d, rckmer_d);
      sk_add_v(sketch, ckmer_d, seg.idx, 0,
               ref_vs.find(seg.idx) != ref_vs.end());
      // all other kmers
      for (p = klen; p < seg.l; ++p) {
        c = to_int[seg.seq[p]] - 1; // A is 1 but it should be 0
        kmer_d = lsappend(kmer_d, c, klen);
        rckmer_d = rsprepend(rckmer_d, reverse_char(c), klen);
        ckmer_d = std::min(kmer_d, rckmer_d);
        sk_add_v(sketch, ckmer_d, seg.idx, p - klen + 1,
                 ref_vs.find(seg.idx) != ref_vs.end());
      }
    }
  }
  gzclose(fp);
  ks_destroy(ks);
  fprintf(stderr, "[M::%s] sketched %ld kmers in %.3f sec\n", __func__,
          sketch->n, realtime() - rt);

  // Write sketch
  rt = realtime();
  if (txt_flag)
    sk_dump(sketch, out_fn.c_str());
  else
    sk_store(sketch, out_fn.c_str());
  fprintf(stderr, "[M::%s] dumped sketch in %.3f sec\n", __func__,
          realtime() - rt);

  // Clean everything
  free(kmer);
  free(s.s);
  rb3_fmi_free(&fmd);
  free(qints);
  sk_destroy(sketch);

  fprintf(stderr, "[M::%s] completed in %.3f sec\n", __func__,
          realtime() - rt0);

  return 0;
}
