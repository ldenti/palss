#include <assert.h>
#include <cstdint>
#include <iostream>
#include <map>
#include <vector>
#include <zlib.h>

#include "fm-index.h"
#include "ketopt.h"
#include "kseq.h"

#include "sketch.hpp"
extern "C" {
#include "graph.h"
#include "utils.h"
}

KSTREAM_INIT(gzFile, gzread, 65536)

using namespace std;

/* Backward search */
int search(const rb3_fmi_t *index, uint8_t *kmer, int k) {
  rb3_sai_t ik;
  int begin = k - 1;
  rb3_fmd_set_intv(index, kmer[begin], &ik);
  while (ik.size != 0 && begin > 0) {
    --begin;
    rb3_sai_t ok[RB3_ASIZE]; // output SA intervals (one for each symbol
                             // between 0 and 5)
    rb3_fmd_extend(index, &ik, ok, 1);
    ik = ok[kmer[begin]];
  }
  return ik.size;
}

/* Sketch a set of segments */
void run_sketching(seg_t **segs, int ns, uint8_t klen, rb3_fmi_t *fmd, int ng,
                   int threads, map<uint64_t, uint64_t> &sketch) {
  double rt = realtime();
  vector<map<uint64_t, uint64_t>> sketches(ns);

#pragma omp parallel for num_threads(threads)
  for (int i = 0; i < ns; ++i) {
    char *kmer = (char *)malloc(sizeof(char) *
                                (klen + 1)); // first kmer on sequence (plain)
    uint8_t *s;
    uint64_t kmer_d = 0;   // kmer
    uint64_t rckmer_d = 0; // reverse and complemented kmer
    uint64_t ckmer_d = 0;  // canonical kmer
    uint8_t c;             // new character to append
    int p = 0;             // current position on segment
    int hits = 0;          // hits in the FMD index
    seg_t *seg = segs[i];
    strncpy(kmer, seg->seq, klen);
    kmer_d = k2d(kmer, klen);
    rckmer_d = rc(kmer_d, klen);
    ckmer_d = std::min(kmer_d, rckmer_d);

    s = (uint8_t *)kmer;
    rb3_char2nt6(klen, s);
    hits = search(fmd, s, klen);
    assert(hits > 0);
    sk_add(sketches[i], ckmer_d, seg->idx, 0, hits == ng);

    for (p = klen; p < seg->l; ++p) {
      c = to_int[seg->seq[p]] - 1; // A is 1 but it should be 0
      kmer_d = lsappend(kmer_d, c, klen);
      rckmer_d = rsprepend(rckmer_d, reverse_char(c), klen);
      ckmer_d = std::min(kmer_d, rckmer_d);

      d23(kmer_d, klen, kmer);
      s = (uint8_t *)kmer;
      rb3_char2nt6(klen, s);
      hits = search(fmd, s, klen);
      assert(hits >= 0);
      sk_add(sketches[i], ckmer_d, seg->idx, p - klen + 1, hits == ng);
    }
    free(kmer);
  }
  uint64_t nkmers = 0;
  for (int i = 0; i < ns; ++i) {
    for (auto &pair : sketches[i]) {
      ++nkmers;
      auto x = sketch.find(pair.first);
      if (x == sketch.end()) {
        sketch[pair.first] = pair.second;
      } else {
        pair.second = pair.second & ~1;
      }
    }
  }

  fprintf(stderr, "[M::%s] sketched %d segments (%ld kmers) in %.3f sec\n",
          __func__, ns, nkmers, realtime() - rt);
}

int main_sketch(int argc, char *argv[]) {
  double rt0 = realtime();
  double rt = rt0;

  int klen = 27;  // kmer size
  int ng = 1;     // expected number of genomes
  int vpb = 5000; // number of vertices to load per batch
  int threads = 4;
  static ko_longopt_t longopts[] = {{NULL, 0, 0}};
  ketopt_t opt = KETOPT_INIT;
  int _c;
  while ((_c = ketopt(&opt, argc, argv, 1, "k:g:v:@:", longopts)) >= 0) {
    if (_c == 'k')
      klen = atoi(opt.arg);
    else if (_c == 'g')
      ng = atoi(opt.arg);
    else if (_c == 'v')
      vpb = atoi(opt.arg);
    else if (_c == '@')
      threads = atoi(opt.arg);
  }

  if (argc - opt.ind != 2) {
    fprintf(stderr, "Argh");
    return 1;
  }
  char *gfa_fn = argv[opt.ind++];
  char *fmd_fn = argv[opt.ind++];

  // FMD-index loading
  rb3_fmi_t fmd;
  rb3_fmi_restore(&fmd, fmd_fn, 0);
  if (fmd.e == 0 && fmd.r == 0) {
    fprintf(stderr, "Error restoring index");
    return 1;
  }
  fprintf(stderr, "[M::%s] restored FMD index in %.3f sec\n", __func__,
          realtime() - rt);
  rt = realtime();
  // ---

  // Sketching
  kstring_t s = {0, 0, 0};
  int dret;
  gzFile fp = gzopen(gfa_fn, "r");
  if (fp == 0)
    return 0;
  kstream_t *ks = ks_init(fp);

  sketch_t sketch;
  seg_t **segs = (seg_t **)malloc(vpb * sizeof(seg_t *));
  for (int i = 0; i < vpb; ++i)
    segs[i] = init_seg();
  int si = 0; // current segment
  int nvertices = 0;
  while (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) >= 0) {
    if (s.s[0] == 'S') {
      ++nvertices;
      gfa_parse_S(s.s, segs[si]);
      if (segs[si]->l < klen)
        continue;
      ++si;

      if (si == vpb) {
        run_sketching(segs, si, klen, &fmd, ng, threads, sketch);
        si = 0;
      }
    }
  }
  if (si < vpb) {
    run_sketching(segs, si, klen, &fmd, ng, threads, sketch);
  }

  // Clean everything
  for (int i = 0; i < vpb; ++i)
    destroy_seg(segs[i]);
  free(segs);
  free(s.s);
  ks_destroy(ks);
  gzclose(fp);
  rb3_fmi_free(&fmd);

  fprintf(stderr, "[M::%s] sketched graph in %.3f sec\n", __func__,
          realtime() - rt);
  rt = realtime();
  // ---

  // Write sketch to stdout
  sk_store(sketch, "-");
  // ---

  fprintf(stderr, "[M::%s] completed in %.3f sec\n", __func__,
          realtime() - rt0);

  return 0;
}
