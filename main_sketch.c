#include <assert.h>
#include <omp.h>
#include <stdio.h>
#include <zlib.h>

#include "fm-index.h"
#include "ketopt.h"
#include "kseq.h"
#include "ksort.h"

#include "graph.h"
#include "kmer.h"
#include "misc.h"
#include "sketch.h"
#include "usage.h"

KSORT_INIT_GENERIC(uint64_t)

/* Backward search until qinterval size is > than min_size */
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

/* Backward search a given kmer and set qinterval ik */
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

/* Q-intervals memoization (all mklen-mers) */
rb3_sai_t *memoize(const rb3_fmi_t *fmd, const int mklen) {
  int n = 1 << (mklen * 2);
  rb3_sai_t *qints = malloc(n * sizeof(rb3_sai_t));

  char *kmer = malloc(mklen + 1);
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
  int txt_f = 0;      // dump sketch in txt
  int big_f = 0;      // do we need big sketch?
  int nth = 4;        // number of threads
  static ko_longopt_t longopts[] = {{NULL, 0, 0}};
  ketopt_t opt = KETOPT_INIT;
  int _c;
  while ((_c = ketopt(&opt, argc, argv, 1, "k:m:g:@:tbh", longopts)) >= 0) {
    if (_c == 'k')
      klen = atoi(opt.arg);
    else if (_c == 'm')
      mklen = atoi(opt.arg);
    else if (_c == 'g')
      nh = atoi(opt.arg);
    else if (_c == 't')
      txt_f = 1;
    else if (_c == 'b')
      big_f = 1;
    else if (_c == '@')
      nth = atoi(opt.arg);
    else if (_c == 'h') {
      fprintf(stderr, "%s", SKETCH_USAGE_MESSAGE);
      return 0;
    }
  }

  if (argc - opt.ind != 2) {
    fprintf(stderr, "%s", SKETCH_USAGE_MESSAGE);
    return 1;
  }
  char *gfa_fn = argv[opt.ind++];
  char *fmd_fn = argv[opt.ind++];

  /* FMD-index loading */
  rb3_fmi_t fmd;
  rb3_fmi_restore(&fmd, fmd_fn, 0);
  if (fmd.e == 0 && fmd.r == 0) {
    fprintf(stderr, "[E::%s] FMD-index cannot be restored", __func__);
    exit(EXIT_FAILURE);
  }
  fprintf(stderr, "[M::%s] restored FMD index in %.3f sec\n", __func__,
          realtime() - rt);
  rt = realtime();

  /* Q-intervals memoization */
  int memo_bm = (1 << (mklen * 2)) - 1;
  rb3_sai_t *qints = memoize(&fmd, mklen);
  fprintf(stderr, "[M::%s] memoization (all %d-mers) in %.3f sec\n", __func__,
          mklen, realtime() - rt);
  rt = realtime();

  /* First pass over graph. Store all kmers in an array */
  sketch_t *sketch =
      sk_init((uint64_t)1 << (big_f ? 32 : 30), klen, mklen); // XXX: hardcoded
  uint64_t *kmers = sketch->vls; // we are going to "reuse" values array
  uint64_t totkmers = 0;         // index for insertion in kmers
  seg_t *seg = init_seg();

  char *kmer =
      malloc(sizeof(char) * (klen + 1)); // first kmer on sequence (plain)
  uint64_t kmer_d = 0;                   // kmer
  uint64_t rckmer_d = 0;                 // reverse and complemented kmer
  uint64_t ckmer_d = 0;                  // canonical kmer
  uint8_t c;                             // new character to append
  int p = 0;                             // current position on segment
  int nvertices = 0;                     // number of vertices processed

  kstring_t s = {0, 0, 0};
  int dret;
  gzFile fp = gzopen(gfa_fn, "r");
  if (fp == 0)
    return 0;
  kstream_t *ks = ks_init(fp);

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
      if (seg->l < klen)
        continue;

      strncpy(kmer, seg->seq, klen);
      kmer_d = k2d(kmer, klen);
      rckmer_d = rc(kmer_d, klen);
      ckmer_d = MIN(kmer_d, rckmer_d);
      kmers[totkmers++] = ckmer_d;

      for (p = klen; p < seg->l; ++p) {
        c = to_int[seg->seq[p]] - 1; // A is 1 but it should be 0
        kmer_d = lsappend(kmer_d, c, klen);
        rckmer_d = rsprepend(rckmer_d, reverse_char(c), klen);
        ckmer_d = MIN(kmer_d, rckmer_d);
        kmers[totkmers++] = ckmer_d;
      }
    }
  }
  gzclose(fp);
  ks_destroy(ks);
  fprintf(stderr, "[M::%s] loaded %ld kmers (from %d vertices) in %.3f sec\n",
          __func__, totkmers, nvertices, realtime() - rt);
  rt = realtime();

  /* Sort kmers */
  rt = realtime();
  ks_introsort(uint64_t, totkmers, kmers);
  fprintf(stderr, "[M::%s] sorted %ld kmers in %.3f sec\n", __func__, totkmers,
          realtime() - rt);

  /* Flag all repeated anchors */
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

  fprintf(stderr, "[M::%s] searching kmers in the FMD-index using %d threads\n",
          __func__, nth);
  rt = realtime();

  /* Flag all non-solid anchors by querying the FMD-index */
  char **kmers_t = malloc(nth * sizeof(char *));
  for (int i = 0; i < nth; ++i)
    kmers_t[i] = malloc(sizeof(char) * (klen + 1));
  uint8_t **qkmers_t = malloc(nth * sizeof(uint8_t *));
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

  /* Add solid anchors to sketch */
  for (uint64_t i = 0; i < totkmers; ++i) {
    kmer_d = kmers[i];
    if (kmer_d == 0)
      continue;
    sk_add(sketch, kmer_d);
  }
  fprintf(stderr, "[M::%s] added %ld kmers to sketch in %.3f sec\n", __func__,
          sketch->n, realtime() - rt);

  /* Reiterate over graph to assign values to solid anchors */
  fp = gzopen(gfa_fn, "r");
  ks = ks_init(fp);

  rt1 = realtime();
  nvertices = 0;
  while (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) >= 0) {
    if (s.s[0] == 'S') {
      ++nvertices;
      if (nvertices % 1000000 == 0) {
        fprintf(stderr, "[M::%s] parsed 1000000 vertices in %.3f sec\n",
                __func__, realtime() - rt1);
        rt1 = realtime();
      }

      gfa_parse_S(s.s, seg);
      if (seg->l < klen)
        continue;

      strncpy(kmer, seg->seq, klen);
      kmer_d = k2d(kmer, klen);
      rckmer_d = rc(kmer_d, klen);
      ckmer_d = MIN(kmer_d, rckmer_d);
      sk_add_v(sketch, ckmer_d, seg->idx, 0);

      for (p = klen; p < seg->l; ++p) {
        c = to_int[seg->seq[p]] - 1; // A is 1 but it should be 0
        kmer_d = lsappend(kmer_d, c, klen);
        rckmer_d = rsprepend(rckmer_d, reverse_char(c), klen);
        ckmer_d = MIN(kmer_d, rckmer_d);
        sk_add_v(sketch, ckmer_d, seg->idx, p - klen + 1);
      }
    }
  }
  gzclose(fp);
  ks_destroy(ks);
  fprintf(stderr, "[M::%s] sketched %ld kmers in %.3f sec\n", __func__,
          sketch->n, realtime() - rt);

  /* Dump sketch */
  rt = realtime();
  if (txt_f)
    sk_dump(sketch, "-");
  else
    sk_store(sketch, "-");
  fprintf(stderr, "[M::%s] dumped sketch in %.3f sec\n", __func__,
          realtime() - rt);

  /* Clean everything */
  free(kmer);
  free(s.s);
  destroy_seg(seg);
  rb3_fmi_free(&fmd);
  free(qints);
  sk_destroy(sketch);

  fprintf(stderr, "[M::%s] completed in %.3f sec\n", __func__,
          realtime() - rt0);

  return 0;
}
