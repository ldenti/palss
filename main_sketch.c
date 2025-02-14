#include <assert.h>
#include <zlib.h>

#include "fm-index.h"
#include "ketopt.h"
#include "kseq.h"

#include "graph.h"
#include "kmer.h"
#include "misc.h"
#include "sketch.h"
#include "usage.h"

/* Backward search */
int search_full(const rb3_fmi_t *fmd, uint8_t *kmer, int k, int min_size) {
  rb3_sai_t ik;
  int begin = k - 1;
  rb3_fmd_set_intv(fmd, kmer[begin], &ik);
  /* printf("%ld:%ld:%ld\n", ik.x[0], ik.x[1], ik.size); */
  while (ik.size != 0 && begin > 0) {
    --begin;
    rb3_sai_t ok[RB3_ASIZE]; // output SA intervals (one for each symbol
                             // between 0 and 5)
    rb3_fmd_extend(fmd, &ik, ok, 1);
    ik = ok[kmer[begin]];
    /* printf("%d - %ld:%ld:%ld\n", kmer[begin], ik.x[0], ik.x[1], ik.size); */
    if (ik.size <= min_size)
      return ik.size;
  }
  return ik.size;
}

/* Backward search */
int search(const rb3_fmi_t *fmd, rb3_sai_t ik, uint8_t *kmer, int k,
           int min_size) {
  int begin = k;
  /* printf("%ld:%ld:%ld\n", ik.x[0], ik.x[1], ik.size); */
  while (ik.size != 0 && begin > 0) {
    --begin;
    rb3_sai_t ok[RB3_ASIZE]; // output SA intervals (one for each symbol
                             // between 0 and 5)
    rb3_fmd_extend(fmd, &ik, ok, 1);
    ik = ok[kmer[begin]];
    /* printf("%d - %ld:%ld:%ld\n", kmer[begin], ik.x[0], ik.x[1], ik.size); */
    if (ik.size <= min_size)
      return ik.size;
  }
  return ik.size;
}

/* Backward search */
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

  int klen = 27; // kmer size
  int mklen = 7; // memoization kmer size
  int nh = 1;    // expected number of haplotypes
  int txt_f = 0; // 1: dump sketch in txt
  int threads = 4;
  static ko_longopt_t longopts[] = {{NULL, 0, 0}};
  ketopt_t opt = KETOPT_INIT;
  int _c;
  while ((_c = ketopt(&opt, argc, argv, 1, "k:q:g:@:th", longopts)) >= 0) {
    if (_c == 'k')
      klen = atoi(opt.arg);
    else if (_c == 'q')
      mklen = atoi(opt.arg);
    else if (_c == 'g')
      nh = atoi(opt.arg);
    else if (_c == 't')
      txt_f = 1;
    else if (_c == '@')
      threads = atoi(opt.arg);
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

  // memoization
  rb3_sai_t *qints = memoize(&fmd, mklen);
  fprintf(stderr, "[M::%s] memoization (k=%d) in %.3f sec\n", __func__, mklen,
          realtime() - rt);
  rt = realtime();
  /* for (int i = 0; i < (1 << (mklen * 2)); ++i) */
  /*   printf("%d\t%ld:%ld %ld\n", i, qints[i].x[0], qints[i].x[1],
   * qints[i].size); */
  // ---

  // Sketching
  kstring_t s = {0, 0, 0};
  int dret;
  gzFile fp = gzopen(gfa_fn, "r");
  if (fp == 0)
    return 0;
  kstream_t *ks = ks_init(fp);

  sketch_t *sketch = sk_init();
  seg_t *seg = init_seg();

  char *kmer =
      malloc(sizeof(char) * (klen + 1)); // first kmer on sequence (plain)
  uint8_t *qkmer;
  uint64_t kmer_d = 0;   // kmer
  uint64_t rckmer_d = 0; // reverse and complemented kmer
  uint64_t ckmer_d = 0;  // canonical kmer
  rb3_sai_t qint;        // memoization interval
  uint8_t c;             // new character to append
  int p = 0;             // current position on segment
  int hits = 0;          // hits in the FMD index

  rt1 = realtime();
  int nvertices = 0;
  while (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) >= 0) {
    if (s.s[0] == 'S') {
      ++nvertices;
      if (nvertices % 25000 == 0) {
        fprintf(stderr, "[M::%s] sketched 25000 vertices in %.3f sec\n",
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

      qkmer = (uint8_t *)kmer;
      rb3_char2nt6(klen, qkmer);
      qint = qints[kmer_d & ((1 << (mklen * 2)) - 1)];
      hits = 1;
      if (qint.size > nh) {
        hits = search(&fmd, qint, qkmer, klen - mklen, nh);
        /* int hitsf = search_full(&fmd, qkmer, klen, nh); */
        /* assert(hits == hitsf); */
        assert(hits > 0);
      }

      if (hits <= nh)
        // TODO: keep track of number of skipped anchors
        sk_add(sketch, ckmer_d, seg->idx, 0, 1);

      for (p = klen; p < seg->l; ++p) {
        c = to_int[seg->seq[p]] - 1; // A is 1 but it should be 0
        kmer_d = lsappend(kmer_d, c, klen);
        rckmer_d = rsprepend(rckmer_d, reverse_char(c), klen);
        ckmer_d = MIN(kmer_d, rckmer_d);

        d23(kmer_d, klen, kmer);
        qkmer = (uint8_t *)kmer;
        rb3_char2nt6(klen, qkmer);
        qint = qints[kmer_d & ((1 << (mklen * 2)) - 1)];
        hits = 1;
        if (qint.size > nh) {
          hits = search(&fmd, qint, qkmer, klen - mklen, nh);
          /* hitsf = search_full(&fmd, qkmer, klen, nh); */
          /* assert(hits == hitsf); */
          assert(hits > 0);
        }
        // TODO: keep track of number of skipped anchors
        if (hits <= nh)
          sk_add(sketch, ckmer_d, seg->idx, p - klen + 1, 1);
      }
    }
  }

  fprintf(stderr, "[M::%s] sketched %d anchors in %.3f sec\n", __func__,
          kh_size(sketch), realtime() - rt);
  rt = realtime();

  // ---

  sk_clean(sketch);
  fprintf(stderr, "[M::%s] cleaned sketch in %.3f sec\n", __func__,
          realtime() - rt);
  rt = realtime();
  // Write sketch to stdout
  if (txt_f)
    sk_dump(sketch, "-");
  else
    sk_store(sketch, "-");
  fprintf(stderr, "[M::%s] dumped sketch (%d anchors) in %.3f sec\n", __func__,
          kh_size(sketch), realtime() - rt);

  // ---

  // Clean everything
  free(kmer);
  free(s.s);
  destroy_seg(seg);
  ks_destroy(ks);
  gzclose(fp);
  rb3_fmi_free(&fmd);
  free(qints);
  sk_destroy(sketch);

  fprintf(stderr, "[M::%s] completed in %.3f sec\n", __func__,
          realtime() - rt0);

  return 0;
}
