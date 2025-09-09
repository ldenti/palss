#include <stdint.h>
#include <string.h>
#include <zlib.h>

#include "kseq.h"

#include "kmer.hpp"
#include "misc.hpp"
#include "sketch.hpp"
#include "usage.h"

KSEQ_INIT(gzFile, gzread)

int main_kan(int argc, char *argv[]) {
  int klen = 27; // kmer size
  // int print = 0; // print kmers
  int out_bed = 1;

  int _c;
  while ((_c = getopt(argc, argv, "k:rh")) != -1) {
    switch (_c) {
    case 'k':
      klen = std::stoi(optarg);
      break;
    case 'r':
      out_bed = 0;
      break;
    case 'h':
      fprintf(stderr, "%s", KAN_USAGE_MESSAGE);
      return 0;
    }
  }
  if (argc - optind != 2) {
    fprintf(stderr, "%s", KAN_USAGE_MESSAGE);
    return 1;
  }
  char *skt_fn = argv[optind++];
  char *fx_fn = argv[optind++];

  double rt0, rt;
  rt0 = realtime();
  rt = rt0;

  // Graph sketching
  sketch_t *sketch = sk_load(skt_fn);
  fprintf(stderr, "[M::%s] loaded %lu sketches in %.3f sec\n", __func__,
          sketch->n, realtime() - rt);
  rt = realtime();

  // ---

  gzFile fp = gzopen(fx_fn, "r");
  kseq_t *seq = kseq_init(fp);
  int l;
  char kmer[klen + 1];
  kmer[klen] = '\0';
  hit_t hit;
  uint64_t kmer_d = 0, rckmer_d = 0, ckmer_d = 0;
  uint8_t c; // new character to append
  uint p;    // current position on sequence
  int last_uncovered_p = -1;
  int tot = 0;
  int qidx = 0;
  while ((l = kseq_read(seq)) >= 0) {
    last_uncovered_p = -1;

    strncpy(kmer, seq->seq.s, klen);
    kmer_d = k2d(kmer, klen);
    rckmer_d = rc(kmer_d, klen);
    ckmer_d = std::min(kmer_d, rckmer_d);

    hit = sk_get(sketch, ckmer_d);
    if (hit.first == -1) {
      last_uncovered_p = 0;
      ++tot;
    }
    for (p = klen; p < seq->seq.l; ++p) {
      c = to_int[seq->seq.s[p]] - 1; // A is 1 but it should be 0
      kmer_d = lsappend(kmer_d, c, klen);
      rckmer_d = rsprepend(rckmer_d, reverse_char(c), klen);
      ckmer_d = std::min(kmer_d, rckmer_d);
      hit = sk_get(sketch, ckmer_d);

      if (hit.first == -1) {
        if (last_uncovered_p == -1)
          last_uncovered_p = p - klen + 1;
        ++tot;
      } else {
        if (last_uncovered_p != -1) {
          if (out_bed) {
            // half-open interval, [)
            printf("%s\t%d\t%d\n", seq->name.s, last_uncovered_p, p - klen + 1);
          }
        }
        last_uncovered_p = -1;
      }

      if (out_bed && p % 10000000 == 0)
        fprintf(stderr, "Read %d bases from %s\r", p, seq->name.s);
    }
    if (out_bed)
      fprintf(stderr, "Read %d bases from %s\n", seq->seq.l, seq->name.s);
    else {
      printf("%s\t%d\t%d\n", seq->name.s, tot, l - klen + 1);
      tot = 0;
      ++qidx;
      if (qidx % 10000 == 0) {
        fprintf(stderr, "[M::%s] parsed 10000 reads %.3f sec\n", __func__,
                realtime() - rt);
        rt = realtime();
      }
    }
  }

  kseq_destroy(seq);
  gzclose(fp);
  sk_destroy(sketch);
  // ---

  fprintf(stderr, "[M::%s] completed in %.3f sec\n", __func__,
          realtime() - rt0);

  return 0;
}
