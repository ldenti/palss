#include <stdint.h>
#include <string>
#include <zlib.h>

extern "C" {
#include "kmer.h"
#include "kseq.h"
#include "sketch.h"
}
#include "misc.hpp"
#include "usage.hpp"

KSEQ_INIT(gzFile, gzread)

int main_kan(int argc, char *argv[]) {
  int klen = 27;      // kmer size
  bool print = false; // print kmers
  bool out_bed = true;

  int _c;
  while ((_c = getopt(argc, argv, "k:rph")) != -1) {
    switch (_c) {
    case 'k':
      klen = std::stoi(optarg);
      break;
    case 'p':
      print = true;
      break;
    case 'r':
      out_bed = false;
      break;
    case 'h':
      fprintf(stderr, "%s", KAN_USAGE_MESSAGE);
      return 0;
    }
  }
  if (argc - optind != 2) {
    fprintf(stderr, KAN_USAGE_MESSAGE);
    return 1;
  }
  char *skt_fn = argv[optind++];
  char *fx_fn = argv[optind++];

  double rt = realtime();

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
  int p;     // current position on sequence
  int last_uncovered_p = -1;
  int tot = 0;
  int qidx = 0;
  while ((l = kseq_read(seq)) >= 0) {
    last_uncovered_p = -1;

    strncpy(kmer, seq->seq.s, klen);
    kmer_d = k2d(kmer, klen);
    rckmer_d = rc(kmer_d, klen);
    ckmer_d = MIN(kmer_d, rckmer_d);

    hit = sk_get(sketch, ckmer_d);
    if (hit.first == -1) {
      last_uncovered_p = 0;
      ++tot;
    }
    for (p = klen; p < seq->seq.l; ++p) {
      c = to_int[seq->seq.s[p]] - 1; // A is 1 but it should be 0
      kmer_d = lsappend(kmer_d, c, klen);
      rckmer_d = rsprepend(rckmer_d, reverse_char(c), klen);
      ckmer_d = MIN(kmer_d, rckmer_d);
      hit = sk_get(sketch, ckmer_d);

      if (hit.first == -1) {
        if (last_uncovered_p == -1)
          last_uncovered_p = p - klen + 1;
        ++tot;
      } else {
        if (last_uncovered_p != -1) {
          if (out_bed) {
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

  return 0;
}
