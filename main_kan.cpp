#include <stdint.h>
#include <string.h>
#include <zlib.h>

#include "kseq.h"

#include "kmer.hpp"
#include "misc.hpp"
#include "sketch.hpp"
#include "usage.hpp"

KSEQ_INIT(gzFile, gzread)

int main_kan(int argc, char *argv[]) {
  // int print = 0; // print kmers
  bool in_fa = true;
  bool reference_only = false;

  int _c;
  while ((_c = getopt(argc, argv, "qh")) != -1) {
    switch (_c) {
    // case 'r':
    //   reference_only = true;
    //   break;
    case 'q':
      in_fa = false;
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

  double rt;
  rt = realtime();

  // Sketch
  rt = realtime();
  sketch_t *sketch = sk_load(skt_fn);
  size_t klen = sketch->k;
  fprintf(stderr, "[M::%s] Restored sketch in %.3f sec\n", __func__,
          realtime() - rt);

  // ---

  rt = realtime();

  gzFile fp = gzopen(fx_fn, "r");
  kseq_t *seq = kseq_init(fp);
  int l;
  char kmer[klen + 1];
  kmer[klen] = '\0';

  hit_t hit; // hit from sketch
  uint64_t kmer_d = 0, rckmer_d = 0, ckmer_d = 0;
  uint8_t c; // new character to append
  uint p;    // current position on sequence
  int last_uncovered_p = -1;
  int missed = 0;
  int qidx = 0;
  while ((l = kseq_read(seq)) >= 0) {
    last_uncovered_p = -1;

    strncpy(kmer, seq->seq.s, klen);
    kmer_d = k2d(kmer, klen);
    rckmer_d = rc(kmer_d, klen);
    ckmer_d = std::min(kmer_d, rckmer_d);

    hit = sk_get(sketch, ckmer_d, reference_only);
    if (hit.value == -1UL) {
      last_uncovered_p = 0;
      ++missed;
    }
    for (p = klen; p < seq->seq.l; ++p) {
      c = to_int[(uint8_t)seq->seq.s[p]] - 1; // A is 1 but it should be 0
      kmer_d = lsappend(kmer_d, c, klen);
      rckmer_d = rsprepend(rckmer_d, reverse_char(c), klen);
      ckmer_d = std::min(kmer_d, rckmer_d);
      hit = sk_get(sketch, ckmer_d, reference_only);
      if (hit.value == -1UL) {
        if (last_uncovered_p == -1)
          last_uncovered_p = p - klen + 1;
        ++missed;
      } else {
        if (last_uncovered_p != -1) {
          if (in_fa) {
            // half-open interval, [)
            printf("%s\t%d\t%d\n", seq->name.s, last_uncovered_p, p - klen + 1);
          }
        }
        last_uncovered_p = -1;
      }

      if (in_fa && p % 10000000 == 0) {
        fprintf(stderr, "Read %d bases from %s in %.3f sec\r", p, seq->name.s,
                realtime() - rt);
        rt = realtime();
      }
    }
    if (in_fa) {
      fprintf(stderr, "Read %d bases from %s in %.3f sec\n", p, seq->name.s,
              realtime() - rt);
      rt = realtime();
    } else {
      printf("%s\t%d\t%d\n", seq->name.s, missed, l - klen + 1);
      missed = 0;
      ++qidx;
      if (qidx % 10000 == 0) {
        fprintf(stderr, "[M::%s] parsed 10000 reads in %.3f sec\n", __func__,
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
