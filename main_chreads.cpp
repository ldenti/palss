#include <zlib.h>

#include "ketopt.h"
#include "kseq.h"

#include "sketch.hpp"
#include "utils.h"

KSEQ_INIT(gzFile, gzread)

using namespace std;

int main_chreads(int argc, char *argv[]) {
  double rt0 = realtime();
  double rt = rt0, rt1;

  int klen = 27; // kmer size
  static ko_longopt_t longopts[] = {};
  ketopt_t opt = KETOPT_INIT;
  int _c;
  while ((_c = ketopt(&opt, argc, argv, 1, "k:", longopts)) >= 0) {
    if (_c == 'k')
      klen = atoi(opt.arg);
  }
  if (argc - opt.ind != 2) {
    fprintf(stderr, "Argh");
    return 1;
  }
  char *skt_fn = argv[opt.ind++];
  char *fq_fn = argv[opt.ind++];

  // Graph sketching and path extraction
  sketch_t sketch;
  sk_load(sketch, skt_fn);
  fprintf(stderr, "[M::%s] loaded %ld sketches in %.3f sec\n", __func__,
          sketch.size(), realtime() - rt);
  rt = realtime();
  rt1 = rt;
  // ---

  // Specific strings computation and anchoring
  gzFile fp = gzopen(fq_fn, "r");
  kseq_t *seq = kseq_init(fp);
  int l;
  uint qidx = 0;

  char *kmer = (char *)malloc(sizeof(char) *
                              (klen + 1)); // first kmer on sequence (plain)
  uint64_t kmer_d = 0;                     // kmer
  uint64_t rckmer_d = 0;                   // reverse and complemented kmer
  uint64_t ckmer_d = 0;                    // canonical kmer
  uint8_t c;                               // new character to append
  int p = 0;                               // current position on read
  int n = 0;                               // number of anchors on read
  hit_t hit;
  while ((l = kseq_read(seq)) >= 0) {
    strncpy(kmer, seq->seq.s, klen);
    kmer_d = k2d(kmer, klen);
    rckmer_d = rc(kmer_d, klen);
    ckmer_d = std::min(kmer_d, rckmer_d);

    n = sk_get(sketch, ckmer_d).first != -1;
    for (p = klen; p < l; ++p) {
      c = to_int[seq->seq.s[p]] - 1; // A is 1 but it should be 0
      kmer_d = lsappend(kmer_d, c, klen);
      rckmer_d = rsprepend(rckmer_d, reverse_char(c), klen);
      ckmer_d = std::min(kmer_d, rckmer_d);
      n += sk_get(sketch, ckmer_d).first != -1;
    }
    printf("%s\t%d\t%d\n", seq->name.s, n, l - klen + 1);
    ++qidx;
    if (qidx % 10000 == 0) {
      fprintf(stderr, "[M::%s] parsed %d reads %.3f sec\n", __func__, qidx,
              realtime() - rt1);
      rt1 = realtime();
    }
  }
  free(kmer);
  kseq_destroy(seq);
  gzclose(fp);
  fprintf(stderr, "[M::%s] done in %.3f sec\n", __func__, realtime() - rt0);

  return 0;
}
