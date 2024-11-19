#include <cstdint>
#include <cstring>
#include <zlib.h>

#include "ketopt.h"
#include "kseq.h"

#include "sketch.hpp"
#include "utils.h"

KSEQ_INIT(gzFile, gzread)

using namespace std;

int main_kan(int argc, char *argv[]) {
  int klen = 27; // kmer size
  int print = 0; // print kmers
  int verbose = 0;
  static ko_longopt_t longopts[] = {{NULL, 0, 0}};
  ketopt_t opt = KETOPT_INIT;
  int _c;
  while ((_c = ketopt(&opt, argc, argv, 1, "k:pv", longopts)) >= 0) {
    if (_c == 'k')
      klen = atoi(opt.arg);
    else if (_c == 'p')
      print = 1;
    else if (_c == 'v')
      verbose = 1;
  }

  if (argc - opt.ind != 2) {
    fprintf(stderr, "Argh");
    return 1;
  }
  char *skt_fn = argv[opt.ind++];
  char *fa_fn = argv[opt.ind++];

  double rt0, rt;
  rt0 = realtime();
  rt = rt0;

  // Graph sketching and path extraction
  sketch_t sketch;
  sk_load(sketch, skt_fn);
  fprintf(stderr, "[M::%s] loaded %lu sketches in %.3f sec\n", __func__,
          sketch.size(), realtime() - rt);
  rt = realtime();

  // ---

  gzFile fp = gzopen(fa_fn, "r");
  kseq_t *seq = kseq_init(fp);
  int l;
  char kmer[klen + 1];
  kmer[klen] = '\0';
  hit_t hit;
  uint64_t kmer_d = 0, rckmer_d = 0, ckmer_d = 0;
  uint8_t c; // new character to append
  int p;     // current position on chromosome
  while ((l = kseq_read(seq)) >= 0) {
    p = 0;
    strncpy(kmer, seq->seq.s, klen);
    kmer_d = k2d(kmer, klen);
    rckmer_d = rc(kmer_d, klen);
    ckmer_d = std::min(kmer_d, rckmer_d);
    hit = sk_get(sketch, ckmer_d);
    if (hit.first != -1)
      printf("%s:%d-%d %lu:%d %s\n", seq->name.s, p + 1, p + klen, hit.first,
             hit.second, print ? d2s(kmer_d, klen).c_str() : "");
    for (p = klen; p < seq->seq.l; ++p) {
      c = to_int[seq->seq.s[p]] - 1; // A is 1 but it should be 0
      kmer_d = lsappend(kmer_d, c, klen);
      rckmer_d = rsprepend(rckmer_d, reverse_char(c), klen);
      ckmer_d = std::min(kmer_d, rckmer_d);
      hit = sk_get(sketch, ckmer_d);
      if (hit.first != -1)
        printf("%s:%d-%d %lu:%d %s\n", seq->name.s, p - klen + 1 + 1, p + 1,
               hit.first, hit.second, print ? d2s(kmer_d, klen).c_str() : "");

      if (verbose && p % 500000 == 0)
        fprintf(stderr, "Read %d bases from %s\r", p, seq->name.s);
    }
    if (verbose)
      fprintf(stderr, "Read %d bases from %s\n", seq->seq.l, seq->name.s);
  }
  kseq_destroy(seq);
  gzclose(fp);

  // ---

  fprintf(stderr, "[M::%s] completed in %.3f sec\n", __func__,
          realtime() - rt0);

  return 0;
}
