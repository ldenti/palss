#include <cstdint>
#include <iostream>

#include "fm-index.h"
#include "ketopt.h"

#include "gsketch.hpp"
#include "utils.h"

/* Backward search */
int search(const rb3_fmi_t *index, char *kmer, int k) {
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

int main_sketch(int argc, char *argv[]) {
  int k = 27; // kmer size
  int G = 1;  // expected number of genomes
  int fa = 0; // fasta output
  static ko_longopt_t longopts[] = {{NULL, 0, 0}};
  ketopt_t opt = KETOPT_INIT;
  int _c;
  while ((_c = ketopt(&opt, argc, argv, 1, "k:g:f", longopts)) >= 0) {
    if (_c == 'k')
      k = atoi(opt.arg);
    else if (_c == 'g')
      G = atoi(opt.arg);
    else if (_c == 'f')
      fa = 1;
  }

  if (argc - opt.ind != 2) {
    fprintf(stderr, "Argh");
    return 1;
  }
  char *gfa_fn = argv[opt.ind++];
  char *fmd_fn = argv[opt.ind++];

  double rt = realtime();
  GSK gsk(gfa_fn, k);
  gsk.build_sketch();
  fprintf(stderr, "[M::%s] sketched graph with %d vertices in %.3f sec\n",
          __func__, gsk.nvertices, realtime() - rt);
  rt = realtime();

  // FMD-index loading
  rb3_fmi_t f;
  rb3_fmi_restore(&f, fmd_fn, 0);
  if (f.e == 0 && f.r == 0) {
    fprintf(stderr, "Error restoring index");
    return 1;
  }
  fprintf(stderr, "[M::%s] restored FMD index in %.3f sec\n", __func__,
          realtime() - rt);
  rt = realtime();
  // ---

  // Retain only really unique kmers (wrt genomes)
  char *kmer = (char *)malloc((k + 1) * sizeof(char));
  kmer[k] = '\0';
  int hits;
  int not_unique = 0;
  for (auto &it : gsk.sketch) {
    d23(it.first, k, kmer);
    hits = search(&f, kmer, k);
    assert(hits > 0);
    if (hits != G) {
      it.second = it.second & ~1;
      ++not_unique;
    }
  }
  fprintf(stderr,
          "[M::%s] backward searched %ld kmers (%d not unique) in %.3f sec\n",
          __func__, gsk.sketch.size(), not_unique, realtime() - rt);
  rt = realtime();

  gsk.store_sketch(stdout, fa);
  fprintf(stderr, "[M::%s] dumped sketch (%ld unique kmers) in %.3f sec\n",
          __func__, gsk.sketch.size() - not_unique, realtime() - rt);

  rb3_fmi_free(&f);

  return 0;
}
