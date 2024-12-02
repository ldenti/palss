#include <assert.h>
#include <cstdint>
#include <map>
#include <utility>
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

using namespace std;

int main_call(int argc, char *argv[]) {
  double rt0 = realtime();
  double rt = rt0;

  int klen = 27;  // kmer size
  int nh = 1;     // expected number of haplotypes
  int vpb = 5000; // number of vertices to load per batch
  int threads = 4;
  static ko_longopt_t longopts[] = {{NULL, 0, 0}};
  ketopt_t opt = KETOPT_INIT;
  int _c;
  while ((_c = ketopt(&opt, argc, argv, 1, "k:g:v:@:", longopts)) >= 0) {
    if (_c == 'k')
      klen = atoi(opt.arg);
    else if (_c == 'g')
      nh = atoi(opt.arg);
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

  // // Sketching
  // kstring_t s = {0, 0, 0};
  // int dret;
  // gzFile fp = gzopen(gfa_fn, "r");
  // if (fp == 0)
  //   return 0;
  // kstream_t *ks = ks_init(fp);

  // sketch_t sketch;
  // seg_t **segs = (seg_t **)malloc(vpb * sizeof(seg_t *));
  // for (int i = 0; i < vpb; ++i)
  //   segs[i] = init_seg();
  // int si = 0; // current segment
  // int nvertices = 0;
  // while (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) >= 0) {
  //   if (s.s[0] == 'S') {
  //     ++nvertices;
  //     gfa_parse_S(s.s, segs[si]);
  //     if (segs[si]->l < klen)
  //       continue;
  //     ++si;

  //     if (si == vpb) {
  //       run_sketching(segs, si, klen, &fmd, nh, threads, sketch);
  //       si = 0;
  //     }
  //   }
  // }
  // if (si < vpb) {
  //   run_sketching(segs, si, klen, &fmd, nh, threads, sketch);
  // }

  // // Clean everything
  // for (int i = 0; i < vpb; ++i)
  //   destroy_seg(segs[i]);
  // free(segs);
  // free(s.s);
  // ks_destroy(ks);
  // gzclose(fp);
  // rb3_fmi_free(&fmd);

  // fprintf(stderr, "[M::%s] sketched graph in %.3f sec\n", __func__,
  //         realtime() - rt);
  // rt = realtime();
  // // ---

  // // Write sketch to stdout
  // sk_store(as_const(sketch), "-");
  // // ---

  fprintf(stderr, "[M::%s] completed in %.3f sec\n", __func__,
          realtime() - rt0);

  return 0;
}
