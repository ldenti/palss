#include <cstdint>
#include <cstring>
#include <map>
#include <zlib.h>

#include "gsketch.hpp"
#include "ketopt.h"
#include "kseq.h"
#include "utils.h"

// KSEQ_INIT(gzFile, gzread) // we already init kstream in gsketch
__KSEQ_TYPE(gzFile)
__KSEQ_BASIC(static, gzFile)
__KSEQ_READ(static)

using namespace std;

int main_map(int argc, char *argv[]) {
  // int klen = 27; // kmer size
  // char pidx[128] = "";
  // static ko_longopt_t longopts[] = {{NULL, 0, 0}};
  // ketopt_t opt = KETOPT_INIT;
  // int _c;
  // while ((_c = ketopt(&opt, argc, argv, 1, "k:p:", longopts)) >= 0) {
  //   if (_c == 'k')
  //     klen = atoi(opt.arg);
  //   else if (_c == 'p')
  //     strncpy(pidx, opt.arg, strlen(opt.arg));
  // }

  // if (argc - opt.ind != 3) {
  //   fprintf(stderr, "Argh");
  //   return 1;
  // }
  // char *gfa_fn = argv[opt.ind++];
  // char *skt_fn = argv[opt.ind++];
  // char *sfs_fn = argv[opt.ind++];

  // double rt0, rt;
  // rt0 = realtime();
  // rt = rt0;

  // // Graph sketching and path extraction
  // GSK gsk(gfa_fn, klen);
  // gsk.load_vertices();
  // fprintf(stderr, "[M::%s] loaded %d vertices in %.3f sec\n", __func__,
  //         gsk.nvertices, realtime() - rt);
  // rt = realtime();

  // gsk.load_sketch(skt_fn);
  // fprintf(stderr, "[M::%s] loaded %lu sketches in %.3f sec\n", __func__,
  //         gsk.sketch.size(), realtime() - rt);
  // rt = realtime();

  // gsk.load_paths();
  // fprintf(stderr, "[M::%s] loaded %ld paths in %.3f sec\n", __func__,
  //         gsk.paths.size(), realtime() - rt);
  // rt = realtime();

  // // ---

  // map<int, int> positions = gsk.get_positions(pidx);
  // int totl = 0;
  // for (const auto &p : positions)
  //   totl += gsk.get_vl(p.first);
  // fprintf(stderr, "[M::%s] extracted path '%s' in %.3f sec\n", __func__,
  // pidx,
  //         realtime() - rt);
  // rt = realtime();

  // // ---
  // gzFile fp = gzopen(sfs_fn, "r");
  // kseq_t *seq = kseq_init(fp);
  // int l;
  // char k1[klen + 1], k2[klen + 1];
  // k1[klen] = '\0';
  // k2[klen] = '\0';

  // uint64_t k1_d, k2_d, kmer_d, rckmer_d;
  // pair<int64_t, int16_t> p1, p2;
  // uint64_t v1, v2;
  // int pos1, pos2, d;

  // printf("@HD\tVN:1.6\tSO:coordinate\n");
  // printf("@SQ\tSN:%s\tLN:%d\n", "chr20", totl);
  // while ((l = kseq_read(seq)) >= 0) {
  //   strncpy(k1, seq->seq.s, klen);
  //   kmer_d = k2d(k1, klen);
  //   rckmer_d = rc(kmer_d, klen);
  //   k1_d = std::min(kmer_d, rckmer_d);

  //   strncpy(k2, seq->seq.s + (l - klen), klen);
  //   kmer_d = k2d(k2, klen);
  //   rckmer_d = rc(kmer_d, klen);
  //   k2_d = std::min(kmer_d, rckmer_d);

  //   p1 = gsk.get(k1_d);
  //   v1 = p1.first;
  //   p2 = gsk.get(k2_d);
  //   v2 = p2.first;

  //   if (v1 == -1 || v2 == -1)
  //     // unanchored specific string
  //     continue;

  //   pos1 = positions.at(v1) + p1.second;
  //   pos2 = positions.at(v2) + p2.second;

  //   d = pos2 - (pos1 + klen);

  //   printf("%s\t0\t%s\t%d\t60\t%dM%dN%dM\t*\t0\t0\t%s%s\t*\n", seq->name.s,
  //          "chr20", pos1 + 1, klen, d, klen, k1, k2);
  // }
  // kseq_destroy(seq);
  // gzclose(fp);

  // // ---

  // gsk.destroy_paths();
  // fprintf(stderr, "[M::%s] completed in %.3f sec\n", __func__,
  //         realtime() - rt0);

  return 0;
}
