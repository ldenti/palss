#include <cstdint>
#include <cstring>

#include "ketopt.h"
#include "kseq.h"

#include "sketch.hpp"
#include "utils.h"

using namespace std;

int main_dump(int argc, char *argv[]) {
  int klen = 27; // kmer size
  static ko_longopt_t longopts[] = {{NULL, 0, 0}};
  ketopt_t opt = KETOPT_INIT;
  int _c;
  while ((_c = ketopt(&opt, argc, argv, 1, "k:", longopts)) >= 0) {
    if (_c == 'k')
      klen = atoi(opt.arg);
  }

  if (argc - opt.ind != 1) {
    fprintf(stderr, "Argh");
    return 1;
  }
  char *skt_fn = argv[opt.ind++];

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

  for (const auto &pair : sketch)
    printf("%s\t%d\t%ld:%d\n", d2s(pair.first, klen).c_str(),
           sk_decode_unique(pair.second), sk_decode_v(pair.second),
           sk_decode_off(pair.second));

  // ---

  fprintf(stderr, "[M::%s] completed in %.3f sec\n", __func__,
          realtime() - rt0);

  return 0;
}
