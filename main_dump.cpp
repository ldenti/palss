#include <stdint.h>
#include <string.h>

#include "graph.hpp"
#include "kmer.hpp"
#include "sketch.hpp"

/* TODO: print kmers as strings */

int main_dump(int argc, char *argv[]) {
  if (argc != 3) {
    fprintf(stderr, "./palss dump <.skt> <.gbz>\n");
    return 1;
  }
  char *skt_fn = argv[1];
  char *gbz_fn = argv[2];

  Graph graph(gbz_fn, "");
  graph.load();

  sketch_t *sketch = sk_load(skt_fn);
  int klen = sketch->k;
  fprintf(stderr, "Total number of sketches: %ld\n", sketch->n);

  char kmer[klen + 1];
  printf("kmer\ti1>i2\tv1>v2\trc\tref\n");
  size_t ref = 0;
  for (int i = 0; i < sketch->n; ++i) {
    d2s(sketch->sxs[i], klen, kmer);

    uint32_t v1 = (uint32_t)(sketch->vls[i] >> 32);
    uint32_t v2 = (uint32_t)sketch->vls[i];

    uint32_t pos1 = (sketch->info[i] >> 17) & 0x7FFF;
    uint32_t pos2 = (sketch->info[i] >> 2) & 0x7FFF;
    bool has_both = (sketch->info[i] >> 1) & 1;
    bool is_reference = sketch->info[i] & 1;
    printf("%s\t%d:%d>%d:%d\t%s>%s\t%d\t%d\n", kmer, v1, pos1, v2, pos2,
           graph.get_gfa_name(v1 >> 1).c_str(),
           graph.get_gfa_name(v2 >> 1).c_str(), has_both, is_reference);
    ref += is_reference;
  }
  fprintf(stderr, "%ld/%ld (%.3f) anchors on reference path\n", ref, sketch->n,
          (float)ref / sketch->n);

  return 0;
}
