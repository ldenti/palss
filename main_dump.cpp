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
  int solid = 0;
  for (int i = 0; i < sketch->n; ++i) {
    d2s(sketch->sxs[i], klen, kmer);
    uint32_t v1 = (uint32_t)(sketch->vls[i] >> 33);
    uint32_t v2 = (uint32_t)sketch->vls[i] >> 1;
    printf("%s\t%d>%d\t%s>%s\t%ld\t%ld\n", kmer, v1, v1,
           graph.get_gfa_name(v1 >> 1).c_str(),
           graph.get_gfa_name(v2 >> 1).c_str(), sketch->vls[i] & 1,
           sketch->sxs[i]);
    solid += sketch->vls[i] != -1U;
  }
  fprintf(stderr, "Total number of solid anchors: %d\n", solid);

  return 0;
}
