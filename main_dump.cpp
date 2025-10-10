#include <stdint.h>
#include <string.h>

#include "kmer.hpp"
#include "sketch.hpp"

/* TODO: print kmers as strings */

int main_dump(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "./palss dump <.skt>\n");
    return 1;
  }
  char *skt_fn = argv[1];

  sketch_t *sketch = sk_load(skt_fn);
  int klen = sketch->k;
  fprintf(stderr, "Total number of sketches: %ld\n", sketch->n);

  char kmer[klen + 1];
  int solid = 0;
  for (int i = 0; i < sketch->n; ++i) {
    d2s(sketch->sxs[i], klen, kmer);
    printf("%s\t%d\t%ld\n", kmer, sketch->vls[i], sketch->sxs[i]);
    solid += sketch->vls[i] != -1U;
  }
  fprintf(stderr, "Total number of solid anchors: %d\n", solid);

  return 0;
}
