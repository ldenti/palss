#include <bitset>
#include <iostream>
#include <stdint.h>
#include <string.h>

#include "sketch.hpp"

int main_stats(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "./palss stats <.skt>\n");
    return 1;
  }
  char *skt_fn = argv[1];

  sketch_t *sketch = sk_load(skt_fn);
  printf("Total number of sketches: %ld\n", sketch->n);

  // for (int i = 0; i < sketch->n; ++i)
  //   printf("%ld\n", sketch->sxs[i]);

  sk_destroy(sketch);

  return 0;
}
