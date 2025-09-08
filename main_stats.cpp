#include <stdint.h>
#include <string.h>

#include "misc.hpp"
#include "sketch.hpp"

int main_stats(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "./palss stats <.skt>\n");
    return 1;
  }
  char *skt_fn = argv[1];

  double rt0, rt;
  rt0 = realtime();
  rt = rt0;

  rt = realtime();
  sketch_t *sketch = sk_load(skt_fn);

  printf("Total number of sketches: %ld\n", sketch->n);

  int ref = 0;
  for (int v = 0; v < sketch->n; ++v) {
    ref += sk_decode_isref(sketch->vls[v]);
    if (!sk_decode_isref(sketch->vls[v]))
      printf("%ld\n", sk_decode_v(sketch->vls[v]));
  }
  printf("Sketches on reference path(s): %d\n", ref);

  // sk_dump(sketch, "-");

  printf("Loading time: %.3f sec\n", realtime() - rt);

  return 0;
}
