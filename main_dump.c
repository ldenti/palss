#include <stdint.h>
#include <string.h>

#include "misc.h"
#include "sketch.h"

/* TODO: print kmers as strings */

int main_dump(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "./palss dump <.skt>\n");
    return 1;
  }
  char *skt_fn = argv[1];

  double rt0, rt;
  rt0 = realtime();
  rt = rt0;

  rt = realtime();
  sketch_t *sketch = sk_load(skt_fn);
  fprintf(stderr, "[M::%s] loaded %lu sketches in %.3f sec\n", __func__,
          sketch->n, realtime() - rt);

  sk_dump(sketch, "-");

  fprintf(stderr, "[M::%s] completed in %.3f sec\n", __func__,
          realtime() - rt0);

  return 0;
}
