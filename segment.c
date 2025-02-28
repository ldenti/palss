#include "segment.h"

seg_t *init_seg() {
  seg_t *seg = malloc(1 * sizeof(seg_t));
  seg->l = 0;
  // seg->idx = (char *)malloc(1024);
  seg->seq = NULL; /* malloc(4096 * sizeof(char)); */
  seg->c = 0;
  kv_init(seg->paths);
  return seg;
}

void destroy_seg(seg_t *seg) {
  if (seg->seq != NULL)
    free(seg->seq);
  kv_destroy(seg->paths);
  free(seg);
}
