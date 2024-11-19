#include "graph.h"

seg_t *init_seg() {
  seg_t *seg = malloc(1 * sizeof(seg_t));
  seg->l = 0;
  // seg->idx = (char *)malloc(1024);
  seg->seq = malloc(4096 * sizeof(char));
  seg->c = 4096;
  return seg;
}
void destroy_seg(seg_t *seg) {
  free(seg->seq);
  free(seg);
}

void gfa_parse_S(char *s, seg_t *ret) {
  int i, is_ok = 0;
  char *p, *q, *seg = 0, *seq = 0, *rest = 0;
  uint32_t sid, len = 0;
  for (i = 0, p = q = s + 2;; ++p) {
    if (*p == 0 || *p == '\t') {
      int c = *p;
      *p = 0;
      if (i == 0) {
        ret->idx = atoi(q);
        // strcpy(ret->idx, q);
      } else if (i == 1) {
        // TODO: reallocate if vertex is longer than 4096
        // right now we assume to have a vg chopped graph
        strcpy(ret->seq, q);
        ret->l = p - q;
        is_ok = 1, rest = c ? p + 1 : 0;
        break;
      }
      ++i, q = p + 1;
      if (c == 0)
        break;
    }
  }
  // if (!is_ok) { // something is missing
}

