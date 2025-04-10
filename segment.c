#include "segment.h"

seg_t *init_seg() {
  seg_t *seg = malloc(1 * sizeof(seg_t));
  seg->l = 0;
  // seg->idx = (char *)malloc(1024);
  seg->seq = NULL; /* malloc(4096 * sizeof(char)); */
  seg->c = 0;

  /* seg->paths_p = 0; */
  /* seg->paths_b = 0; */
  /* seg->paths_c = 512; */
  /* seg->paths = calloc(seg->paths_c, 1); */
  /* memset(seg->cnts, 0, 48); */

  /* kv_init(seg->pord); */
  /* seg->cpord = NULL; */
  /* kv_init(seg->starts); */

  /* seg->lasth = -1; */
  /* seg->rl1 = 0; */

  return seg;
}

void update_seg(seg_t *seg, int h, int o) {
  /* if (seg->lasth != h) { */
  /*   if (seg->lasth == -1) { */
  /*     /\* We are inserting the first path *\/ */
  /*     if (h != 0) { */
  /*       /\* First path is not 0, we have to fill the initial gap of 0s *\/ */
  /*       int64_t tmp_cnt[6]; */
  /*       memset(tmp_cnt, 0, 48); */

  /*       seg->paths_p = 0; */
  /*       if (seg->paths_b + 8 >= seg->paths_c) { */
  /*         fprintf(stderr, "Reallocating (a)\n"); */
  /*         seg->paths = realloc(seg->paths, seg->paths_c * 2); */
  /*         seg->paths_c *= 2; */
  /*       } */
  /*       seg->paths_b = */
  /*           rle_insert(seg->paths, seg->paths_p, 0, h, tmp_cnt, seg->cnts);
   */
  /*       /\* printf("%d/%d: inserted %d 0s (a). New size: %d\n", seg->idx, h,
   * h, */
  /*        *\/ */
  /*       /\*        seg->paths_b); *\/ */
  /*       seg->paths_p = h; */
  /*       seg->cnts[0] = h; */
  /*     } else { */
  /*       seg->paths_p = 0; */
  /*     } */
  /*     seg->rl1 = 1; */
  /*     /\* now, rl is 1 since we know we have seen only 1 path. lasth is h and
   */
  /*      * insertion position is h  *\/ */
  /*   } else { */
  /*     /\* We are not inserting the first path *\/ */
  /*     if (h - seg->lasth == 1) { */
  /*       /\* If it's consecutive to last one, just increase run length of 1s
   */
  /*        *\/ */
  /*       ++seg->rl1; */
  /*     } else { */
  /*       /\* Insert runs of 1s *\/ */
  /*       int64_t tmp_cnt[6]; */
  /*       memset(tmp_cnt, 0, 48); */
  /*       if (seg->paths_b + 16 >= seg->paths_c) { */
  /*         seg->paths = realloc(seg->paths, seg->paths_c * 2); */
  /*         seg->paths_c *= 2; */
  /*       } */
  /*       seg->paths_b = rle_insert(seg->paths, seg->paths_p, 1, seg->rl1, */
  /*                                 tmp_cnt, seg->cnts); */
  /*       /\* printf("%d/%d: inserted %d 1s (b). New size: %d\n", seg->idx, h,
   * *\/ */
  /*       /\*        seg->rl1, seg->paths_b); *\/ */
  /*       seg->paths_p += seg->rl1; */
  /*       seg->cnts[1] += seg->rl1; */

  /*       /\* Insert runs of 0s *\/ */
  /*       memset(tmp_cnt, 0, 48); */
  /*       int rl = h - seg->paths_p; // FIXME: do we need a -1? */
  /*       seg->paths_b = */
  /*           rle_insert(seg->paths, seg->paths_p, 0, rl, tmp_cnt, seg->cnts);
   */
  /*       /\* printf("%d/%d: inserted %d 0s (b). New size: %d\n", seg->idx, h,
   * rl, */
  /*        *\/ */
  /*       /\*        seg->paths_b); *\/ */
  /*       seg->paths_p += rl; */
  /*       seg->cnts[0] += rl; */

  /*       seg->rl1 = 1; */
  /*     } */
  /*   } */
  /*   kv_push(int, seg->starts, kv_size(seg->pord)); */
  /*   seg->lasth = h; */
  /* } */
  /* kv_push(int, seg->pord, o); */
}

void compress_seg(seg_t *seg) {
  /* /\* Insert last run of 1s *\/ */
  /* int64_t tmp_cnt[6]; */
  /* memset(tmp_cnt, 0, 48); */
  /* seg->paths_b = */
  /*     rle_insert(seg->paths, seg->paths_p, 1, seg->rl1, tmp_cnt, seg->cnts);
   */
  /* seg->paths_p += seg->rl1; */
  /* seg->cnts[1] += seg->rl1; */

  /* seg->cpord = malloc(vsbound32(kv_size(seg->pord))); */
  /* v8enc32((uint *)seg->pord.a, kv_size(seg->pord), seg->cpord); */
  /* kv_destroy(seg->pord); */
  /* seg->pord.a = NULL; */
}

void destroy_seg(seg_t *seg) {
  if (seg->seq != NULL)
    free(seg->seq);
  /* free(seg->paths); */
  /* if (seg->pord.a != NULL) */
  /*   kv_destroy(seg->pord); */
  /* if (seg->cpord != NULL) */
  /*   free(seg->cpord); */
  /* kv_destroy(seg->starts); */
  free(seg);
}
