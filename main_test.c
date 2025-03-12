#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <unistd.h>

#include "graph.h"
#include "path.h"
#include "segment.h"

#include "rle.h"

static inline double realtime() {
  struct timeval tp;
  struct timezone tzp;
  gettimeofday(&tp, &tzp);
  return (double)tp.tv_sec + (double)tp.tv_usec * 1e-6;
}

int main_test(int argc, char *argv[]) {
  /**
  uint8_t *block = calloc(16, 1);
  int64_t tmp_cnt[6];
  memset(tmp_cnt, 0, 48);
  int64_t cnt[6];
  memset(cnt, 0, 48);

  int b = 0;
  int p = 0;
  b = rle_insert(block, p, 0, 2, tmp_cnt, cnt);
  p += 3;
  cnt[0] += 3;
  b = rle_insert(block, p, 1, 23, tmp_cnt, cnt);
  p += 3;
  cnt[1] += 3;
  rle_print(block, 1);

  printf("%d\n", b);

  free(block);
  **/

  char *gfa_fn = argv[1];

  double rt = realtime();
  graph_t *graph = init_graph(gfa_fn);
  load_vertices(graph, 1);
  fprintf(stderr, "loaded %d vertices in %.3f secs\n", graph->nv,
          realtime() - rt);
  rt = realtime();
  load_paths(graph);
  fprintf(stderr, "loaded %d paths in %.3f secs\n", graph->np, realtime() - rt);

  /* for (int v = 0; v < graph->nv; ++v) { */
  /*   printf("%d %d %s %d %d %d\n", v, graph->vertices[v]->idx, */
  /*          graph->vertices[v]->seq, *rle_nptr(graph->vertices[v]->paths), */
  /*          graph->vertices[v]->paths_b, kv_size(graph->vertices[v]->starts));
   */

  /*   rle_print(graph->vertices[v]->paths, 0); */
  /* } */
  /* dump_gfa(graph, "-"); */

  /* int vidx; */
  /* path_t *path; */
  /* for (int v = 0; v < graph->nv; ++v) { */
  /*   vidx = graph->vertices[v]->idx; */
  /*   printf("=== %d (%d) ===\n", vidx, v); */
  /*   for (int p = 0; p < graph->np; ++p) { */
  /*     path = graph->paths[p]; */
  /*     int *ords = p_get_vord(path, vidx); */
  /*     int n = p_get_vocc(path, vidx); */
  /*     printf("%d %s %d (%d)", p, graph->paths[p]->idx, p_has_v(path, vidx),
   * n); */
  /*     for (int i = 0; i < n; ++i) { */
  /*       printf(" %d%c", ords[i], p_get_vs(path, ords[i]) ? '+' : '-'); */
  /*     } */
  /*     printf("\n"); */
  /*     free(ords); */
  /*   } */
  /* } */

  /* /\* path_t *path; *\/ */
  /* for (int p = 0; p < graph->np; ++p) { */
  /*   compress_path(graph->paths[p]); */
  /*   /\* printf("%d %d %f %f\n", path->l, path->capacity, *\/ */
  /*   /\*        (float)path->l * 32 / 8 / 1024, *\/ */
  /*   /\*        (float)path->capacity * 32 / 8 / 1024); *\/ */
  /* } */

  destroy_graph(graph);

  return 0;
}
