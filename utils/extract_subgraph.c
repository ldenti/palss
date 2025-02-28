#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "graph.h"
#include "misc.h"
#include "path.h"
#include "segment.h"

int parse_region(char *region, char *chr, int *s, int *e) {
  int i;
  char *p, *q;
  for (i = 0, p = q = region;; ++p) {
    if (*p == ':' || *p == '-' || *p == 0) {
      int c = *p;
      *p = 0;
      if (i == 0 && c == ':') {
        strncpy(chr, q, p - q);
        chr[p - q] = '\0';
      } else if (i == 1 && c == '-') {
        *s = atoi(q);
      } else if (i == 2 && c == 0) {
        *e = atoi(q);
        return 0;
      } else {
        return 1;
      }
      ++i;
      q = p + 1;
    }
  }
}

int main(int argc, char *argv[]) {
  char *gfa_fn = argv[1];
  char *region = argv[2];

  double rt = realtime();
  graph_t *graph = init_graph(gfa_fn);
  load_vertices(graph, 0);
  fprintf(stderr, "loaded %d vertices in %.3f secs\n", graph->nv,
          realtime() - rt);
  rt = realtime();
  load_paths(graph);
  fprintf(stderr, "loaded %d paths in %.3f secs\n", graph->np, realtime() - rt);

  char *chr = malloc((strlen(region) + 1) * sizeof(char));
  int s, e;
  if (parse_region(region, chr, &s, &e)) {
    fprintf(stderr, "Wrong region\n");
    return 1;
  }

  int p;
  for (p = 0; p < graph->np; ++p) {
    if (strcmp(graph->paths[p]->idx, chr) == 0)
      break;
  }
  if (p == graph->np) {
    fprintf(stderr, "Path not found\n");
    return 1;
  }

  path_t *path = graph->paths[p];
  int offset = 0;
  int sv = -1, ev = -1;
  for (int v = 0; v < path->l; ++v) {
    int piv = p_get_v(path, v);
    int giv = get_iidx(graph, piv);
    offset += graph->vertices[giv]->l;
    if (sv == -1 && offset >= s)
      sv = graph->vertices[giv]->idx;
    if (sv != -1 && ev == -1 && offset >= e) {
      ev = graph->vertices[giv]->idx;
      break;
    }
  }
  if (sv == -1 || ev == -1) {
    fprintf(stderr, "Something is wrong\n");
    return 1;
  }

  fprintf(stderr, "Extracting subgraph %d>%d\n", sv, ev);
  int *ps, *pe;
  for (p = 0; p < graph->np; ++p) {
    path = graph->paths[p];
    ps = p_get_pv(path, sv, 1);
    pe = p_get_pv(path, ev, 0);
    if (ps == NULL || pe == NULL)
      continue;
    while (ps <= pe) {
      int piv = *ps >> 1;
      int giv = get_iidx(graph, piv);
      printf("%d %d\n", graph->vertices[giv]->idx, graph->vertices[giv]->l);
      ++ps;
    }
  }

  free(chr);
  destroy_graph(graph);

  return 0;
}
