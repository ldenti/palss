#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "ksort.h"

#include "graph.h"
#include "misc.h"
#include "path.h"
#include "segment.h"

KSORT_INIT_GENERIC(uint64_t)

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

void get_path(char *gfa_fn, path_t *path, char *idx) {
  clear_path(path);
  kstring_t s = {0, 0, 0};
  int dret;
  gzFile fp = gzopen(gfa_fn, "r");
  if (fp == 0) {
    exit(EXIT_FAILURE);
  }
  kstream_t *ks = ks_init(fp);
  while (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) >= 0) {
    if (s.s[0] == 'P' || s.s[0] == 'W') {
      if (s.s[0] == 'P')
        gfa_parse_P(s.s, path);
      else
        gfa_parse_W(s.s, path);
      if (strcmp(path->idx, idx) == 0) {
        break;
      }
      clear_path(path);
    }
  }
  free(s.s);
  ks_destroy(ks);
  gzclose(fp);
}

int main(int argc, char *argv[]) {
  char *gfa_fn = argv[1];
  char *region = argv[2];

  double rt0 = realtime(), rt = realtime(), rt1;
  graph_t *graph = init_graph(gfa_fn);
  load_vertices(graph, 1);
  fprintf(stderr, "loaded %d vertices in %.3f secs\n", graph->nv,
          realtime() - rt);
  rt = realtime();

  char chr[32];
  int sp, ep;
  path_t *path = init_path(1 << 15);  /* Current "reference" path */
  path_t *path2 = init_path(1 << 15); /* Current path we are visiting */

  parse_region(region, chr, &sp, &ep);

  /* Get the path from GFA if path idx changed */
  get_path(gfa_fn, path, chr);
  if (path->l == 0) {
    fprintf(stderr, "Path %s not found\n", chr);
    exit(EXIT_FAILURE);
  }

  /* Get source/sink vertices by iterating over the path */
  rt = realtime();
  int offset = 0;
  int source = -1, sink = -1;
  for (int v = 0; v < path->l; ++v) {
    int piv = p_get_v(path, v);
    int giv = get_iidx(graph, piv);
    offset += graph->vertices[giv]->l;
    if (source == -1 && offset >= sp)
      source = graph->vertices[giv]->idx;
    if (source != -1 && sink == -1 && offset >= ep) {
      sink = graph->vertices[giv]->idx;
      break;
    }
  }
  if (source == -1 || sink == -1) {
    fprintf(stderr, "Something bad happened\n");
    exit(EXIT_FAILURE);
  }
  fprintf(stderr, "found source (%d) and sink (%d) in %.3f secs\n", source,
          sink, realtime() - rt);

  /* Source and sink are in GFA space */
  kvec_t(uint64_t) vertices; // these could be uint32_t but we initialized a
                             // single ksort using uint64_t
  kv_init(vertices);
  kvec_t(uint64_t) edges;
  kv_init(edges);
  int paths_c = 1024;
  int paths_n = 0;
  path_t **paths = malloc(paths_c * sizeof(path_t *));

  int v1, v2; // vertices we may need

  int *ps, *pe; // pointers to source/sink along path
  kstring_t s = {0, 0, 0};
  int dret;
  gzFile fp = gzopen(gfa_fn, "r");
  if (fp == 0) {
    exit(EXIT_FAILURE);
  }
  kstream_t *ks = ks_init(fp);
  int pp = 0;
  rt = realtime();
  while (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) >= 0) {
    if (s.s[0] == 'P' || s.s[0] == 'W') {
      ++pp;
      if (pp % 50000 == 0) {
        rt1 = realtime() - rt;
        fprintf(stderr,
                "Iterated over %d paths in %.3f secs (avg %.3f per path)\r", pp,
                rt1, rt1 / (float)pp);
      }
      if (s.s[0] == 'P')
        gfa_parse_P(s.s, path2);
      else
        gfa_parse_W(s.s, path2);
      /* Path is in GFA space */
      ps = p_get_pv(path2, source, 1);
      pe = p_get_pv(path2, sink, 0);
      if (ps != NULL && pe != NULL && ps <= pe) {
        if (paths_n == paths_c) {
          paths = realloc(paths, paths_c * 2 * sizeof(path_t *));
          paths_c *= 2;
        }
        paths[paths_n] = init_path(pe - ps + 1);
        strcpy(paths[paths_n]->idx, path2->idx); /* XXX: is this always
      safe? */

        v1 = get_iidx(graph, *ps >> 1);
        kv_push(uint64_t, vertices, v1);
        p_add_v(paths[paths_n], *ps >> 1, *ps & 1);
        while (ps < pe) {
          v2 = get_iidx(graph, *(ps + 1) >> 1);
          p_add_v(paths[paths_n], *(ps + 1) >> 1, *(ps + 1) & 1);
          kv_push(uint64_t, vertices, v2);
          kv_push(uint64_t, edges, (uint64_t)v1 << 32 | v2);
          v1 = v2;
          ++ps;
        }
        ++paths_n;
      } else {
        /* fprintf(stderr, "Skipping path %s since source comes after sink
      along
         * the path\n", path2->idx); */
      }
      clear_path(path2);
    }
  }
  fprintf(stderr, "Iterated over %d paths in %.3f secs\n", pp, realtime() - rt);
  free(s.s);
  ks_destroy(ks);
  gzclose(fp);
  destroy_path(path);
  destroy_path(path2);

  fprintf(stderr, "Dumping %ld vertices (at most)\n", kv_size(vertices));
  ks_introsort(uint64_t, kv_size(vertices), vertices.a);
  v1 = -1;
  for (int i = 0; i < kv_size(vertices); ++i) {
    v2 = kv_A(vertices, i);
    if (v2 != v1) {
      printf("S\t%d\t%s\n", graph->vertices[v2]->idx, graph->vertices[v2]->seq);
      v1 = v2;
    }
  }
  kv_destroy(vertices);

  fprintf(stderr, "Dumping %ld edges (at most)\n", kv_size(edges));
  ks_introsort(uint64_t, kv_size(edges), edges.a);
  uint64_t e1 = -1; /* XXX: this choice could be not robust enough */
  uint64_t e2;
  for (int i = 0; i < kv_size(edges); ++i) {
    e2 = kv_A(edges, i);
    if (e2 != e1) {
      v1 = e2 >> 32;
      v2 = e2 & 0xFFFFFFFF;
      printf("L\t%d\t+\t%d\t+\t*\n", graph->vertices[v1]->idx,
             graph->vertices[v2]->idx);
      e1 = e2;
    }
  }
  kv_destroy(edges);

  /* print paths */
  fprintf(stderr, "Dumping %d paths\n", paths_n);
  char end;
  for (int i = 0; i < paths_n; ++i) {
    path = paths[i];
    printf("P\t%s\t", path->idx);
    end = ',';
    for (int j = 0; j < path->l; ++j) {
      if (j == path->l - 1)
        end = '\t';
      printf("%d%c%c",
             path->vertices[j] >>
                 1 /*graph->vertices[path->vertices[j] >> 1]->idx*/,
             (path->vertices[j] & 1) == 1 ? '+' : '-', end);
    }
    printf("*\n");
    destroy_path(path);
  }
  free(paths);
  destroy_graph(graph);

  fprintf(stderr, "Extracted subgraph in %.3f secs\n", realtime() - rt0);

  return 0;
}
