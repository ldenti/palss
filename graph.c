#include "graph.h"

graph_t *init_graph(char *fn) {
  graph_t *g = malloc(1 * sizeof(graph_t));
  g->fn = fn;
  // vertices
  g->nv = 0;
  g->cv = 16384;
  g->vertices = malloc(g->cv * sizeof(seg_t *));
  for (int i = 0; i < g->cv; ++i)
    g->vertices[i] = init_seg();
  g->v_map = kh_init(im);

  // paths
  g->np = 0;
  g->cp = 128;
  g->paths = malloc(g->cp * sizeof(path_t *));
  for (int i = 0; i < g->cp; ++i)
    g->paths[i] = init_path(16384); // XXX: find a better value
  return g;
}

void destroy_graph(graph_t *g) {
  for (int i = 0; i < g->cv; ++i)
    destroy_seg(g->vertices[i]);
  free(g->vertices);
  kh_destroy(im, g->v_map);
  for (int i = 0; i < g->cp; ++i)
    destroy_path(g->paths[i]);
  free(g->paths);
  free(g);
}

int load_vertices(graph_t *g) {
  kstring_t s = {0, 0, 0};
  int dret;
  gzFile fp = gzopen(g->fn, "r");
  if (fp == 0)
    return 0;
  kstream_t *ks = ks_init(fp);
  int hret;
  khint_t k;
  while (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) >= 0) {
    if (s.s[0] == 'S') {
      if (g->nv == g->cv) {
        g->vertices = realloc(g->vertices, g->cv * 2 * sizeof(seg_t *));
        g->cv *= 2;
        for (int64_t i = g->nv; i < g->cv; ++i)
          g->vertices[i] = init_seg();
      }
      gfa_parse_S(s.s, g->vertices[g->nv]);

      k = kh_put(im, g->v_map, g->vertices[g->nv]->idx, &hret);
      kh_value(g->v_map, k) = g->nv;
      ++g->nv;
    }
  }
  free(s.s);
  ks_destroy(ks);
  gzclose(fp);
  return 0;
}

seg_t *get_vertex(graph_t *g, int idx) {
  // XXX: check if key is not in map
  return g->vertices[kh_value(g->v_map, idx)];
}

int compatible(graph_t *g, int x, int y) {
  if (x > y) {
    int tmp = x;
    x = y;
    y = tmp;
  }
  for (int p = 0; p < g->np; ++p) {
    int f = 0; //, ok = 0;
    for (int i = 0; i < g->paths[p]->l; ++i) {
      if (g->paths[p]->vertices[i] == x)
        f = 1;
      if (f && g->paths[p]->vertices[i] == y)
        return f;
    }
  }
  return 0;
}

// Init a path with initial capacity c
path_t *init_path(int c) {
  path_t *path = malloc(1 * sizeof(path_t));
  path->idx = malloc(128 * sizeof(char));
  path->vertices = malloc(c * sizeof(int));
  path->l = 0;
  path->capacity = c;
  return path;
}

// Add a vertex to the path, reallocating if needed
void add_vertex(path_t *path, int v) {
  if (path->l >= path->capacity) {
    path->vertices =
        (int *)realloc(path->vertices, path->capacity * 2 * sizeof(int));
    path->capacity *= 2;
  }
  path->vertices[path->l] = v;
  ++path->l;
}

/* Load paths from GFA file */
int load_paths(graph_t *g) {
  kstring_t s = {0, 0, 0};
  int dret;
  gzFile fp = gzopen(g->fn, "r");
  if (fp == 0) {
    exit(1); // return paths;
  }
  kstream_t *ks = ks_init(fp);

  while (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) >= 0) {
    if (s.s[0] == 'P' || s.s[0] == 'W') {
      if (g->np == g->cp) {
        g->paths = realloc(g->paths, g->cp * 2 * sizeof(path_t *));
        g->cp *= 2;
        for (int i = g->np; i < g->cp; ++i)
          g->paths[i] = init_path(16384); // XXX: find a better value
      }
      if (s.s[0] == 'P')
        gfa_parse_P(s.s, g->paths[g->np]);
      else
        gfa_parse_W(s.s, g->paths[g->np]);
      ++g->np;
    }
  }
  free(s.s);
  ks_destroy(ks);
  gzclose(fp);
  return 0;
}

// Destroy the path
void destroy_path(path_t *p) {
  free(p->idx);
  free(p->vertices);
  free(p);
}

seg_t *init_seg() {
  seg_t *seg = malloc(1 * sizeof(seg_t));
  seg->l = 0;
  // seg->idx = (char *)malloc(1024);
  seg->seq = malloc(4096 * sizeof(char));
  seg->c = 4096; // XXX: assuming chopped graph
  return seg;
}

void destroy_seg(seg_t *seg) {
  free(seg->seq);
  free(seg);
}

void gfa_parse_S(char *s, seg_t *ret) {
  int i; // , is_ok = 0;
  char *p, *q; // *seg = 0, *seq = 0, *rest = 0;
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
        // is_ok = 1, rest = c ? p + 1 : 0;
        break;
      }
      ++i, q = p + 1;
      if (c == 0)
        break;
    }
  }
  // if (!is_ok) { // something is missing
}

void gfa_parse_P(char *s, path_t *path) {
  // int x = 0; // current index for insertion
  int i;
  char *p, *q, *qq;
  for (i = 0, p = q = s + 2;; ++p) {
    if (*p == 0 || *p == '\t') {
      *p = 0;
      if (i == 0) {
        // TODO: check for duplicates
        strcpy(path->idx, q);
      } else if (i == 1) {
        char strand = *(p - 1);
        qq = q;
        for (qq = q;; ++qq) {
          if (*qq == 0 || *qq == ',') {
            int c = *qq;
            *qq = 0;
            if (*(qq - 1) != strand) {
              fprintf(stderr, "Mixed +/- strands in path %s. Aborting...\n",
                      path->idx);
              exit(1);
            }
            *(qq - 1) = 0;
            add_vertex(path, atoi(q));
            q = qq + 1;
            if (c == 0)
              break;
          }
        }
        break;
      }
      ++i, q = p + 1;
    }
  }
}

void gfa_parse_W(char *s, path_t *path) {
  // int x = 0; // current index for insertion
  int i;
  char *p, *q, *qq;
  for (i = 0, p = q = s + 2;; ++p) {
    if (*p == 0 || *p == '\t') {
      *p = 0;
      if (i == 0) {
        strcpy(path->idx, q);
        *(p + 2) = 0;
        path->idx[p - q] = '#';
        path->idx[p - q + 1] = *(p + 1);
        path->idx[p - q + 2] = '\0';
      } else if (i == 5) {
        ++q;
        for (qq = q;; ++qq) {
          if (*qq == 0 || *qq == '>') {
            int c = *qq;
            *qq = 0;
            add_vertex(path, atoi(q));
            q = qq + 1;
            if (c == 0)
              break;
          }
        }
        break;
      }
      ++i, q = p + 1;
    }
  }
}
