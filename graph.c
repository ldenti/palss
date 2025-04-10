#include "graph.h"

graph_t *init_graph(char *fn) {
  graph_t *g = malloc(1 * sizeof(graph_t));
  g->fn = fn;

  /* vertices */
  g->nv = 0;
  g->cv = 16384;
  g->vertices = malloc(g->cv * sizeof(seg_t *));

  g->v_map = kh_init(im);

  /* edges */
  /* g->ne = 0; */
  /* g->out_edges = malloc(16384 * sizeof(kvec_t(int))); */
  /* g->oe_c = 16384; */
  /* for (int i = 0; i < g->oe_c; ++i) */
  /*   kv_init(g->out_edges[i]); */
  /* g->in_edges = malloc(16384 * sizeof(kvec_t(int))); */
  /* g->ie_c = 16384; */
  /* for (int i = 0; i < g->ie_c; ++i) */
  /*   kv_init(g->in_edges[i]); */

  /* paths */
  g->np = 0;
  g->cp = 128;
  g->paths = malloc(g->cp * sizeof(path_t *));
  for (int i = 0; i < g->cp; ++i)
    g->paths[i] = init_path(4096); /* XXX: find a better value */
  return g;
}

void destroy_graph(graph_t *g) {
  for (int i = 0; i < g->nv; ++i)
    destroy_seg(g->vertices[i]);
  free(g->vertices);
  kh_destroy(im, g->v_map);

  /* for (int i = 0; i < g->oe_c; ++i) */
  /*   kv_destroy(g->out_edges[i]); */
  /* free(g->out_edges); */
  /* for (int i = 0; i < g->ie_c; ++i) */
  /*   kv_destroy(g->in_edges[i]); */
  /* free(g->in_edges); */

  for (int i = 0; i < g->cp; ++i)
    destroy_path(g->paths[i]);
  free(g->paths);
  free(g);
}

void dump_gfa(graph_t *g, char *fn) {
  FILE *fp = strcmp(fn, "-") ? fopen(fn, "w") : fdopen(fileno(stdout), "w");
  if (fp == 0) {
    fprintf(stderr, "Cannot open output file\n");
    exit(1);
  }
  // H line
  fprintf(fp, "H\tVN:Z:1.1\n");

  // S lines
  seg_t *s;
  for (int i = 0; i < g->nv; ++i) {
    s = g->vertices[i];
    fprintf(fp, "S\t%d\t%s\n", s->idx, s->seq);
  }
  // L lines
  for (int i = 0; i < g->nv; ++i) {
    for (int x = 0; x < kv_size(g->out_edges[i]); ++x) {
      fprintf(fp, "L\t%d\t+\t%d\t+\t0M\n", g->vertices[i]->idx,
              g->vertices[kv_A(g->out_edges[i], x)]->idx);
    }
  }

  // P lines
  char sep = ',';
  for (int i = 0; i < g->np; ++i) {
    fprintf(fp, "P\t%s\t", g->paths[i]->idx);
    sep = ',';
    for (int y = 0; y < g->paths[i]->l; ++y) {
      if (y == g->paths[i]->l - 1)
        sep = '\t';
      fprintf(fp, "%d%c%c", p_get_v(g->paths[i], y),
              p_get_vs(g->paths[i], y) ? '+' : '-', sep);
    }
    fprintf(fp, "*\n");
  }
}

int load_vertices(graph_t *g, int wseq) {
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
      }
      g->vertices[g->nv] = init_seg();
      gfa_parse_S(s.s, g->vertices[g->nv], wseq);

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

int get_iidx(graph_t *g, int v) {
  return kh_value(g->v_map, kh_get(im, g->v_map, v));
}

/* int load_edges(graph_t *g) { */
/*   if (g->nv == 0) */
/*     load_vertices(g, 0); */

/*   /\* allocate more space if needed *\/ */
/*   if (g->nv > g->oe_c) { */
/*     g->out_edges = realloc(g->out_edges, g->nv * sizeof(kvec_t(int))); */
/*     for (int i = g->oe_c; i < g->nv; ++i) */
/*       kv_init(g->out_edges[i]); */
/*     g->oe_c = g->nv; */
/*   } */
/*   if (g->nv > g->ie_c) { */
/*     g->in_edges = realloc(g->in_edges, g->nv * sizeof(kvec_t(int))); */
/*     for (int i = g->ie_c; i < g->nv; ++i) */
/*       kv_init(g->in_edges[i]); */
/*     g->ie_c = g->nv; */
/*   } */

/*   /\* load edges *\/ */
/*   kstring_t s = {0, 0, 0}; */
/*   int dret; */
/*   gzFile fp = gzopen(g->fn, "r"); */
/*   if (fp == 0) */
/*     return 0; */
/*   kstream_t *ks = ks_init(fp); */
/*   int v1 = -1, v2 = -1; */
/*   while (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) >= 0) { */
/*     if (s.s[0] == 'L') { */
/*       gfa_parse_L(s.s, &v1, &v2); */
/*       v1 = get_iidx(g, v1); */
/*       v2 = get_iidx(g, v2); */
/*       kv_push(int, g->out_edges[v1], v2); */
/*       kv_push(int, g->in_edges[v2], v1); */
/*       ++g->ne; */
/*     } */
/*   } */

/*   free(s.s); */
/*   ks_destroy(ks); */
/*   gzclose(fp); */

/*   return 0; */
/* } */

int update_segments(graph_t *g) {
  path_t *path;
  int iv;
  int min_v = 0, bsize = 2000000;
  while (min_v < g->nv) {
    fprintf(stderr, "Updating vertices [%d,%d)\n", min_v, min_v + bsize);
    for (int p = 0; p < g->np; ++p) {
      path = g->paths[p];
      for (int v = 0; v < path->l; ++v) {
        iv = get_iidx(g, vbgetx32(path->cvertices, v) >> 1);
        /* iv = get_iidx(g, path->vertices[v] >> 1); */
        if (iv >= min_v && iv < min_v + bsize)
          update_seg(g->vertices[iv], p, v);
      }
    }
    for (int v = min_v; v < MIN(g->nv, min_v + bsize); ++v) {
      compress_seg(g->vertices[v]);
    }
    min_v += bsize;
  }
  return 0;
}

int load_paths(graph_t *g) {
  if (g->nv == 0)
    load_vertices(g, 0);

  kstring_t s = {0, 0, 0};
  int dret;
  gzFile fp = gzopen(g->fn, "r");
  if (fp == 0) {
    exit(1);
  }
  kstream_t *ks = ks_init(fp);

  while (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) >= 0) {
    if (s.s[0] == 'P' || s.s[0] == 'W') {
      /* current path: g->np */
      if (g->np == g->cp) {
        g->paths = realloc(g->paths, g->cp * 2 * sizeof(path_t *));
        g->cp *= 2;
        for (int i = g->np; i < g->cp; ++i)
          g->paths[i] = init_path(4096); // XXX: find a better value
      }
      if (s.s[0] == 'P')
        gfa_parse_P(s.s, g->paths[g->np]);
      else
        gfa_parse_W(s.s, g->paths[g->np]);
      compress_path(g->paths[g->np]);
      ++g->np;
    }
  }
  free(s.s);
  ks_destroy(ks);
  gzclose(fp);
  return 0;
}

/** === GFA PARSING ======= **/
void gfa_parse_S(char *s, seg_t *ret, int wseq) {
  int i;       // , is_ok = 0;
  char *p, *q; // *seg = 0, *seq = 0, *rest = 0;
  for (i = 0, p = q = s + 2;; ++p) {
    if (*p == 0 || *p == '\t') {
      int c = *p;
      *p = 0;
      if (i == 0) {
        ret->idx = atoi(q);
        // strcpy(ret->idx, q);
      } else if (i == 1) {
        ret->l = p - q;
        if (wseq == 1) {
          if (ret->l > ret->c) {
            char *temp = realloc(ret->seq, (ret->l * 2) * sizeof(char));
            if (temp == NULL) {
              free(ret->seq);
              fprintf(stderr,
                      "Error while reallocating memory for segment %d\n",
                      ret->idx);
              exit(2);
            } else {
              ret->seq = temp;
            }
            ret->c = ret->l * 2;
          }
          strcpy(ret->seq, q);
        }
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

/* void gfa_parse_L(char *s, int *idx1, int *idx2) { */
/*   int i; */
/*   char *p, *q; */
/*   for (i = 0, p = q = s + 2;; ++p) { */
/*     if (*p == 0 || *p == '\t') { */
/*       int c = *p; */
/*       *p = 0; */
/*       if (i == 0) { */
/*         *idx1 = atoi(q); */
/*       } else if (i == 1) { */
/*         /\* if (*q != '+' && *q != '-') *\/ */
/*         /\*   ; /\\* return -2; *\\/ *\/ */
/*         /\* oriv = (*q != '+'); *\/ */
/*       } else if (i == 2) { */
/*         *idx2 = atoi(q); */
/*       } else if (i == 3) { */
/*         /\* if (*q != '+' && *q != '-') *\/ */
/*         /\* ; /\\* return -2; *\\/ *\/ */
/*         /\* oriw = (*q != '+'); *\/ */
/*       } */
/*       ++i, q = p + 1; */
/*       if (c == 0) */
/*         break; */
/*     } */
/*   } */
/*   /\* return 0; *\/ */
/* } */

void gfa_parse_P(char *s, path_t *path) {
  int i;
  char *p, *q, *qq;
  for (i = 0, p = q = s + 2;; ++p) {
    if (*p == 0 || *p == '\t') {
      *p = 0;
      if (i == 0) {
        /* TODO: check for duplicates */
        strcpy(path->idx, q);
      } else if (i == 1) {
        char strand = *(p - 1);
        qq = q;
        for (qq = q;; ++qq) {
          if (*qq == 0 || *qq == ',') {
            int c = *qq;
            *qq = 0;
            strand = *(qq - 1);
            *(qq - 1) = 0;
            p_add_v(path, atoi(q), strand == '+');
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
  int i;
  char *p, *q, *qq;
  for (i = 0, p = q = s + 2;; ++p) {
    if (*p == 0 || *p == '\t') {
      *p = 0;
      if (i < 2) {
        *p = '#';
        ++i;
        continue;
      } else if (i == 2) {
        strcpy(path->idx, q);
      } else if (i == 5) {
        char strand = *q;
        ++q;
        for (qq = q;; ++qq) {
          if (*qq == 0 || *qq == '>' || *qq == '<') {
            int c = *qq;
            *qq = 0;
            p_add_v(path, atoi(q), strand == '>');
            q = qq + 1;
            if (c == 0)
              break;
            strand = c;
          }
        }
        break;
      }
      ++i, q = p + 1;
    }
  }
}

int get_sequence(graph_t *graph, path_t *path, char **pseq, int *pseq_c) {
  int l = 0;
  for (int i = 0; i < path->l; ++i)
    l += graph->vertices[path->vertices[i]]->l;
  if (l + 1 > *pseq_c) {
    char *temp = (char *)realloc(*pseq, (l + 1) * 2 * sizeof(char));
    if (temp == NULL) {
      free(pseq);
      fprintf(stderr, "Error while reallocating memory for path string\n");
      exit(2);
    } else {
      *pseq = temp;
    }
    *pseq_c = (l + 1) * 2;
  }
  int p = 0;
  for (int i = 0; i < path->l; ++i) {
    l = graph->vertices[path->vertices[i]]->l;
    strncpy(*pseq + p, graph->vertices[path->vertices[i]]->seq, l);
    p += l;
  }
  (*pseq)[p] = '\0';

  return p;
}

/**
graph_t *canonicalize(graph_t *g) {
  int total_size = 0;
  for (int v = 0; v < g->nv; ++v)
    total_size += g->vertices[v]->l;

  graph_t *cg = malloc(1 * sizeof(graph_t));
  g->fn = NULL;

  // vertices
  cg->nv = 0;
  cg->cv = total_size;
  cg->vertices = malloc(cg->cv * sizeof(seg_t *));
  for (int i = 0; i < cg->cv; ++i)
    cg->vertices[i] = init_seg();
  cg->v_map = kh_init(im);

  // edges
  cg->ne = 0;
  cg->out_edges = malloc(total_size * sizeof(kvec_t(int)));
  cg->oe_c = total_size;
  for (int i = 0; i < cg->oe_c; ++i)
    kv_init(cg->out_edges[i]);
  cg->in_edges = malloc(total_size * sizeof(kvec_t(int)));
  cg->ie_c = total_size;
  for (int i = 0; i < cg->ie_c; ++i)
    kv_init(cg->in_edges[i]);

  // paths
  cg->np = 0;
  cg->cp = g->np;
  cg->paths = malloc(cg->cp * sizeof(path_t *));
  for (int i = 0; i < cg->cp; ++i)
    cg->paths[i] = init_path(16384); // XXX: find a better value

  // for each vertex of original graph (internal id in [0,n)),
  // store starting/ending canonical vertices in new graph
  int *starting = malloc(g->nv * sizeof(int));
  int *ending = malloc(g->nv * sizeof(int));

  int hret;
  khint_t k;
  seg_t *seg;
  seg_t *new_seg;
  for (int v = 0; v < g->nv; ++v) {
    seg = g->vertices[v];

    starting[v] = cg->nv;
    for (int c = 0; c < seg->l; ++c) {
      new_seg = cg->vertices[cg->nv];
      new_seg->idx = cg->nv;
      // XXX: what if we realloc to 1?
      // ns->seq = realloc(ns->seq, 2 * sizeof(char));
      new_seg->seq[0] = seg->seq[c];
      new_seg->seq[1] = '\0';
      new_seg->l = 1;

      k = kh_put(im, cg->v_map, new_seg->idx, &hret);
      kh_value(cg->v_map, k) = new_seg->idx;

      // internal link
      if (c > 0) {
        kv_push(int, cg->out_edges[cg->nv - 1], cg->nv);
        ++cg->ne;
        kv_push(int, cg->in_edges[cg->nv], cg->nv - 1);
      }
      ++cg->nv;
    }
    ending[v] = cg->nv - 1;
  }

  // links between original vertices
  int vv, xx;
  for (int v = 0; v < g->nv; ++v) {
    vv = ending[v];
    for (int x = 0; x < kv_size(g->out_edges[v]); ++x) {
      xx = starting[kv_A(g->out_edges[v], x)];
      kv_push(int, cg->out_edges[vv], xx);
      ++cg->ne;
      kv_push(int, cg->in_edges[xx], vv);
    }
  }

  // XXX: check this
  path_t *path, *new_path;
  for (int p = 0; p < g->np; ++p) {
    path = g->paths[p];
    new_path = cg->paths[p];
    strcpy(new_path->idx, path->idx); // XXX: check this
    for (int v = 0; v < path->l; ++v) {
      for (int x = starting[v]; x <= ending[v]; ++x)
        add_vertex(new_path, x);
    }
    ++cg->np;
  }
  assert(cg->np == cg->cp);

  free(starting);
  free(ending);

  return cg;
}

graph_t *extract_subgraph(graph_t *g, int v1, int v2) {
  // v1 and v2 are in GFA space
  graph_t *sg = init_graph(NULL);

  // XXX: assuming topological sorting for now
  sg->vertices = realloc(sg->vertices, (v2 - v1 + 1) * sizeof(seg_t *));
  sg->cv = v2 - v1 + 1;
  for (int i = sg->nv; i < sg->cv; ++i)
    sg->vertices[i] = init_seg();

  sg->out_edges = malloc((v2 - v1 + 1) * sizeof(kvec_t(int)));
  sg->oe_c = v2 - v1 + 1;
  for (int i = 0; i < sg->oe_c; ++i)
    kv_init(sg->out_edges[i]);
  sg->in_edges = malloc((v2 - v1 + 1) * sizeof(kvec_t(int)));
  sg->ie_c = v2 - v1 + 1;
  for (int i = 0; i < sg->ie_c; ++i)
    kv_init(sg->in_edges[i]);

  sg->paths = realloc(sg->paths, g->np * sizeof(path_t *));
  sg->cp = g->np;
  // they will be allocated later while extracting the subpaths
  // for (int i = sg->np; i < sg->cp; ++i)
  //   sg->paths[i] = init_path(16384); // XXX: find a better value

  // vertices
  int hret;
  khint_t k;
  for (int v = v1; v <= v2; ++v) {
    seg_t *s = get_vertex(g, v);
    seg_t *ns = sg->vertices[sg->nv];
    ns->idx = s->idx;
    ns->seq = realloc(ns->seq, (s->l + 1) * sizeof(char));
    strncpy(ns->seq, s->seq, s->l + 1);
    ns->l = s->l;
    ns->c = ns->l + 1;

    k = kh_put(im, sg->v_map, ns->idx, &hret);
    kh_value(sg->v_map, k) = sg->nv;

    // printf("%d:%d\n", ns->idx, sg->nv);
    ++sg->nv;
  }

  // edges
  int iv1, siv1;
  int iv2, siv2;
  for (int v = v1; v <= v2; ++v) {
    iv1 = kh_value(g->v_map, kh_get(im, g->v_map, v));
    siv1 = kh_value(sg->v_map, kh_get(im, sg->v_map, v));

    if (v < v2) {
      for (int x = 0; x < kv_size(g->out_edges[iv1]); ++x) {
        iv2 = kv_A(g->out_edges[iv1], x);
        siv2 =
            kh_value(sg->v_map, kh_get(im, sg->v_map, g->vertices[iv2]->idx));
        kv_push(int, sg->out_edges[siv1], siv2);
        ++sg->ne;
      }
    }
    if (v > v1) {
      for (int x = 0; x < kv_size(g->in_edges[iv1]); ++x) {
        iv2 = kv_A(g->in_edges[iv1], x);
        siv2 =
            kh_value(sg->v_map, kh_get(im, sg->v_map, g->vertices[iv2]->idx));
        kv_push(int, sg->in_edges[siv2], siv1);
      }
    }
  }

  // paths
  for (int i = 0; i < g->np; ++i) {
    sg->paths[i] = extract_subpath(g, g->paths[i], v1, v2);
    ++sg->np;
  }
  assert(sg->np == sg->cp);

  // FIXME: we may end up with vertices/links not in the paths

  return sg;
}

seg_t *get_vertex(graph_t *g, int v) {
  // XXX: check if key is not in map
  return g->vertices[get_iidx(g, v)];
}

int is_source(graph_t *g, int v) {
  return kv_size(g->in_edges[get_iidx(g, v)]) == 0;
}

int is_sink(graph_t *g, int v) {
  return kv_size(g->out_edges[get_iidx(g, v)]) == 0;
}

int contains(graph_t *g, path_t *p, int v) {
  for (int i = 0; i < p->l; ++i) {
    if (p->vertices[i] == get_iidx(g, v))
      return 1;
  }
  return 0;
}

int get_father(graph_t *g, path_t *p, int v) {
  for (int i = 1; i < p->l; ++i) {
    if (p->vertices[i] == get_iidx(g, v))
      return g->vertices[p->vertices[i - 1]]->idx;
  }
  return -1;
}

// XXX: this should return a string, but I'll use it only for 1bp vertices
char get_label(graph_t *g, int v) {
  return g->vertices[v]->seq[0];
}

int compatible(graph_t *g, int x, int y) {
  // x and y are in GFA space
  if (x > y) {
    int tmp = x;
    x = y;
    y = tmp;
  }
  for (int xi = 0; xi < kv_size(g->vertices[get_iidx(g, x)]->paths); ++xi) {
    for (int yi = 0; yi < kv_size(g->vertices[get_iidx(g, y)]->paths); ++yi) {
      if (kv_A(g->vertices[get_iidx(g, x)]->paths, xi) ==
          kv_A(g->vertices[get_iidx(g, y)]->paths, yi))
        return 1;
    }
  }
  return 0;
}

path_t *extract_subpath(graph_t *g, path_t *path, int x, int y) {
  // x and y are in GFA space
  path_t *subpath = init_path(y - x + 1);
  int f = 0, ok = 0;
  strcpy(subpath->idx, path->idx); // XXX: check this
  for (int i = 0; i < path->l; ++i) {
    if (g->vertices[path->vertices[i]]->idx == x)
      f = 1;
    if (g->vertices[path->vertices[i]]->idx == y)
      ok = f;
    if (f) {
      add_vertex(subpath, path->vertices[i]);
    }
    if (g->vertices[path->vertices[i]]->idx > y || ok)
      break;
  }
  if (!ok) {
    destroy_path(subpath);
    subpath = NULL;
  }
  return subpath;
}
**/
