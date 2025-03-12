#include "path.h"

path_t *init_path(int c) {
  path_t *path = malloc(1 * sizeof(path_t));
  path->idx = malloc(128 * sizeof(char));
  path->vertices = malloc(c * sizeof(uint));
  path->cvertices = NULL;
  path->l = 0;
  path->capacity = c;
  return path;
}

void compress_path(path_t *path) {
  path->cvertices = malloc(vsbound32(path->l));
  v8enc32(path->vertices, path->l, path->cvertices);
  free(path->vertices);
  path->vertices = NULL;
}

void clear_path(path_t *path) { path->l = 0; }

void destroy_path(path_t *path) {
  free(path->idx);
  if (path->vertices != NULL)
    free(path->vertices);
  if (path->cvertices != NULL)
    free(path->cvertices);
  free(path);
}

void p_add_v(path_t *path, int v, int strand) {
  if (path->l >= path->capacity) {
    path->vertices = realloc(path->vertices, path->capacity * 2 * sizeof(uint));
    path->capacity *= 2;
  }
  path->vertices[path->l] = (v << 1) | strand;
  ++path->l;
}

/* int p_has_v(path_t *path, int v) { */
/*   return kh_get(im, path->occ, v) != kh_end(path->occ); */
/* } */

int p_get_v(path_t *path, int i) { return path->vertices[i] >> 1; }

int p_get_vs(path_t *path, int i) { return path->vertices[i] & 1; }

/* int p_get_vocc(path_t *path, int v) { */
/*   return kh_value(path->occ, kh_get(im, path->occ, v)); */
/* } */

/* int *p_get_vords(path_t *path, int v) { */
/*   int n = p_get_vocc(path, v); */
/*   int *r = malloc(n * sizeof(int)); */
/*   int ev; */
/*   for (int i = 0; i < n; ++i) { */
/*     ev = p_encode(v, i); */
/*     r[i] = kh_value(path->ord, kh_get(im, path->ord, ev)); */
/*   } */
/*   return r; */
/* } */

/* int p_get_vord(path_t *path, int v, int first) { */
/*   int x = first ? 1 : p_get_vocc(path, v); */
/*   int ev = p_encode(v, x - 1); */
/*   return kh_value(path->ord, kh_get(im, path->ord, ev)); */
/* } */

/* int *p_get_pv(path_t *path, int v, int first) { */
/*   return p_has_v(path, v) ? path->vertices + p_get_vord(path, v, first) :
 * NULL; */
/* } */
