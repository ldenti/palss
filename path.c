#include "path.h"

path_t *init_path(int c) {
  path_t *path = malloc(1 * sizeof(path_t));
  path->idx = malloc(128 * sizeof(char));
  path->vertices = malloc(c * sizeof(int));
  path->ord = kh_init(im);
  path->occ = kh_init(im);
  path->l = 0;
  path->capacity = c;
  return path;
}

void destroy_path(path_t *path) {
  free(path->idx);
  free(path->vertices);
  kh_destroy(im, path->ord);
  kh_destroy(im, path->occ);
  free(path);
}

void p_add_v(path_t *path, int v, int strand) {
  if (path->l >= path->capacity) {
    path->vertices =
        (int *)realloc(path->vertices, path->capacity * 2 * sizeof(int));
    path->capacity *= 2;
  }
  path->vertices[path->l] = (v << 1) | strand;

  int hret;
  khint_t k;

  k = kh_put(im, path->occ, v, &hret);
  assert(hret < 2);
  if (hret == 1) {
    kh_value(path->occ, k) = 0;
  }
  int ev = p_encode(v, kh_value(path->occ, k));
  ++kh_value(path->occ, k);

  k = kh_put(im, path->ord, ev, &hret);
  assert(hret == 1);
  kh_value(path->ord, k) = path->l;

  ++path->l;
}

int p_has_v(path_t *path, int v) {
  return kh_get(im, path->occ, v) != kh_end(path->occ);
}

int p_get_v(path_t *path, int i) { return path->vertices[i] >> 1; }

int p_get_vs(path_t *path, int i) { return path->vertices[i] & 1; }

int p_get_vocc(path_t *path, int v) {
  return kh_value(path->occ, kh_get(im, path->occ, v));
}

int *p_get_vords(path_t *path, int v) {
  int n = p_get_vocc(path, v);
  int *r = malloc(n * sizeof(int));
  int ev;
  for (int i = 0; i < n; ++i) {
    ev = p_encode(v, i);
    r[i] = kh_value(path->ord, kh_get(im, path->ord, ev));
  }
  return r;
}

int p_get_vord(path_t *path, int v, int first) {
  int x = first ? 1 : p_get_vocc(path, v);
  int ev = p_encode(v, x - 1);
  return kh_value(path->ord, kh_get(im, path->ord, ev));
}

int *p_get_pv(path_t *path, int v, int first) {
  return p_has_v(path, v) ? path->vertices + p_get_vord(path, v, first) : NULL;
}
