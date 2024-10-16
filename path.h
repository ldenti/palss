#ifndef PSPATH_HPP
#define PSPATH_HPP

struct path_t {
  char *idx;
  int *vertices;
  int l;
  int capacity;
};

// init a path with initial capacity c
inline path_t *init_path(int c) {
  path_t *path = (path_t *)malloc(1 * sizeof(path_t));
  path->idx = (char *)malloc(128 * sizeof(char));
  path->vertices = (int *)malloc(c * sizeof(int));
  path->l = 0;
  path->capacity = c;
  return path;
}

// add a vertex to the path
inline void add_vertex(path_t *path, int v) {
  if (path->l >= path->capacity) {
    path->vertices =
        (int *)realloc(path->vertices, path->capacity * 2 * sizeof(int));
    path->capacity *= 2;
  }
  path->vertices[path->l] = v;
  ++path->l;
}

// destroy the path
inline void destroy_path(path_t *p) {
  free(p->idx);
  free(p->vertices);
  free(p);
}

#endif
