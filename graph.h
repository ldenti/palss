#ifndef PSGRAPH_HPP
#define PSGRAPH_HPP

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct path_t {
  char *idx;
  int *vertices;
  int l;
  int capacity;
} path_t;

path_t *init_path(int c);
void add_vertex(path_t *path, int v);
int compatible(int x, int y);
void destroy_path(path_t *p);

typedef struct {
  int idx;   // identifier (as in gfa)
  char *seq; // sequence
  int l;     // actual length
  int c;     // capacity
} seg_t;

typedef struct {
  int idx1;
  int idx2;
} link_t;

seg_t *init_seg();
void destroy_seg(seg_t *seg);
void gfa_parse_S(char *s, seg_t *ret);
void gfa_parse_P(char *s, path_t *path);
void gfa_parse_W(char *s, path_t *path);

#endif
