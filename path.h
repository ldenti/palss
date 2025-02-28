#ifndef PS_PATH_H
#define PS_PATH_H

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "khash.h"

KHASH_MAP_INIT_INT(im, int)

static inline int p_encode(int v, int n) { return (v << 4) | (n & 0xF); }

typedef struct {
  char *idx;     /* path identifier (as in gfa) */
  int *vertices; /* identifiers of the vertices in path order:  31 bits for id
                  * in GFA space, 1 bit for strand */
  khash_t(im) *
      occ; /* number of times each vertex occur, key is just the vertex */
  khash_t(im) *
      ord; /* vertex ordering, 28 bits for vertex id 4 bits for cardinality */
  int l;   /* actual length */
  int capacity; /* capacity */
} path_t;

/* Initialize a path with capacity c*/
path_t *init_path(int c);

/* Destroy a path */
void destroy_path(path_t *path);

/* Add vertex v to path, reallocating it if needed */
void p_add_v(path_t *path, int v, int strand);

/* Check if vertex v (GFA space) is in the graph */
int p_has_v(path_t *path, int v);

/* Return i-th vertex in the path */
int p_get_v(path_t *path, int i);

/* Return strand of i-th vertex in the path */
int p_get_vs(path_t *path, int i);

/* Return #occurences of vertex v */
int p_get_vocc(path_t *path, int v);

/* Return order of vertex v in path ordering (array since we may have cycles).
 * Array must be freed by caller */
int *p_get_vords(path_t *path, int v);

/* Return first (or last) order of vertex v in path ordering */
int p_get_vord(path_t *path, int v, int first);

/* Return pointer to first (or last) occurrence of v in vertices */
int *p_get_pv(path_t *path, int v, int first);

#endif
