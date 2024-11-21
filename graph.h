#ifndef PSGRAPH_HPP
#define PSGRAPH_HPP

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "khash.h"
#include "kseq.h"

KSTREAM_INIT(gzFile, gzread, 65536)
KHASH_MAP_INIT_INT(im, int)

/* Some assumptions I made:
 * - graph is chopped (size of segments < 4096)
 * - graph has less than INT_MAX nodes
 */

typedef struct {
  int idx;   // identifier (as in gfa) - assuming integer
  char *seq; // actual sequence
  int l;     // actual length
  int c;     // capacity
} seg_t;

// typedef struct {
//   int idx1;
//   int idx2;
// } link_t;

typedef struct path_t {
  char *idx;     // identifier (as in gfa)
  int *vertices; // identifiers of the vertices along the path
  int l;         // actual length
  int capacity;  // capacity
} path_t;

typedef struct graph_t {
  char *fn;         // file name
  int nv;           // number of vertices
  int cv;           // allocated vertices
  seg_t **vertices; // vertices
  khash_t(im) *
      v_map; // mapping between gfa idx and internal idx (position on vertices)
  int np;    // number of paths
  int cp;    // allocated paths
  path_t **paths; // paths
} graph_t;

/* Initialize a graph */
graph_t *init_graph(char *fn);

/* Destroy the graph */
void destroy_graph(graph_t *);

/* Load all vertices from gfa */
int load_vertices(graph_t *);

/* Load all paths from gfa */
int load_paths(graph_t *);

/* Return the segment given the gfa idx */
seg_t *get_vertex(graph_t *, int);

// int compatible(int x, int y);

/* Initialize a segment */
seg_t *init_seg();

/* Destroy the segment */
void destroy_seg(seg_t *seg);

/* Initialize a path with capacity c*/
path_t *init_path(int c);

/* Add vertex v to path */
void add_vertex(path_t *path, int v);

/* Destroy a path */
void destroy_path(path_t *p);

/* GFA reading utilities */
void gfa_parse_S(char *s, seg_t *ret);
void gfa_parse_P(char *s, path_t *path);
void gfa_parse_W(char *s, path_t *path);

#endif
