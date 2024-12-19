#ifndef PSGRAPH_H
#define PSGRAPH_H

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "khash.h"
#include "kseq.h"
#include "kvec.h"

KSTREAM_INIT(gzFile, gzread, 65536)
KHASH_MAP_INIT_INT(im, int)

/* Some assumptions I made:
 * - graph is chopped (size of segments < 4096)
 * - graph has less than INT_MAX nodes
 * - graph is topological sorted
 */

typedef struct {
  int idx;   // identifier (as in gfa) - assuming integer
  char *seq; // actual sequence
  int l;     // actual length
  int c;     // capacity
  kvec_t(int) paths;
} seg_t;

typedef struct {
  int v1;
  int v2;
} link_t;

typedef struct path_t {
  char *idx;     // identifier (as in gfa)
  int *vertices; // identifiers of the vertices along the path (graph space)
  int l;         // actual length
  int capacity;  // capacity
} path_t;

typedef struct graph_t {
  char *fn; // file name

  int nv;           // number of vertices
  int cv;           // allocated vertices
  seg_t **vertices; // vertices
  khash_t(im) *
      v_map; // mapping between gfa idx and internal idx (position on vertices)

  int ne;                  // total number of edges
  kvec_t(int) * out_edges; // outgoing edges
  int oe_c;
  kvec_t(int) * in_edges; // incoming edges
  int ie_c;

  int np;         // number of paths
  int cp;         // allocated paths
  path_t **paths; // paths
} graph_t;

/* Initialize a graph */
graph_t *init_graph(char *fn);

/* Dump graph in GFA format*/
void dump_gfa(graph_t *g, char *fn);

/* Destroy the graph */
void destroy_graph(graph_t *g);

/* Load all vertices from gfa */
int load_vertices(graph_t *g);

/* Load all edges from gfa */
int load_edges(graph_t *g);

/* Load all paths from gfa */
int load_paths(graph_t *g);

/* Extract subgraph in between vertices v1 and v2 */
// XXX: assuming topological sorting
graph_t *extract_subgraph(graph_t *g, int v1, int v2);

/* Make graph canonical */
graph_t *canonicalize(graph_t *g);

/* Get internal identifier from GFA identifier */
int get_iidx(graph_t *g, int v);

/* Return the segment given the gfa idx */
seg_t *get_vertex(graph_t *g, int v);

int is_source(graph_t *g, int v);
int is_sink(graph_t *g, int v);
int contains(graph_t *g, path_t *p, int v);
int get_father(graph_t *g, path_t *p, int v);
char get_label(graph_t *g, int v);

/* Extract subpath */
// XXX: this function init the returned path. Maybe we should init it outside
path_t *extract_subpath(graph_t *g, path_t *path, int x, int y);

/* Check if x and y are on at least one path */
int compatible(graph_t *g, int x, int y);

/* Initialize a segment */
seg_t *init_seg();

/* Destroy the segment */
void destroy_seg(seg_t *seg);

/* Initialize a path with capacity c*/
path_t *init_path(int c);

/* Add vertex v to path, reallocating it if needed */
void add_vertex(path_t *path, int v);

/* Get path sequence */
int get_sequence(graph_t *graph, path_t *path, char **pseq, int *pseq_c);

/* Destroy a path */
void destroy_path(path_t *p);

/* GFA reading utilities */
void gfa_parse_S(char *s, seg_t *ret);
void gfa_parse_L(char *s, int *idx1, int *idx2);
void gfa_parse_P(char *s, path_t *path);
void gfa_parse_W(char *s, path_t *path);

#endif
