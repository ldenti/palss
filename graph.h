#ifndef PS_GRAPH_H
#define PS_GRAPH_H

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "khash.h"
#include "kseq.h"
#include "kvec.h"

#include "misc.h"
#include "path.h"
#include "segment.h"

KSTREAM_INIT(gzFile, gzread, 65536)
/* INT2INT map is already init in path.h */
/* KHASH_MAP_INIT_INT(im, int) */

/* Some assumptions I made:
 * - vertex identifiers in GFA are integers
 * - graph has less than INT_MAX nodes
 */

typedef struct {
  char *fn; // file name
  int nv;   // number of vertices
  int cv;   // allocated vertices

  char *labels;
  int labels_c;
  int labels_n;
  int *labels_ofx;
  int labels_ofx_c;

  seg_t **vertices; // vertices
  khash_t(im) *
      v_map; // mapping between gfa idx and internal idx (position on vertices)

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
int load_vertices(graph_t *g, int wseq);

/* Load all edges from gfa */
int load_edges(graph_t *g);

/* Load all paths from gfa */
int load_paths(graph_t *g);

/* Get internal identifier from GFA identifier */
int get_iidx(graph_t *g, int v);

/* Update vertices with path information */
int update_segments(graph_t *g);

/* Extract subgraph in between vertices v1 and v2 */
// XXX: assuming topological sorting
/* graph_t *extract_subgraph(graph_t *g, int v1, int v2); */

/* Make graph canonical */
/* graph_t *canonicalize(graph_t *g); */

/* Return the segment given the gfa idx */
/* seg_t *get_vertex(graph_t *g, int v); */

/* int is_source(graph_t *g, int v); */
/* int is_sink(graph_t *g, int v); */
/* int contains(graph_t *g, path_t *p, int v); */
/* int get_father(graph_t *g, path_t *p, int v); */
/* char get_label(graph_t *g, int v); */

/* Extract subpath */
// XXX: this function init the returned path. Maybe we should init it
// outside
/* path_t *extract_subpath(graph_t *g, path_t *path, int x, int y); */

/* Check if x and y are on at least one path */
/* int compatible(graph_t *g, int x, int y); */

/* Get path sequence */
/* int get_sequence(graph_t *graph, path_t *path, char **pseq, int *pseq_c);
 */

/* GFA reading utilities */
void gfa_parse_S(char *s, seg_t *ret, int wseq);
void gfa_parse_L(char *s, int *idx1, int *idx2);
void gfa_parse_P(char *s, path_t *path);
void gfa_parse_W(char *s, path_t *path);

#endif
