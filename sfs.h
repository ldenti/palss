#ifndef PS_SFS_H
#define PS_SFS_H

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* #include "utils.h" */

typedef struct {
  int64_t v;    // vertex on graph
  int offset;   // offset on vertex
  int p;        // position on query
  uint64_t seq; // kmer
  int closest;  // is this anchor the closest one to the specific string?
} anchor_t;

typedef struct {
  int qidx;     // read index
  int s;        // start on query
  int l;        // length
  anchor_t *a;  // left anchor
  anchor_t *b;  // right anchor
  int strand;   // inferred strand
  uint64_t esk; // expected starting kmer (from cluster)
  uint64_t eek; // expected ending kmer (from cluster)
  int good;     // is it good for calling step?
  char *rname;  // plain read name
  uint8_t *seq; // 1-4encoded sequence
} sfs_t;

sfs_t *init_sfs();
void destroy_sfs(sfs_t *s);

/* anchor_t parse_anchor(char *line); */
/* void parse_sfs_line(char *line, sfs_t *s); */

#endif
