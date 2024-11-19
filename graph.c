#include "graph.h"

// Init a path with initial capacity c
path_t *init_path(int c) {
  path_t *path = (path_t *)malloc(1 * sizeof(path_t));
  path->idx = (char *)malloc(128 * sizeof(char));
  path->vertices = (int *)malloc(c * sizeof(int));
  path->l = 0;
  path->capacity = c;
  return path;
}

// Add a vertex to the path, reallocating if needed
void add_vertex(path_t *path, int v) {
  if (path->l >= path->capacity) {
    path->vertices =
        (int *)realloc(path->vertices, path->capacity * 2 * sizeof(int));
    path->capacity *= 2;
  }
  path->vertices[path->l] = v;
  ++path->l;
}

// Destroy the path
void destroy_path(path_t *p) {
  free(p->idx);
  free(p->vertices);
  free(p);
}

seg_t *init_seg() {
  seg_t *seg = malloc(1 * sizeof(seg_t));
  seg->l = 0;
  // seg->idx = (char *)malloc(1024);
  seg->seq = malloc(4096 * sizeof(char));
  seg->c = 4096;
  return seg;
}
void destroy_seg(seg_t *seg) {
  free(seg->seq);
  free(seg);
}

void gfa_parse_S(char *s, seg_t *ret) {
  int i, is_ok = 0;
  char *p, *q, *seg = 0, *seq = 0, *rest = 0;
  uint32_t sid, len = 0;
  for (i = 0, p = q = s + 2;; ++p) {
    if (*p == 0 || *p == '\t') {
      int c = *p;
      *p = 0;
      if (i == 0) {
        ret->idx = atoi(q);
        // strcpy(ret->idx, q);
      } else if (i == 1) {
        // TODO: reallocate if vertex is longer than 4096
        // right now we assume to have a vg chopped graph
        strcpy(ret->seq, q);
        ret->l = p - q;
        is_ok = 1, rest = c ? p + 1 : 0;
        break;
      }
      ++i, q = p + 1;
      if (c == 0)
        break;
    }
  }
  // if (!is_ok) { // something is missing
}

void gfa_parse_P(char *s, path_t *path) {
  int x = 0;        // current index for insertion
  int i;            // , oriv = -1, oriw = -1, is_ok = 0;
  char *p, *q, *qq; //, *segv = 0, *segw = 0, *rest = 0;
  // int32_t ov = INT32_MAX, ow = INT32_MAX;
  for (i = 0, p = q = s + 2;; ++p) {
    if (*p == 0 || *p == '\t') {
      *p = 0;
      if (i == 0) {
        // TODO: check for duplicates
        strcpy(path->idx, q);
      } else if (i == 1) {
        char strand = *(p - 1);
        qq = q;
        for (qq = q;; ++qq) {
          if (*qq == 0 || *qq == ',') {
            int c = *qq;
            *qq = 0;
            if (*(qq - 1) != strand) {
              fprintf(stderr, "Mixed +/- strands in path %s. Aborting...\n",
                      path->idx);
              exit(1);
            }
            *(qq - 1) = 0;
            /*cerr << "Adding " << q << " to " << path->idx << endl;*/
            add_vertex(path, atoi(q));
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
  int x = 0; // current index for insertion
  int i;
  char *p, *q, *qq;
  for (i = 0, p = q = s + 2;; ++p) {
    if (*p == 0 || *p == '\t') {
      *p = 0;
      if (i == 0) {
        strcpy(path->idx, q);
        *(p + 2) = 0;
        path->idx[p - q] = '#';
        path->idx[p - q + 1] = *(p + 1);
        path->idx[p - q + 2] = '\0';
      } else if (i == 5) {
        ++q;
        for (qq = q;; ++qq) {
          if (*qq == 0 || *qq == '>') {
            int c = *qq;
            *qq = 0;
            add_vertex(path, atoi(q));
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
