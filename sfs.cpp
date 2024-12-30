#include "sfs.h"

anchor_t parse_anchor(char *line) {
  anchor_t a;
  int i;
  char *p, *q;
  for (i = 0, p = q = line;; ++p) {
    if (*p == 0 || *p == ':') {
      int c = *p;
      *p = 0;
      /* 0: vertex (as in GFA)
         1: offset
         2: kmer
       */
      if (i == 0) {
        a.v = atoi(q);
      } else if (i == 1) {
        a.offset = atoi(q);
      } else if (i == 2) {
        // char *pEnd;
        a.seq = strtoull(q, NULL /*&pEnd*/, 10);
      }
      ++i;
      q = p + 1;
      if (c == 0)
        break;
    }
  }
  return a;
}

sfs_t parse_sfs_line(char *line) {
  sfs_t s;
  int i;
  char *p, *q;
  for (i = 0, p = q = line;; ++p) {
    if (*p == 0 || *p == ' ') {
      int c = *p;
      *p = 0;
      /* 0: read name
         1: start on read
         2: length
         3: strand
         4: if we kept it
         5: sequence
         6: left anchor
         7: right anchor
       */
      if (i == 0) {
        s.rname = (char *)malloc(p - q + 1);
        strncpy(s.rname, q, p - q);
        s.rname[p - q] = '\0';
      } else if (i == 1) {
        s.s = atoi(q);
      } else if (i == 2) {
        s.l = atoi(q);
      } else if (i == 5) {
        s.seq = (uint8_t *)malloc(s.l + 1);
        memcpy(s.seq, q, s.l);
        s.seq[s.l] = '\0';
        for (int _i = 0; _i < s.l; ++_i)
          s.seq[_i] = s.seq[_i] < 128 ? to_int[s.seq[_i]] : 5;
      } else if (i == 6) {
        s.a = parse_anchor(q);
        s.a.p = s.s;
      } else if (i == 7) {
        s.b = parse_anchor(q);
        s.b.p = s.s;
      }
      ++i;
      q = p + 1;
      if (c == 0)
        break;
    }
  }
  return s;
}
