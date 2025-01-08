#include <algorithm>
#include <bit>
#include <cstdint>
#include <iostream>
#include <vector>
#include <zlib.h>

#include "abpoa.h"
#include "fm-index.h"
#include "ketopt.h"
#include "kseq.h"
#include "ksw2.h"

#include "utils.h"

int main_sketch(int argc, char *argv[]);
int main_dump(int argc, char *argv[]);
int main_kan(int argc, char *argv[]);
int main_search(int argc, char *argv[]);
int main_map(int argc, char *argv[]);
int main_chreads(int argc, char *argv[]);
int main_call(int argc, char *argv[]);

using namespace std;

int align(int argc, char *argv[]) {
  // https://github.com/lh3/minimap2/blob/69e36299168d739dded1c5549f662793af10da83/options.c#L144
  int sc_mch = 1, sc_mis = -9, gapo = 16, gape = 2, gapo2 = 41, gape2 = 1;
  int zdrop = -1;
  int flag = 0x40;
  static ko_longopt_t longopts[] = {};
  ketopt_t opt = KETOPT_INIT;
  int _c;
  while ((_c = ketopt(&opt, argc, argv, 1, "m:x:o:e:z:g", longopts)) >= 0) {
    if (_c == 'm')
      sc_mch = atoi(opt.arg);
    else if (_c == 'x')
      sc_mis = -atoi(opt.arg);
    else if (_c == 'o')
      gapo = atoi(opt.arg);
    else if (_c == 'e')
      gape = atoi(opt.arg);
    else if (_c == 'z')
      zdrop = atoi(opt.arg);
    else if (_c == 'g')
      flag = 0;
  }

  char *tseq = argv[opt.ind++];
  char *qseq = argv[opt.ind++];

  int i;
  int8_t a = sc_mch, b = sc_mis < 0 ? sc_mis : -sc_mis; // a>0 and b<0
  int8_t mat[25] = {a, b, b, b, 0, b, a, b, b, 0, b, b, a,
                    b, 0, b, b, b, a, 0, 0, 0, 0, 0, 0};
  int tl = strlen(tseq), ql = strlen(qseq);
  uint8_t *ts, *qs, c[256];
  ksw_extz_t ez;

  memset(&ez, 0, sizeof(ksw_extz_t));
  memset(c, 4, 256);
  c['A'] = c['a'] = 0;
  c['C'] = c['c'] = 1;
  c['G'] = c['g'] = 2;
  c['T'] = c['t'] = 3; // build the encoding table
  ts = (uint8_t *)malloc(tl);
  qs = (uint8_t *)malloc(ql);
  for (i = 0; i < tl; ++i)
    ts[i] = c[(uint8_t)tseq[i]]; // encode to 0/1/2/3
  for (i = 0; i < ql; ++i)
    qs[i] = c[(uint8_t)qseq[i]];
  // ksw_extz(0, ql, qs, tl, ts, 5, mat, gapo, gape, -1, -1, 0, &ez);

  ksw_extd2_sse(0, ql, qs, tl, ts, 5, mat, gapo, gape, gapo2, gape2, -1, -1, -1,
                0, &ez);
  // ksw_extz2_sse(0, ql, qs, tl, ts, 5, mat, gapo, gape, -1, -1, zdrop, flag,
  // &ez);

  for (i = 0; i < ez.n_cigar; ++i) // print CIGAR
    printf("%d%c", ez.cigar[i] >> 4, "MID"[ez.cigar[i] & 0xf]);
  putchar('\n');

  printf("%d\n", ez.score);
  free(ez.cigar);
  free(ts);
  free(qs);
  return 0;
}

int main(int argc, char *argv[]) {
  if (strcmp(argv[1], "sketch") == 0)
    return main_sketch(argc - 1, argv + 1);
  else if (strcmp(argv[1], "dump") == 0)
    return main_dump(argc - 1, argv + 1);
  else if (strcmp(argv[1], "kan") == 0)
    return main_kan(argc - 1, argv + 1);
  else if (strcmp(argv[1], "chreads") == 0)
    return main_chreads(argc - 1, argv + 1);
  else if (strcmp(argv[1], "map") == 0)
    return main_map(argc - 1, argv + 1);
  else if (strcmp(argv[1], "align") == 0)
    return align(argc - 1, argv + 1);
  else if (strcmp(argv[1], "search") == 0)
    return main_search(argc - 1, argv + 1);
  else if (strcmp(argv[1], "call") == 0)
    return main_call(argc - 1, argv + 1);
  return 0;
}
