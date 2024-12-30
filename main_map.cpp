#include <cstdint>
#include <cstring>
#include <map>
#include <zlib.h>

#include "ketopt.h"
#include "sfs.h"
#include "utils.h"

#include "sketch.hpp"
extern "C" {
#include "graph.h"
}

using namespace std;

/* Get positions of vertices in a given path with name pidx */
map<int, int> get_positions(graph_t *graph, char *pidx) {
  map<int, int> positions;
  path_t *path;
  int cp = 0;
  int get_first = strcmp(pidx, "") == 0;
  for (int p = 0; p < graph->np; ++p) {
    path = graph->paths[p];
    if (get_first || strcmp(path->idx, pidx) == 0) {
      seg_t *s;
      for (int i = 0; i < path->l; ++i) {
        s = graph->vertices[path->vertices[i]];
        // fprintf(stderr, "%d : %d\n", path->vertices[i], s->idx);
        positions[s->idx] = cp;
        cp += s->l;
      }
      break;
    }
  }
  return positions;
}

int main_map(int argc, char *argv[]) {
  // XXX: this works for only one path at the time
  int klen = 27; // kmer size
  char pidx[128] = "";
  static ko_longopt_t longopts[] = {{NULL, 0, 0}};
  ketopt_t opt = KETOPT_INIT;
  int _c;
  while ((_c = ketopt(&opt, argc, argv, 1, "k:p:", longopts)) >= 0) {
    if (_c == 'k')
      klen = atoi(opt.arg);
    else if (_c == 'p')
      strncpy(pidx, opt.arg, strlen(opt.arg));
  }
  if (argc - opt.ind != 3) {
    fprintf(stderr, "Argh");
    return 1;
  }
  char *gfa_fn = argv[opt.ind++];
  char *skt_fn = argv[opt.ind++];
  char *sfs_fn = argv[opt.ind++];

  double rt0, rt;
  rt0 = realtime();
  rt = rt0;

  // Load sketch and graph
  sketch_t sketch;
  sk_load(sketch, skt_fn);
  fprintf(stderr, "[M::%s] loaded %ld sketches in %.3f sec\n", __func__,
          sketch.size(), realtime() - rt);
  rt = realtime();

  graph_t *graph = init_graph(gfa_fn);
  load_vertices(graph);
  fprintf(stderr, "[M::%s] loaded %d vertices in %.3f sec\n", __func__,
          graph->nv, realtime() - rt);
  rt = realtime();
  load_paths(graph);
  fprintf(stderr, "[M::%s] loaded %d paths in %.3f sec\n", __func__, graph->np,
          realtime() - rt);
  rt = realtime();

  map<int, int> positions = get_positions(graph, pidx);
  int totl = 0;
  for (const auto &p : positions)
    totl += get_vertex(graph, p.first)->l;
  fprintf(stderr, "[M::%s] extracted path '%s' (size: %d) in %.3f sec\n",
          __func__, pidx, totl, realtime() - rt);
  rt = realtime();

  // ---

  // Iterate over specific strings to extract alignments

  printf("@HD\tVN:1.6\tSO:coordinate\n");
  printf("@SQ\tSN:%s\tLN:%d\n", pidx, totl);

  FILE *fp;
  char *line = NULL;
  size_t len = 0;
  ssize_t read;

  fp = fopen(sfs_fn, "r");
  if (fp == NULL)
    exit(EXIT_FAILURE);

  sfs_t ss;
  pair<int64_t, int16_t> p1, p2;
  uint64_t v1, v2;
  int pos1, pos2, d;
  while ((read = getline(&line, &len, fp)) != -1) {
    if (line[0] == 'X')
      continue;
    ss = parse_sfs_line(line + 2);
    if (ss.a.v > ss.b.v) {
      fprintf(stderr, "%s %d %d - %d - %d %d %d %d\n", ss.rname, ss.s, ss.l,
              ss.strand, ss.a.v, ss.a.offset, ss.b.v, ss.b.offset);
      break;
    }
    p1 = sk_get(sketch, ss.a.seq);
    v1 = p1.first;
    p2 = sk_get(sketch, ss.b.seq);
    v2 = p2.first;

    if (v1 != -1 && v2 != -1) {
      // anchored specific string
      if (positions.find(v1) == positions.end() ||
          positions.find(v2) == positions.end()) {
        // TODO: find best path by intersecting paths passing trough the two
        // vertices and report over it
        printf("%s\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n", ss.rname);
      } else {
        pos1 = positions.at(v1) + p1.second;
        pos2 = positions.at(v2) + p2.second;
        d = pos2 - (pos1 + klen);

        printf("%s\t0\t%s\t%d\t60\t%dM%dN%dM\t*\t0\t0\t%s%s\t*\n", ss.rname,
               pidx, pos1 + 1, klen, d, klen, d2s(ss.a.seq, klen).c_str(),
               d2s(ss.b.seq, klen).c_str());
      }
    }
    free(ss.rname);
    free(ss.seq);
  }
  fclose(fp);

  // ---

  destroy_graph(graph);
  fprintf(stderr, "[M::%s] completed in %.3f sec\n", __func__,
          realtime() - rt0);

  return 0;
}
