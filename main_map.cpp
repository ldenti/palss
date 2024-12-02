#include <cstdint>
#include <cstring>
#include <map>
#include <zlib.h>

// #include "gsketch.hpp"
#include "ketopt.h"
#include "kseq.h"
#include "utils.h"

#include "sketch.hpp"
extern "C" {
#include "graph.h"
}

// KSEQ_INIT(gzFile, gzread) // we already init kstream in graph.h
// XXX: there should be a better way to do this
__KSEQ_TYPE(gzFile)
__KSEQ_BASIC(static, gzFile)
__KSEQ_READ(static)

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
      for (int i = 0; i < path->l; ++i) {
        positions[path->vertices[i]] = cp;
        cp += get_vertex(graph, path->vertices[i])->l;
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

  // Iterate over fasta to extract alignments
  gzFile fp = gzopen(sfs_fn, "r");
  kseq_t *seq = kseq_init(fp);
  int l;
  char k1[klen + 1], k2[klen + 1];
  k1[klen] = '\0';
  k2[klen] = '\0';

  uint64_t k1_d, k2_d, kmer_d, rckmer_d;
  pair<int64_t, int16_t> p1, p2;
  uint64_t v1, v2;
  int pos1, pos2, d;

  printf("@HD\tVN:1.6\tSO:coordinate\n");
  printf("@SQ\tSN:%s\tLN:%d\n", pidx, totl);
  while ((l = kseq_read(seq)) >= 0) {
    strncpy(k1, seq->seq.s, klen);
    kmer_d = k2d(k1, klen);
    rckmer_d = rc(kmer_d, klen);
    k1_d = std::min(kmer_d, rckmer_d);

    strncpy(k2, seq->seq.s + (l - klen), klen);
    kmer_d = k2d(k2, klen);
    rckmer_d = rc(kmer_d, klen);
    k2_d = std::min(kmer_d, rckmer_d);

    p1 = sk_get(sketch, k1_d);
    v1 = p1.first;
    p2 = sk_get(sketch, k2_d);
    v2 = p2.first;

    if (v1 == -1 || v2 == -1)
      // unanchored specific string
      continue;

    pos1 = positions.at(v1) + p1.second;
    pos2 = positions.at(v2) + p2.second;

    d = pos2 - (pos1 + klen);

    printf("%s\t0\t%s\t%d\t60\t%dM%dN%dM\t*\t0\t0\t%s%s\t*\n", seq->name.s,
           pidx, pos1 + 1, klen, d, klen, k1, k2);
  }
  kseq_destroy(seq);
  gzclose(fp);

  // ---

  destroy_graph(graph);
  fprintf(stderr, "[M::%s] completed in %.3f sec\n", __func__,
          realtime() - rt0);

  return 0;
}
