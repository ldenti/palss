#include <cstdint>
#include <cstring>
#include <map>
#include <vector>
#include <zlib.h>

#include "graph.hpp"
extern "C" {
#include "sfs.h"
#include "sketch.h"
}

using namespace std;

typedef struct {
  map<string, int> references;
  map<int, string> v2ref;
  map<int, int> offsets;
} positions_t;

// Get positions of vertices in a given path with name pidx
positions_t get_positions(const Graph &graph, const std::string &pidx) {
  positions_t output;
  uint npaths = graph.gbz.index.metadata.paths();
  for (uint pp = 0; pp < npaths; ++pp) {
    int p = gbwt::Path::encode(pp, 0);
    // here we need the path identifier without strand
    std::string sample_name = graph.gbz.index.metadata.fullPath(pp).sample_name;
    std::string contig_name = graph.gbz.index.metadata.fullPath(pp).contig_name;
    // size_t haplotype = graph.gbz.index.metadata.fullPath(pp).haplotype;
    // std::cerr << sample_name << " " << contig_name << " " << haplotype
    //           << std::endl;
    int cp = 0;
    // here we need the encoded path identifier (with strand)
    if (sample_name.compare("_gbwt_ref") == 0) {
      gbwt::vector_type path = graph.gbz.index.extract(p);
      // std::cerr << path.size() << std::endl;
      for (const gbwt::node_type v : path) {
        // v has strand information
        gbwtgraph::nid_t vv = gbwt::Node::id(v);
        // vv is internal node_id
        gbwtgraph::handle_t hh = graph.gbz.graph.get_handle(vv);
        // string name = graph.gbz.graph.get_segment_name(hh);
        int l = graph.gbz.graph.get_length(hh);
        output.offsets[vv] = cp;
        output.v2ref[vv] = contig_name;
        cp += l;
      }
      output.references[contig_name] = cp;
    }
  }
  return output;
}

string d2s(const uint64_t kmer, int k) {
  string seq(k, 'N');
  for (int i = 1; i <= k; ++i)
    seq[i - 1] = "ACGT"[(kmer >> (k - i) * 2) & 3];
  return seq;
}

// XXX: this works for only one path at the time
int main_map(int argc, char *argv[]) {
  double rt = realtime();

  int klen = 27; // kmer size
  std::string pidx = "";

  int _c;
  while ((_c = getopt(argc, argv, "k:p:")) != -1) {
    switch (_c) {
    case 'k':
      klen = std::stoi(optarg);
      break;
    case 'p':
      pidx = optarg;
      break;
      // case 'h':
      //   fprintf(stderr, "%s", SEARCH_USAGE_MESSAGE);
      //   return 0;
    default:
      fprintf(stderr, "Error\n");
      return 1;
    }
  }
  if (argc - optind != 2) {
    fprintf(stderr, "Error!\n");
    return 1;
  }
  std::string gbz_fn = argv[optind++];
  std::string sfs_fn = argv[optind++];

  // Graph sketching and extraction
  sketch_t *sketch =
      sk_load((gbz_fn + ".k" + std::to_string(klen) + ".skt").c_str());
  fprintf(stderr, "[M::%s] loaded %ld sketches in %.3f sec\n", __func__,
          sketch->n, realtime() - rt);
  rt = realtime();

  // Graph
  Graph graph(gbz_fn);
  graph.load();
  // graph.print_stats();

  positions_t positions = get_positions(graph, pidx);
  fprintf(stderr, "[M::%s] extracted %ld paths in %.3f sec\n", __func__,
          positions.references.size(), realtime() - rt);
  rt = realtime();

  // ---

  // Iterate over specific strings to extract alignments
  printf("@HD\tVN:1.6\tSO:coordinate\n");
  for (auto const &[key, val] : positions.references) {
    printf("@SQ\tSN:%s\tLN:%d\n", key.c_str(), val);
  }

  FILE *fp;
  char *line = NULL;
  size_t len = 0;
  ssize_t read;

  fp = fopen(sfs_fn.c_str(), "r");
  if (fp == NULL)
    exit(EXIT_FAILURE);

  int v1, v2;
  string ref1, ref2;
  int pos1, pos2, d;
  while ((read = getline(&line, &len, fp)) != -1) {
    if (line[0] != 'O')
      continue;
    sfs_t ss = parse_sfs_line(line + 2);
    if (ss.a.v > ss.b.v) {
      fprintf(stderr, "%s %d %d - %d - %ld %d %ld %d\n", ss.rname, ss.s, ss.l,
              ss.strand, ss.a.v, ss.a.offset, ss.b.v, ss.b.offset);
      break;
    }
    assert(ss.a.v != -1 && ss.b.v != -1);
    v1 = ss.a.v;
    v2 = ss.b.v;

    ref1 = positions.v2ref[v1];
    ref2 = positions.v2ref[v2];

    assert(ref1.compare(ref2) == 0);

    if (positions.offsets.find(v1) == positions.offsets.end() ||
        positions.offsets.find(v2) == positions.offsets.end()) {
      // TODO: find best path by intersecting paths passing trough the two
      // vertices and report over it
      printf("%s\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n", ss.rname);
    } else {
      pos1 = positions.offsets.at(v1) + ss.a.offset;
      pos2 = positions.offsets.at(v2) + ss.b.offset;
      d = pos2 - (pos1 + klen);
      if (d < 0)
        continue;
      printf("%s\t0\t%s\t%d\t60\t%dM%dN%dM\t*\t0\t0\t%s%s\t*\n", ss.rname,
             ref1.c_str(), pos1 + 1, klen, d, klen, d2s(ss.a.seq, klen).c_str(),
             d2s(ss.b.seq, klen).c_str());
    }
  }
  fclose(fp);

  // ---

  return 0;
}
