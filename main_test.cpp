#include <string>

#include "graph.hpp"
#include "poa.hpp"

int main_test(int argc, char *argv[]) {
  double rt = realtime();

  char *gfa_fn = argv[1];
  int v1 = std::stoi(argv[2]);
  int v2 = std::stoi(argv[3]);
  std::string seq = argv[4];

  Graph graph(gfa_fn);
  rt = realtime();
  graph.load_vertices();
  graph.load_edges();
  fprintf(stderr, "loaded %ld vertices and %d edges in %.3f secs\n",
          graph.vertices.size(), graph.ne, realtime() - rt);

  Graph subgraph = graph.subgraph(v1, v2, 1, 1);
  fprintf(stderr, "extracted %ld vertices and %d edges in %.3f secs\n",
          subgraph.vertices.size(), subgraph.ne, realtime() - rt);
  fprintf(stderr, "%s\n", subgraph.cyclic() ? "cyclic" : "acyclic");

  // subgraph.to_gfa("-");

  Graph cangraph = subgraph.canonicalize();
  fprintf(stderr, "extracted %ld vertices and %d edges in %.3f secs\n",
          cangraph.vertices.size(), cangraph.ne, realtime() - rt);

  fprintf(stderr, "%s\n", cangraph.cyclic() ? "cyclic" : "acyclic");

  //   cangraph.to_gfa("-");

  align(cangraph, seq);

  return 0;
}
