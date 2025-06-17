#include <iostream>
#include <set>
#include <stdint.h>
#include <stdio.h>

#include "graph.hpp"

bool dfs(const Graph &graph, int v1, int v2, int d,
         std::map<int, std::set<int>> &visited) {
  if (visited.find(v1) != visited.end())
    return true; // visited[v1].size() > 0;
  if (v1 == v2)
    return true;
  if (d > 20) // XXX: hardcoded
    return false;

  bool res = false;
  for (const auto &vv : graph.out_edges[v1]) {
    std::cerr << "Visiting " << v1 << ">" << vv << " : ";
    bool res1 = dfs(graph, vv, v2, d + 1, visited);
    std::cerr << res1 << std::endl;
    if (res1) {
      visited[v1].insert(vv);
    }
    res |= res1;
  }

  return res;
}

void subgraph(const Graph &graph, int v1, int v2) {
  std::map<int, std::set<int>> visited;

  dfs(graph, v1, v2, 0, visited);

  for (const auto &it : visited) {
    std::cout << it.first << ":";
    for (const auto &v : it.second)
      std::cout << " " << v;
    std::cout << std::endl;
  }

  for (const auto &it : visited) {
    std::cout << graph.vertices[it.first].idx << ":";
    for (const auto &v : it.second)
      std::cout << " " << graph.vertices[v].idx;
    std::cout << std::endl;
  }

  return;
}

uint64_t bfs(const Graph &graph, int v, int out) {
  std::vector<int> curr_wave;
  std::set<int> next_wave;
  std::set<int> visited;

  curr_wave.push_back(v);
  int distance = 0;
  int rv = -1;
  while (1) {
    while (!curr_wave.empty()) {
      v = curr_wave.back();
      visited.insert(v);
      curr_wave.pop_back();
      if (graph.vertices[v].path != -1) {
        // we reached a reference vertex
        rv = v;
        break;
      }
      // add new vertices to next wave
      if (out) {
        for (const auto v1 : graph.out_edges[v]) {
          if (visited.find(v1) == visited.end())
            next_wave.insert(v1);
        }
      } else {
        for (const auto v1 : graph.in_edges[v]) {
          if (visited.find(v1) == visited.end())
            next_wave.insert(v1);
        }
      }
    }
    if (rv != -1)
      // we have found a reference vertex
      break;
    // Prepare next wave
    for (const auto &v1 : next_wave)
      curr_wave.push_back(v1);
    next_wave.clear();
    if (curr_wave.empty()) {
      distance = -1;
      break;
    }
    ++distance;
  }

  return ((uint64_t)rv << 32) | distance;
}

int distance(const Graph &graph, int v1, int v2,
             std::vector<uint64_t> in_distances,
             std::vector<uint64_t> out_distances) {
  // assuming v1 precedes v2
  int ref1 = v1, d1 = 0;
  if (graph.vertices[v1].path == -1) {
    uint64_t res = out_distances[v1];
    d1 = (uint32_t)res;
    ref1 = (uint32_t)(res >> 32);
  }

  int ref2 = v2, d2 = 0;
  if (graph.vertices[v2].path == -1) {
    uint64_t res = in_distances[v2];
    d2 = (uint32_t)res;
    ref2 = (uint32_t)(res >> 32);
  }

  int path1 = graph.vertices[ref1].path;
  int path2 = graph.vertices[ref2].path;

  if (path1 != path2) {
    std::cerr << v1 << " and " << v2 << " are on different path (" << path1
              << ", " << path2 << ")" << std::endl;
    return -1;
  }
  // === distance between vertices ref1->ref2 on reference path path1
  int d = std::abs(graph.vertices[ref1].pos - graph.vertices[ref2].pos);
  // ===

  printf("d(%d,%d)=d(%d,%d)=%d + %d + %d\n", ref1, ref2,
         graph.vertices[ref1].idx, graph.vertices[ref2].idx, d, d1, d2);
  return d + d1 + d2;
}

int main_test(int argc, char *argv[]) {
  double rt = realtime();

  char *gfa_fn = argv[1];
  // int n = std::stoi(argv[2]);
  int v1 = std::stoi(argv[2]);
  int v2 = std::stoi(argv[3]);

  Graph graph(gfa_fn);
  rt = realtime();
  graph.load_vertices();
  fprintf(stderr, "loaded %ld vertices in %.3f secs\n", graph.vertices.size(),
          realtime() - rt);

  rt = realtime();
  graph.load_edges();
  fprintf(stderr, "loaded %d edges in %.3f secs\n", graph.ne, realtime() - rt);

  subgraph(graph, v1, v2);

  // rt = realtime();
  // graph.load_paths("CHM13");
  // fprintf(stderr, "loaded %ld paths in %.3f secs\n", graph.paths.size(),
  //         realtime() - rt);

  // rt = realtime();
  // std::vector<uint64_t> out_distances(graph.vertices.size());
  // std::vector<uint64_t> in_distances(graph.vertices.size());
  // int index_size = 0;
  // for (uint v = 0; v < graph.vertices.size(); ++v) {
  //   // check if vertex is on reference path
  //   if (graph.vertices[v].path == -1) {
  //     // fprintf(stderr, "visiting vertex %d (%d)\n", v,
  //     graph.vertices[v].idx);
  //     ++index_size;
  //     uint64_t res_out = bfs(graph, v, 1);
  //     out_distances[v] = res_out;
  //     uint64_t res_in = bfs(graph, v, 0);
  //     in_distances[v] = res_in;
  //     // printf("> %d (%d): distance %d from %d (%d)\n", v,
  //     // graph.vertices[v].idx,
  //     //        (uint32_t)res_out, (uint32_t)(res_out >> 32),
  //     //        graph.vertices[res_out >> 32].idx);

  //     // printf("< %d (%d): distance %d from %d (%d)\n", v,
  //     // graph.vertices[v].idx,
  //     //        (uint32_t)res_in, (uint32_t)(res_in >> 32),
  //     //        graph.vertices[res_in >> 32].idx);
  //   }
  // }
  // fprintf(stderr, "built distance index (%d vertices) in %.3f secs\n",
  //         index_size, realtime() - rt);

  // rt = realtime();
  // int i = 0;
  // while (i < n) {
  //   int v1 = rand() % graph.vertices.size();
  //   int v2 = rand() % graph.vertices.size();
  //   if (v1 > v2) {
  //     std::swap(v1, v2);
  //   }
  //   int d = distance(graph, v1, v2, in_distances, out_distances);
  //   printf("d(%d,%d)=d(%d,%d)=%d\n\n", v1, v2, graph.vertices[v1].idx,
  //          graph.vertices[v2].idx, d);
  //   ++i;
  // }
  // fprintf(stderr, "queried %d distances in %.3f secs\n", n, realtime() - rt);

  // // int v1 = 9;
  // // int v2 = 46;
  // // int d = distance(graph, v1, v2, in_distances, out_distances);

  return 0;
}
