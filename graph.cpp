#include "graph.hpp"

KSTREAM_INIT(gzFile, gzread, 65536)

Graph::Graph() {}
Graph::Graph(const std::string &_fn) { fn = _fn; }

int Graph::load_vertices() {
  kstring_t s = {0, 0, 0};
  int dret;
  gzFile fp = gzopen(fn.c_str(), "r");
  if (fp == 0)
    return 0;
  kstream_t *ks = ks_init(fp);
  seg_t seg;
  while (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) >= 0) {
    if (s.s[0] == 'S') {
      gfa_parse_S(s.s, seg);
      v_map[seg.idx] = vertices.size();
      vertices.push_back(seg);
      // sgms_add(vertices, seg.idx, seg.seq.c_str(), seg.l);
    }
  }
  // // TODO: make this something like "finalize()"
  // sgms_add(vertices, -1, "\0", 1);
  // --vertices->n;

  free(s.s);
  ks_destroy(ks);
  gzclose(fp);

  return 0;
}

bool Graph::is_on_ref(int v) {
  v = get_iidx(v);
  return vertices[v].path >= 0;
}

int Graph::load_edges() {
  if (vertices.empty())
    load_vertices();

  ne = 0;
  out_edges.resize(vertices.size());
  in_edges.resize(vertices.size());

  // load edges
  kstring_t s = {0, 0, 0};
  int dret;
  gzFile fp = gzopen(fn.c_str(), "r");
  if (fp == 0)
    return 0;
  kstream_t *ks = ks_init(fp);
  int v1 = 0, v2 = 0;
  while (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) >= 0) {
    if (s.s[0] == 'L') {
      gfa_parse_L(s.s, &v1, &v2);
      v1 = get_iidx(v1);
      v2 = get_iidx(v2);
      out_edges[v1].push_back(v2);
      in_edges[v2].push_back(v1);
      ++ne;
    }
  }

  free(s.s);
  ks_destroy(ks);
  gzclose(fp);

  return 0;
}

int Graph::load_paths(const std::string &ref) {
  if (vertices.empty())
    load_vertices();

  kstring_t s = {0, 0, 0};
  int dret;
  gzFile fp = gzopen(fn.c_str(), "r");
  if (fp == 0) {
    exit(1); // return paths;
  }
  kstream_t *ks = ks_init(fp);

  while (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) >= 0) {
    if (s.s[0] == 'P' || s.s[0] == 'W') {
      path_t path;
      if (s.s[0] == 'P') {
        gfa_parse_P(s.s, path, ref);
      } else {
        gfa_parse_W(s.s, path, ref);
      }

      if (path.vertices.empty())
        continue;

      // path must use "internal" idx
      for (uint v = 0; v < path.vertices.size(); ++v) {
        // printf("%d %d\n", path.vertices[v], get_iidx(path.vertices[v]));
        path.vertices[v] = get_iidx(path.vertices[v]);
        vertices[path.vertices[v]].path = paths.size();
        vertices[path.vertices[v]].pos = v;
      }
      paths.push_back(path);
    }
  }

  free(s.s);
  ks_destroy(ks);
  gzclose(fp);
  return 0;
}

uint64_t Graph::bfs(int v, int out) const {
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
      if (vertices[v].path != -1) {
        // we reached a reference vertex
        rv = v;
        break;
      }
      // add new vertices to next wave
      if (out) {
        for (const auto v1 : out_edges[v]) {
          if (visited.find(v1) == visited.end())
            next_wave.insert(v1);
        }
      } else {
        for (const auto v1 : in_edges[v]) {
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

// we do not need to compute the distance of all vertices, we will do it online
// and memoize the result
/**
int Graph::build_distance_index() {
  out_distances.resize(vertices.size());
  in_distances.resize(vertices.size());

  int index_size = 0;
  // TODO: multithreading
  for (uint v = 0; v < vertices.size(); ++v) {
    // check if vertex is on reference path
    if (vertices[v].path == -1) {
      // fprintf(stderr, "visiting vertex %d (%d)\n", v, vertices[v].idx);
      ++index_size;
      uint64_t res_out = bfs(v, 1);
      out_distances[v] = res_out;
      uint64_t res_in = bfs(v, 0);
      in_distances[v] = res_in;

      // printf("> %d (%d): distance %d from %d (%d)\n", v, vertices[v].idx,
      //        (uint32_t)res_out, (uint32_t)(res_out >> 32),
      //        vertices[res_out >> 32].idx);

      // printf("< %d (%d): distance %d from %d (%d)\n", v, vertices[v].idx,
      //        (uint32_t)res_in, (uint32_t)(res_in >> 32),
      //        vertices[res_in >> 32].idx);
    }
  }
  return 0;
}
**/

// bool Graph::cyclic_helper(int v, std::vector<bool> &recStack,
//                           std::vector<bool> &visited) const {
//   // If the node is already in the recursion stack, a cycle is detected
//   if (recStack[v])
//     return true;

//   // If the node is already visited and not in recursion stack, no need to
//   check
//   // again
//   if (visited[v])
//     return false;

//   // Mark the current node as visited and add it to the recursion stack
//   visited[v] = true;
//   recStack[v] = true;

//   // Recur for all neighbors
//   for (int x : out_edges.at(v)) {
//     if (cyclic_helper(x, visited, recStack))
//       return true;
//   }

//   // Remove the node from the recursion stack
//   recStack[v] = false;
//   return false;
// }

// bool Graph::cyclic() const {
//   int n = vertices.size();
//   std::vector<bool> visited(n, false);
//   std::vector<bool> recStack(n, false);
//   // Check for cycles starting from every unvisited node
//   for (int i = 0; i < n; ++i) {
//     if (!visited[i] && cyclic_helper(i, visited, recStack))
//       return true;
//   }
//   return false;
// }

uint Graph::distance(int v1, int v2) {
  // v1 and v2 are in graph space
  int ref1 = v1, d1 = 0;
  if (vertices[v1].path == -1) {
    // need to compute the distance to reference path (if not already done)
    if (out_distances.find(v1) == out_distances.end()) {
      out_distances[v1] = bfs(v1, 1);
      in_distances[v1] = bfs(v1, 0);
    }

    uint64_t res = out_distances[v1];
    d1 = (uint32_t)res;
    ref1 = (uint32_t)(res >> 32);
  }
  if (d1 == -1)
    return -1;

  int ref2 = v2, d2 = 0;
  if (vertices[v2].path == -1) {
    // need to compute the distance to reference path (if not already done)
    if (out_distances.find(v2) == out_distances.end()) {
      out_distances[v2] = bfs(v2, 1);
      in_distances[v2] = bfs(v2, 0);
    }

    uint64_t res = in_distances[v2];
    d2 = (uint32_t)res;
    ref2 = (uint32_t)(res >> 32);
  }

  if (d2 == -1)
    return -1;

  int path1 = vertices[ref1].path;
  int path2 = vertices[ref2].path;

  if (path1 != path2) {
    // std::cerr << v1 << " and " << v2 << " are on different path (" << path1
    //           << ", " << path2 << ")" << std::endl;
    return -1;
  }
  // === distance between vertices ref1<->ref2 on reference path path1 +
  // additional steps if original vertices are not on reference path
  return std::abs(vertices[ref1].pos - vertices[ref2].pos) + d1 + d2;
}

// bool Graph::dfs(int v1, int v2, int d,
//                 std::map<int, std::set<int>> &edges) const {
//   if (edges.find(v1) != edges.end())
//     // we already visited this vertex, that is a good vertex iif we have
//     stored
//     // at least one outgoing edge
//     return edges[v1].size() > 0;
//   if (v1 == v2)
//     // we have reached the sink
//     return true;
//   if (d > 20) // XXX: hardcoded
//               // we do not want to visit the entire graph
//     return false;

//   bool res = false;
//   edges[v1] = std::set<int>();
//   for (const auto &vv : out_edges[v1]) {
//     bool res1 = dfs(vv, v2, d + 1, edges);
//     if (res1)
//       edges[v1].insert(vv);
//     res |= res1;
//   }

//   return res;
// }

std::string Graph::get_path_sequence(int v1, int v2, int off1, int off2) const {
  // v1 and v2 are in GFA space
  v1 = get_iidx(v1);
  v2 = get_iidx(v2);
  int p1 = vertices[v1].path;
  int pos1 = vertices[v1].pos;
  int p2 = vertices[v2].path;
  int pos2 = vertices[v2].pos;

  if (p1 != p2) {
    fprintf(stderr, "vertices on different paths (%d>%d)", p1, p2);
    return "";
  }
  if (pos1 > pos2) {
    fprintf(stderr, "cycle detected (%d>%d)\n", pos1, pos2);
    return "";
  }
  std::string pseq = "";
  path_t path = paths[p1];

  for (int p = pos1; p <= pos2; ++p) {
    int v = path.vertices[p];
    std::string s;
    std::string seq = vertices[v].seq;
    if (v == v1 && v == v2) {
      s = seq.substr(off1, off2 - off1);
    } else {
      if (v == v1) {
        s = seq.substr(off1, seq.size());
      } else if (v == v2) {
        s = seq.substr(0, off2);
      }
    }
    pseq += s;
  }

  return pseq;
}

// Graph Graph::subgraph(int v1, int v2, int off1, int off2) const {
//   // v1 and v2 are in GFA space
//   v1 = get_iidx(v1);
//   v2 = get_iidx(v2);

//   Graph subgraph;
//   std::map<int, std::set<int>> edges;
//   if (v1 == v2)
//     edges[v1] = std::set<int>();
//   else {
//     dfs(v1, v2, 0, edges);

//     for (auto it = edges.cbegin(); it != edges.cend();) {
//       if (it->second.empty()) {
//         it = edges.erase(it);
//       } else {
//         ++it;
//       }
//     }
//   }

//   // First, add vertices to have size for edges
//   std::queue<int> qv;
//   qv.push(v1);
//   seg_t seg;
//   while (!qv.empty()) {
//     int v = qv.front();
//     qv.pop();
//     seg = vertices.at(v);
//     if (subgraph.v_map.find(seg.idx) == subgraph.v_map.end()) {
//       if (v == v1 && v == v2) {
//         seg.seq = seg.seq.substr(off1, off2 - off1);
//       } else {
//         if (v == v1) {
//           seg.seq = seg.seq.substr(off1, seg.seq.size());
//         } else if (v == v2) {
//           seg.seq = seg.seq.substr(0, off2);
//         }
//       }
//       seg.l = seg.seq.size();
//       subgraph.v_map[seg.idx] = subgraph.vertices.size();
//       subgraph.vertices.push_back(seg);
//     }
//     if (v == v2)
//       break;
//     if (edges.find(v) == edges.end())
//       continue;
//     for (const auto &vv : edges.at(v)) {
//       qv.push(vv);
//     }
//   }
//   assert(subgraph.v_map[vertices[v1].idx] == 0);
//   assert(subgraph.v_map[vertices[v2].idx] == subgraph.vertices.size() - 1);
//   subgraph.out_edges.resize(subgraph.vertices.size());
//   subgraph.in_edges.resize(subgraph.vertices.size());
//   subgraph.ne = 0;

//   // Now, reiterate to add edges
//   int i1, i2;
//   for (auto const &[v, outs] : edges) {
//     seg = vertices.at(v);
//     i1 = subgraph.get_iidx(seg.idx);
//     for (const auto &v1 : outs) {
//       seg = vertices.at(v1);
//       i2 = subgraph.get_iidx(seg.idx);
//       subgraph.out_edges[i1].push_back(i2);
//       subgraph.in_edges[i2].push_back(i1);
//       ++subgraph.ne;
//     }
//   }

//   return subgraph;
// }

// Graph Graph::canonicalize() const {
//   int total_size = 0;
//   for (uint v = 0; v < vertices.size(); ++v)
//     total_size += vertices[v].l;

//   Graph cg;

//   // vertices
//   // cg.vertices.resize(total_size);

//   // edges
//   cg.ne = 0;
//   cg.out_edges.resize(total_size);
//   cg.in_edges.resize(total_size);

//   // TODO: add paths if we need

//   // for each vertex of original graph (internal id in [0,n)), store
//   // starting/ending canonical vertices in new graph
//   std::vector<int> starting(vertices.size(), 0);
//   std::vector<int> ending(vertices.size(), 0);

//   int v = 0;
//   for (const seg_t &seg : vertices) {
//     starting[v] = cg.vertices.size();
//     for (int c = 0; c < seg.l; ++c) {
//       seg_t new_seg;
//       new_seg.idx = cg.vertices.size();
//       new_seg.seq = seg.seq[c];
//       new_seg.l = 1;
//       cg.vertices.push_back(new_seg);
//       cg.v_map[new_seg.idx] = new_seg.idx;

//       // internal link
//       if (c > 0) {
//         cg.out_edges[new_seg.idx - 1].push_back(new_seg.idx);
//         ++cg.ne;
//         cg.in_edges[new_seg.idx].push_back(new_seg.idx - 1);
//       }
//     }
//     ending[v] = cg.vertices.size() - 1;
//     ++v;
//   }

//   // links between original vertices
//   int vv;
//   for (int v = 0; v < vertices.size(); ++v) {
//     vv = ending[v];
//     for (const int &xx : out_edges[v]) {
//       cg.out_edges[vv].push_back(starting[xx]);
//       ++cg.ne;
//       cg.in_edges[starting[xx]].push_back(vv);
//     }
//   }

//   return cg;
// }

// int Graph::to_gfa(const std::string &fn) const {
//   // FILE *fp = fn.compare("-") != 0 ? fopen(fn.c_str(), "w")
//   //                                 : fdopen(fileno(stdout), "w");
//   // if (fp == 0)
//   //   return -1;

//   // fprintf(fp, "H\tVN:Z:1.1\n");
//   printf("H\tVN:Z:1.1\n");
//   for (const seg_t &seg : vertices)
//     printf("S\t%d\t%s\n", seg.idx, seg.seq.c_str());
//   // fprintf(fp, "S\t%d\t%s\n", seg.idx, seg.seq.c_str());
//   for (uint v = 0; v < out_edges.size(); ++v) {
//     for (const auto &v1 : out_edges[v]) {
//       printf("L\t%d\t+\t%d\t+\t0M\n", vertices.at(v).idx,
//       vertices.at(v1).idx);
//       // fprintf(fp, "L\t%d\t+\t%d\t+\t0M\n", vertices.at(v).idx,
//       //         vertices.at(v1).idx);
//     }
//   }
//   // fclose(fp);
//   return 0;
// }

int Graph::get_iidx(int v) const { return v_map.at(v); }

std::string Graph::get_sequence(int v) const {
  // v is in GFA space
  return vertices.at(get_iidx(v)).seq;
}

/** === GFA PARSING ======= **/
void gfa_parse_S(char *s, seg_t &ret) {
  int i;       // , is_ok = 0;
  char *p, *q; // *seg = 0, *seq = 0, *rest = 0;
  for (i = 0, p = q = s + 2;; ++p) {
    if (*p == 0 || *p == '\t') {
      int c = *p;
      *p = 0;
      if (i == 0) {
        ret.idx = atoi(q);
        // strcpy(ret->idx, q);
      } else if (i == 1) {
        ret.l = p - q;
        ret.seq = q;
        // is_ok = 1, rest = c ? p + 1 : 0;
        break;
      }
      ++i, q = p + 1;
      if (c == 0)
        break;
    }
  }
  // if (!is_ok) { // something is missing
}
void gfa_parse_L(char *s, int *idx1, int *idx2) {
  int i;       // , oriv = -1, oriw = -1, is_ok = 0;
  char *p, *q; // , *segv = 0, *segw = 0, *rest = 0;
  // int32_t ov = INT32_MAX, ow = INT32_MAX;
  for (i = 0, p = q = s + 2;; ++p) {
    if (*p == 0 || *p == '\t') {
      int c = *p;
      *p = 0;
      if (i == 0) {
        *idx1 = atoi(q);
      } else if (i == 1) {
        // if (*q != '+' && *q != '-')
        //   ; // return -2;
        // oriv = (*q != '+');
      } else if (i == 2) {
        *idx2 = atoi(q);
      } else if (i == 3) {
        // if (*q != '+' && *q != '-')
        //   ; // return -2;
        // oriw = (*q != '+');
      }
      ++i, q = p + 1;
      if (c == 0)
        break;
    }
  }
  // return 0;
}

void gfa_parse_P(char *s, path_t &path, const std::string &ref) {
  int i;
  char *p, *q, *qq;
  for (i = 0, p = q = s + 2;; ++p) {
    if (*p == 0 || *p == '\t') {
      *p = 0;
      if (i == 0) {
        path.idx = q;
        if (path.idx.find(ref) == std::string::npos &&
            path.idx.compare(0, 3, "chr") != 0)
          return;
      } else if (i == 1) {
        // char strand = *(p - 1);
        qq = q;
        for (qq = q;; ++qq) {
          if (*qq == 0 || *qq == ',') {
            int c = *qq;
            *qq = 0;
            // strand = *(qq - 1);
            *(qq - 1) = 0;
            // v = get_iidx(g, atoi(q));
            path.vertices.push_back(atoi(q));
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

void gfa_parse_W(char *s, path_t &path, const std::string &ref) {
  int i;
  char *p, *q, *qq;
  for (i = 0, p = q = s + 2;; ++p) {
    if (*p == 0 || *p == '\t') {
      *p = 0;
      if (i < 2) {
        *p = '#';
        ++i;
        continue;
      } else if (i == 2) {
        path.idx = q;
        if (path.idx.compare(0, ref.size(), ref) != 0)
          return;
      } else if (i == 5) {
        // char strand = *q;
        ++q;
        for (qq = q;; ++qq) {
          if (*qq == 0 || *qq == '>' || *qq == '<') {
            int c = *qq;
            *qq = 0;
            // v = get_iidx(g, atoi(q));
            // ph_addv(path, v, strand == '>');
            path.vertices.push_back(atoi(q));
            q = qq + 1;
            if (c == 0)
              break;
            // strand = c;
          }
        }
        break;
      }
      ++i, q = p + 1;
    }
  }
}
