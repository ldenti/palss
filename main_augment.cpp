#include <cstdlib>
#include <filesystem>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <stack>
#include <string>
#include <tuple>
#include <utility>

#include "misc.hpp"
#include "usage.hpp"

#include "bdsg/packed_graph.hpp"
#include "gbwtgraph/algorithms.h"
#include "gbwtgraph/gbz.h"
#include "handlegraph/types.hpp"
#include "sdsl/simple_sds.hpp"

std::vector<std::vector<handlegraph::handle_t>>
dp(bdsg::PackedGraph *pg, handlegraph::handle_t source,
   handlegraph::handle_t sink, std::string sequence, size_t maxp) {
  std::vector<std::vector<handlegraph::handle_t>> paths;

  // Define DP[i] = set of vertices v such that there exists a path from source
  // to v that exactly matches S[:i] (so i matches)
  size_t m = sequence.size();
  std::string_view sequence_view(sequence);

  std::vector<std::set<handlegraph::handle_t>> DP(m + 1);
  size_t starting_position = pg->get_length(source);
  DP[starting_position].insert(source);
  // std::cerr << pg->get_id(source) << " " << starting_position << std::endl;
  // std::cerr << sequence_view.substr(0, starting_position) << " "
  //           << pg->get_sequence(source) << std::endl;
  assert(sequence_view.substr(0, starting_position)
             .compare(pg->get_sequence(source)) == 0);

  // store predecessors for backtracking
  // predecessors: pred[(v, j)] -> list of (u, i) where i + len(label[v]) == j
  std::map<std::pair<handlegraph::handle_t, int>,
           std::vector<std::pair<handlegraph::handle_t, int>>>
      predecessors;

  for (size_t i = starting_position; i < m + 1; ++i) {
    for (const auto &h : DP[i]) {

      pg->follow_edges(h, false, [&](const handlegraph::handle_t &out) {
        size_t j = i + pg->get_length(out);
        if (j <= m && sequence_view.substr(i, j - i).compare(
                          pg->get_sequence(out)) == 0) {
          DP[j].insert(out);
          predecessors[std::make_pair(out, j)].push_back(std::make_pair(h, i));
          // if (i == starting_position) {
          //   std::cerr << pg->get_id(h) << " " << i << std::endl;
          //   std::cerr << pg->get_id(out) << " " << j << std::endl;
          //   std::cerr << pg->get_id(
          //                    predecessors[std::make_pair(out, j)][0].first)
          //             << std::endl;
          //   std::cerr << predecessors[std::make_pair(out, j)][0].second
          //             << std::endl;
          // }
        }
      });
    }
  }

  // for (size_t i = 0; i < m + 1; ++i) {
  //   std::cerr << "DP[" << i << "] (" << DP[i].size() << ") [ ";
  //   for (const auto &h : DP[i]) {
  //     std::cerr << pg->get_id(h) << " ";
  //   }
  //   std::cerr << "]" << std::endl;
  // }

  // we should have sink in DP[m]
  // CHECKME: but in some cases we also have other vertices, why?
  if (DP[m].find(sink) == DP[m].end())
    return paths;

  assert(DP[starting_position].size() == 1);

  for (auto &kv : predecessors) {
    auto &v = kv.second;
    std::sort(v.begin(), v.end(),
              [&pg](std::pair<handlegraph::handle_t, int> &a,
                    std::pair<handlegraph::handle_t, int> &b) {
                return pg->get_id(a.first) > pg->get_id(b.first);
              });
  }

  // for (const auto &pred : predecessors) {
  //   std::cerr << "PRED[(" << pg->get_id(pred.first.first) << ","
  //             << pred.first.second << ")] (" << pred.second.size() << ") [ ";
  //   for (const auto &x : pred.second) {
  //     std::cerr << "(" << pg->get_id(x.first) << "," << x.second << ") ";
  //   }
  //   std::cerr << "]" << std::endl;
  // }

  // XXX: backtracking can be done way better. Avoid all those copies
  std::stack<std::tuple<handlegraph::handle_t, size_t,
                        std::vector<handlegraph::handle_t>>>
      stack;
  stack.push(
      std::make_tuple(sink, m, std::vector<handlegraph::handle_t>({sink})));

  while (!stack.empty() && paths.size() < maxp) {
    auto curr = stack.top();
    stack.pop();

    handlegraph::handle_t vertex = std::get<0>(curr);
    size_t pos = std::get<1>(curr);
    std::vector<handlegraph::handle_t> path = std::get<2>(curr);
    // std::cerr << "STACK " << pg->get_id(vertex) << " " << pos << std::endl;
    if (pos == starting_position) {
      std::reverse(path.begin(), path.end());
      paths.push_back(std::move(path));
    } else {
      for (const auto &pred : predecessors[std::make_pair(vertex, pos)]) {
        std::vector<handlegraph::handle_t> tmp = path;
        tmp.push_back(pred.first);
        stack.push(std::make_tuple(pred.first, pred.second, std::move(tmp)));
      }
    }
  }

  return paths;
}

bool flag_vertex(const handlegraph::handle_t &handle,
                 const std::set<handlegraph::handle_t> &known_vertices,
                 const std::map<handlegraph::handle_t, std::set<std::string>>
                     &novel_vertices,
                 size_t min_supp) {
  if (known_vertices.find(handle) == known_vertices.end()) {
    // if vertex is novel
    if (novel_vertices.find(handle) == novel_vertices.end() ||
        novel_vertices.at(handle).size() < min_supp) {
      // if vertex is not real novel OR its support is low
      return true;
    }
  }

  return false;
}

bool flag_edge(
    const handlegraph::edge_t &edge,
    const std::set<handlegraph::edge_t> &known_edges,
    const std::map<handlegraph::edge_t, std::set<std::string>> &novel_edges,
    size_t min_supp) {
  if (known_edges.find(edge) == known_edges.end()) {
    if (novel_edges.find(edge) == novel_edges.end() ||
        novel_edges.at(edge).size() < min_supp) {
      return true;
    }
  }
  return false;
}

int main_augment(int argc, char *argv[]) {
  double rt;

  bool chunk = false;
  size_t min_supp = 2;
  std::string wd = "/tmp";
  size_t maxp = 1024;
  std::string retained_gaf_fn = "";
  int nthreads = 4;
  // size_t max_plen = 100000;

  int _c;
  while ((_c = getopt(argc, argv, "cn:w:g:zs:@:h")) != -1) {
    switch (_c) {
    case 'c':
      chunk = true;
      break;
    case 'n':
      maxp = std::stoi(optarg);
      break;
    case 'w':
      wd = optarg;
      break;
    case '@':
      nthreads = std::stoi(optarg);
      break;
    case 'g':
      retained_gaf_fn = optarg;
      break;
    case 's':
      min_supp = std::stoi(optarg);
      break;
    case 'h':
      fprintf(stderr, "%s", AUGMENT_USAGE_MESSAGE);
      return 0;
    default:
      fprintf(stderr, "%s", AUGMENT_USAGE_MESSAGE);
      return 1;
    }
  }
  if (argc - optind != 2) {
    fprintf(stderr, "%s", AUGMENT_USAGE_MESSAGE);
    return 1;
  }

  std::string gbz_fn = argv[optind++];
  std::string gaf_fn = argv[optind++];

  std::filesystem::create_directories(wd);

  std::vector<std::pair<std::string, std::string>> file_pairs;

  rt = realtime();
  if (chunk) {
    size_t n_chunks;
    std::vector<std::set<std::string>> segments;
    {
      gbwtgraph::GBZ gbz;
      double rt0 = realtime();
      sdsl::simple_sds::load_from(gbz, gbz_fn);
      fprintf(stderr, "[M::%s] Load graph in %.3f sec\n", __func__,
              realtime() - rt0);

      // Chunking required 7m:35s and 32.8GB for HPRCv2
      rt0 = realtime();
      std::pair<std::vector<gbwtgraph::GBZ>, std::vector<std::string>> chunks =
          chunk_graph(gbz, {});
      n_chunks = chunks.first.size();
      fprintf(stderr, "[M::%s] Computed %ld chunks in %.3f sec\n", __func__,
              n_chunks, realtime() - rt0);
      segments.resize(chunks.first.size());

      // clang-format off
      // XXX: hardcoded to 2 threads to keep RAM usage <90GB
      /**
	 Chr GB  Min
	 ===========
	 1   42  40
	 2   41  44
	 3   31  34
	 4   34  34
	 5   32  32
	 6   30  31
	 7   32  30
	 8   29  27
	 9   23  23
	 10  24  24
	 11  25  24
	 12  24  23
	 13  19  19
	 14  18  19
	 15  17  17
	 16  16  17
	 17  15  15
	 18  12  13
	 19  11  12
	 20  12  12
	 21  7.5 8
	 22  9.4 10
	 X   15  15
	 Y   0   0
	 M   0   0
       **/
      // clang-format on
#pragma omp parallel for num_threads(2) schedule(static, 1)
      for (size_t i = 0; i < n_chunks; ++i) {
        double rt0 = realtime();

        const gbwtgraph::GBZ &sub_gbz = chunks.first[i];
        std::string fn = std::to_string(i); // + "_" + chunks.second[i];
        std::string sub_gbz_fn = wd + "/" + fn + ".gbz";
        std::string sub_pg_fn = wd + "/" + fn + ".pg";

        sdsl::simple_sds::serialize_to(chunks.first[i], sub_gbz_fn);
        std::ostringstream cmd;
        cmd << "vg convert --packed-out " << sub_gbz_fn << " > " << sub_pg_fn;

        // XXX: do this better
        (void)std::system(cmd.str().c_str());

        sub_gbz.graph.for_each_handle([&](const handlegraph::handle_t &handle) {
          segments[i].insert(sub_gbz.graph.get_segment_name(handle));
        });

#pragma omp critical(printf_lock)
        {
          fprintf(stderr, "[M::%s::%d] Processed chunk %ld (%s) in %.3f sec\n",
                  __func__, omp_get_thread_num(), i, chunks.second[i].c_str(),
                  realtime() - rt0);
          fflush(stdout);
        }
      }
    }

    // splitting GAF
    {
      std::vector<std::ofstream> sub_gafs(n_chunks + 1);
      for (size_t c = 0; c < sub_gafs.size(); ++c) {
        sub_gafs[c].open(wd + "/" + std::to_string(c) + ".gaf");
      }
      std::ifstream file(gaf_fn);
      std::string line;
      while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string token;
        int t = 0;

        std::string vertex; // XXX: is first vertex enough?
                            // std::vector<std::string> vertices;
        while (std::getline(iss, token, '\t')) {
          if (t == 5) {
            // size_t last_p = 1;
            size_t p;
            for (p = 1 /* last_p */; p < token.size(); ++p) {
              if (token[p] == '<' || token[p] == '>') {
                // vertices.push_back(token.substr(last_p, p - last_p));
                // last_p = p + 1;
                break;
              }
            }
            vertex = token.substr(1 /*last_p*/, p - 1 /*last_p*/);
            // vertices.push_back(token.substr(last_p, token.size() -
            // last_p));

            size_t c;
            for (c = 0; c < segments.size(); ++c) {
              if (segments[c].find(vertex) != segments[c].end())
                break;
            }
            sub_gafs[c] << line << std::endl;
          }
          ++t;
        }
      }
      for (auto &o : sub_gafs) {
        o.close();
      }
    }

    fprintf(stderr, "[M::%s] Chunked %ld components in %.3f sec\n", __func__,
            n_chunks, realtime() - rt);

    // augment chunks
    rt = realtime();
    {
      for (size_t c = 0; c < n_chunks; ++c) {
        std::string original_pg_fn = wd + "/" + std::to_string(c) + ".pg";
        std::string augmented_pg_fn =
            wd + "/" + std::to_string(c) + ".augmented.pg";
        std::string gaf_fn = wd + "/" + std::to_string(c) + ".gaf";
        std::ostringstream cmd;
        cmd << "vg augment --include-paths --min-coverage 1 --gaf "
            << original_pg_fn << " " << gaf_fn << " > " << augmented_pg_fn;

        // XXX: do this better
        (void)std::system(cmd.str().c_str());
        file_pairs.push_back({augmented_pg_fn, gaf_fn});

        fprintf(stderr, "[M::%s] Augmented chunk %ld in %.3f sec\n", __func__,
                c, realtime() - rt);
      }
    }
    fprintf(stderr, "[M::%s] Augmented %ld chunks in %.3f sec\n", __func__,
            n_chunks, realtime() - rt);
  } else {
    std::string pg_fn = wd + "/graph.pg";
    std::string augmented_pg_fn = wd + "/graph.augmented.pg";

    std::ostringstream cmd1;
    cmd1 << "vg convert --packed-out " << gbz_fn << " > " << pg_fn;
    // XXX: do this better
    (void)std::system(cmd1.str().c_str());

    std::ostringstream cmd2;
    cmd2 << "vg augment --include-paths --min-coverage 1 --gaf " << pg_fn << " "
         << gaf_fn << " > " << augmented_pg_fn;
    // XXX: do this better
    (void)std::system(cmd2.str().c_str());

    file_pairs.push_back({augmented_pg_fn, gaf_fn});
    fprintf(stderr, "[M::%s] Converted ad augmented graph in %.3f sec\n",
            __func__, realtime() - rt);
  }

  // =================================================================
  rt = realtime();
#pragma omp parallel for num_threads(nthreads) schedule(static, 1)
  for (size_t c = 0; c < file_pairs.size(); ++c) {
    double rt0 = realtime();
    std::string pg_fn = file_pairs[c].first;
    std::string gaf_fn = file_pairs[c].second;

    bdsg::PackedGraph *pg = new bdsg::PackedGraph();
    pg->deserialize(pg_fn);

    // std::cerr << "Loaded packed graph" << std::endl;

    // size_t total_paths = pg->get_path_count();

    std::set<handlegraph::handle_t> known_vertices;
    std::set<handlegraph::edge_t> known_edges;
    int pi = 0;
    int np = 0;

    pg->for_each_path_handle([&](const bdsg::path_handle_t p) {
      std::string pname = pg->get_path_name(p);
      std::string_view sv = pname;
      // std::cerr << pi << "/" << total_paths << " paths\r";
      // XXX: improve this, quite inefficient (batch insert?)
      if (sv.substr(0, 5).compare("palss") != 0) {
        handlegraph::handle_t last_h;
        bool first = true;
        // pg->for_each_step_in_path(p, [&](const bdsg::step_handle_t s) {
        //   std::cerr << pg->get_id(pg->get_handle_of_step(s)) << " ";
        // });
        // std::cerr << std::endl;

        for (const handlegraph::handle_t &h : pg->scan_path(p)) {
          // std::cerr << pg->get_id(h) << " ";
          known_vertices.insert(h);
          if (!first)
            known_edges.insert(std::make_pair(last_h, h));
          first = false;
          last_h = h;
        }
        // std::cerr << std::endl;
      } else {
        ++np;
      }
      ++pi;
    });
    // std::cerr << pi << "/" << total_paths << " paths" << std::endl;
    // std::cerr << "Novel paths: " << np << std::endl;
    // std::cerr << "Known vertices: " << known_vertices.size() << std::endl;
    // std::cerr << "Known edges: " << known_edges.size() << std::endl;

    // for (const auto &v : known_vertices)
    //   std::cerr << "NOVEL " << pg->get_id(v) << std::endl;

    std::map<std::string, std::set<std::string>> cluster_support;
    {
      std::ifstream file(gaf_fn);
      std::string line;
      while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string token;
        int i = 0;
        std::string cname;
        while (std::getline(iss, token, '\t')) {
          if (i == 0) {
            cname = token;
          } else if (i == 15) {
            std::string rname;
            std::istringstream iss2(token);
            iss2.ignore(5);
            while (std::getline(iss2, rname, '|')) {
              cluster_support[cname].insert(rname);
            }
          }
          ++i;
        }
      }
    }

    std::map<std::string, std::vector<handlegraph::handle_t>> real_novel;
    pg->for_each_path_handle([&](const bdsg::path_handle_t path) {
      std::string pname = pg->get_path_name(path);
      // if (pname.compare("palss-37356.0") != 0)
      //   return;
      std::string_view sv = pname;
      if (sv.substr(0, 5).compare("palss") == 0) {
        std::string seq;
        int plen = 0;
        for (const handlegraph::handle_t &h : pg->scan_path(path)) {
          ++plen;
          seq += pg->get_sequence(h); // already rc
          // pg->get_is_reverse(h)
          // pg->get_id(h)
        }

        handlegraph::handle_t source =
            pg->get_handle_of_step(pg->path_begin(path));
        handlegraph::handle_t sink =
            pg->get_handle_of_step(pg->path_back(path));

        if (source == sink) {
          // XXX: why do we have single-vertex paths?
          assert(plen == 1);
          return;
        }

        uint32_t best_novel_path_count = -1;

        std::vector<handlegraph::handle_t> best_novel_path;

        std::vector<std::vector<handlegraph::handle_t>> reconstructed_paths =
            dp(pg, source, sink, seq, maxp);

        assert(!reconstructed_paths.empty());

        for (const std::vector<handlegraph::handle_t> &rpath :
             reconstructed_paths) {

          // TODO: check edges?

          std::string rseq;
          size_t n_novel_vertices = 0;
          size_t n_novel_edges = 0;

          n_novel_vertices +=
              known_vertices.find(rpath[0]) == known_vertices.end();
          rseq += pg->get_sequence(rpath[0]);
          // std::cerr << pg->get_id(rpath[0]);
          for (size_t h = 1; h < rpath.size(); ++h) {
            // std::cerr << " " << pg->get_id(rpath[h]);
            n_novel_vertices +=
                known_vertices.find(rpath[h]) == known_vertices.end();
            // std::cerr << pg->get_id(rpath[h - 1]) << " " <<
            // pg->get_id(rpath[h])
            //           << " "
            //           << (known_edges.find(std::make_pair(
            //                   rpath[h - 1], rpath[h])) ==
            //                   known_edges.end())
            //           << std::endl;
            n_novel_edges +=
                known_edges.find(std::make_pair(rpath[h - 1], rpath[h])) ==
                known_edges.end();
            rseq += pg->get_sequence(rpath[h]);
          }
          // std::cerr << std::endl;
          // std::cerr << seq << std::endl;
          // std::cerr << rseq << std::endl;
          assert(seq.compare(rseq) == 0);

          // std::cerr << n_novel_vertices << " " << n_novel_edges <<
          // std::endl;

          if (n_novel_vertices == 0 && n_novel_edges == 0) {
            best_novel_path.clear();
            break;
          }
          // XXX: parsimony based on number of novel vertices, not edges
          if (n_novel_vertices < best_novel_path_count) {
            best_novel_path_count = n_novel_vertices;
            best_novel_path = rpath;
          }
        }
        if (!best_novel_path.empty()) {
          real_novel[pname] = best_novel_path;
          // std::cerr << pname << " " << best_novel_path_count << " "
          //           << best_novel_path.size() << std::endl;
        }
      }
    });

    // std::cerr << real_novel.size() << " novel paths will be retained"
    //           << std::endl;
    // std::cerr << "Selecting segments/links to remove..." << std::endl;

    // 2-pass cleaning. If a vertex is supported by at least one path, we want
    // to keep it even if it's not supported by others
    std::map<handlegraph::handle_t, std::set<std::string>> novel_vertices;
    std::map<handlegraph::edge_t, std::set<std::string>> novel_edges;

    // get novel vertices and edges from real novel paths
    for (const auto &pp : real_novel) {
      std::string cname = pp.first;
      // std::cerr << cname << " >";
      // for (const auto &h : pp.second) {
      //   std::cerr << " " << pg->get_id(h);
      // }
      // std::cerr << std::endl;

      if (known_vertices.find(pp.second[0]) == known_vertices.end())
        novel_vertices[pp.second[0]].insert(cluster_support[cname].begin(),
                                            cluster_support[cname].end());
      // std::cerr << pg->get_id(pp.second[0]) << " "
      //           << (known_vertices.find(pp.second[0]) ==
      //           known_vertices.end())
      //           << std::endl;
      for (size_t h = 1; h < pp.second.size(); ++h) {
        // std::cerr << pg->get_id(pp.second[h]) << " "
        //           << (known_vertices.find(pp.second[h]) ==
        //           known_vertices.end())
        //           << std::endl;
        if (known_vertices.find(pp.second[h]) == known_vertices.end())
          novel_vertices[pp.second[h]].insert(cluster_support[cname].begin(),
                                              cluster_support[cname].end());
        if (known_edges.find(std::make_pair(pp.second[h - 1], pp.second[h])) ==
            known_edges.end()) {
          novel_edges[std::make_pair(pp.second[h - 1], pp.second[h])].insert(
              cluster_support[cname].begin(), cluster_support[cname].end());
        }
      }
    }

    // std::cerr << novel_vertices.size() << " " << novel_edges.size() <<
    // std::endl;

    // iterate over original novel paths to get novel vertices not used by
    // real novel (cleaned) paths
    std::set<handlegraph::handle_t> vertices_to_remove;
    std::set<handlegraph::edge_t> edges_to_remove;

    pg->for_each_path_handle([&](const bdsg::path_handle_t p) {
      std::string pname = pg->get_path_name(p);
      std::string_view sv = pname;
      if (sv.substr(0, 5).compare("palss") == 0) {
        std::vector<handlegraph::handle_t> path;
        // if (real_novel.find(pname) == real_novel.end()) {
        for (const handlegraph::handle_t &h : pg->scan_path(p)) {
          // std::cerr << pg->get_id(h) << " ";
          path.push_back(h);
        }
        // } else {
        //   path = real_novel[pname];
        // }

        handlegraph::handle_t prev_handle;
        handlegraph::handle_t handle = path[0];
        if (flag_vertex(handle, known_vertices, novel_vertices, min_supp) &&
            flag_vertex(pg->flip(handle), known_vertices, novel_vertices,
                        min_supp)) {
          vertices_to_remove.insert(handle);
          vertices_to_remove.insert(pg->flip(handle));
        }

        for (size_t h = 1; h < path.size(); ++h) {
          prev_handle = path[h - 1];
          handle = path[h];

          // vertex
          if (flag_vertex(handle, known_vertices, novel_vertices, min_supp) &&
              flag_vertex(pg->flip(handle), known_vertices, novel_vertices,
                          min_supp)) {
            vertices_to_remove.insert(handle);
            vertices_to_remove.insert(pg->flip(handle));
          }

          // edge
          handlegraph::edge_t edge = std::make_pair(prev_handle, handle);
          handlegraph::edge_t edge2 =
              std::make_pair(pg->flip(handle), pg->flip(prev_handle));
          if (flag_edge(edge, known_edges, novel_edges, min_supp) &&
              flag_edge(edge2, known_edges, novel_edges, min_supp)) {
            edges_to_remove.insert(edge);
            edges_to_remove.insert(edge2);
          }
        }
      }
    });

    std::ofstream outfile(wd + "/" + std::to_string(c) + ".retained.gaf");
    std::ifstream file(gaf_fn);
    std::string line;
    while (std::getline(file, line)) {
      std::istringstream iss(line);
      std::string cname;
      std::getline(iss, cname, '\t');
      if (real_novel.find(cname) != real_novel.end())
        outfile << line << std::endl;
    }
    outfile.close();

    outfile.open(wd + "/" + std::to_string(c) + ".final.gfa");

    // === S LINES
    // =====================================================
    pg->for_each_handle(
        [&](const handlegraph::handle_t &handle) {
          // XXX: assuming here handle is always on + strand
          if (vertices_to_remove.find(handle) == vertices_to_remove.end())
            outfile << "S"
                    << "\t" << pg->get_id(handle) << "\t"
                    << pg->get_sequence(handle) << std::endl;
        },
        false);

    // === L LINES
    // =====================================================
    pg->for_each_edge(
        [&](const handlegraph::edge_t &edge) {
          if (edges_to_remove.find(edge) == edges_to_remove.end())
            outfile << "L\t" << pg->get_id(edge.first) << "\t"
                    << (pg->get_is_reverse(edge.first) ? '-' : '+') << "\t"
                    << pg->get_id(edge.second) << "\t"
                    << (pg->get_is_reverse(edge.second) ? '-' : '+') << "\t"
                    << "0M" << std::endl;
        },
        false);

    // === W LINES
    // =====================================================
    pg->for_each_path_handle([&](const bdsg::path_handle_t &path) {
      std::string pname = pg->get_path_name(path);
      std::string_view sv = pname;
      if (sv.substr(0, 5).compare("palss") != 0) {

        // XXX: this is copied from
        // https://github.com/vgteam/vg/blob/cee90d25878c71cdfb733d801e55548db4828a65/src/gfa.cpp#L190
        // FIXME: however, it seems to not work... With P-lines I get only 0,0
        // and then vg mod complains
        size_t start_offset = 0;
        // size_t end_offset = 0;
        auto subrange = pg->get_subrange(path);
        bool w_line = false;
        if (subrange != handlegraph::PathMetadata::NO_SUBRANGE) {
          start_offset = subrange.first;
          if (subrange.second != handlegraph::PathMetadata::NO_END_POSITION) {
            // end_offset = subrange.second;
            w_line = true;
          }
        }
        size_t path_length = 0;
        pg->for_each_step_in_path(
            path, [&](handlegraph::step_handle_t step_handle) {
              path_length +=
                  pg->get_length(pg->get_handle_of_step(step_handle));
            });
        // if (end_offset != 0 && start_offset + path_length != end_offset) {
        //   std::cerr << "[gfa] warning: incorrect end offset (" << end_offset
        //             << ") extracted from from path name "
        //             << pg->get_path_name(path) << ", using "
        //             << (start_offset + path_length) << " instead" <<
        //             std::endl;
        // }
        // =============================================================

        if (w_line) {
          outfile << "W"
                  << "\t" << pg->get_sample_name(path) << "\t"
                  << pg->get_haplotype(path) << "\t" << pg->get_locus_name(path)
                  << "\t" << start_offset << "\t" << start_offset + path_length
                  << "\t";
          for (const handlegraph::handle_t &handle : pg->scan_path(path)) {
            outfile << (pg->get_is_reverse(handle) ? '<' : '>')
                    << pg->get_id(handle);
          }
          outfile << std::endl;
        } else {
          outfile << "P"
                  << "\t" << pname << "\t";

          std::stringstream sp;
          for (const handlegraph::handle_t &handle : pg->scan_path(path)) {
            sp << pg->get_id(handle) << (pg->get_is_reverse(handle) ? "-" : "+")
               << ",";
          }
          std::string p = sp.str();
          p.pop_back();
          outfile << p << "\t"
                  << "*" << std::endl;
        }
      }
    });
    outfile.close();
    fprintf(stderr, "[M::%s] Refined chunk %ld in %.3f sec\n", __func__, c,
            realtime() - rt0);
  }
  fprintf(stderr, "[M::%s] Refined in %.3f sec\n", __func__, realtime() - rt);

  rt = realtime();
  std::ofstream out_gaf;
  if (!retained_gaf_fn.empty())
    out_gaf.open(retained_gaf_fn, std::ios::out);
  std::cout << "H\tVN:Z:1.1" << std::endl;
  for (size_t c = 0; c < file_pairs.size(); ++c) {
    std::ifstream in;
    in.open(wd + "/" + std::to_string(c) + ".final.gfa");
    std::cout << in.rdbuf();
    in.close();
    if (!retained_gaf_fn.empty()) {
      in.open(wd + "/" + std::to_string(c) + ".retained.gaf");
      out_gaf << in.rdbuf();
      in.close();
    }
    //
  }
  out_gaf.close();

  fprintf(stderr, "[M::%s] Merged in %.3f sec\n", __func__, realtime() - rt);

  return 0;
}
