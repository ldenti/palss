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

void rewrite_P_line(std::string &line,
                    std::map<std::string, std::string> &mapping) {
  size_t p = 0;
  size_t i = 0;
  while (p <= line.size()) {
    std::size_t pos = line.find('\t', p);
    if (pos == std::string::npos)
      pos = line.size();
    if (i == 0) {
      std::cout << "P\t";
    } else if (i == 1) {
      std::cout << line.substr(p, pos - p) << "\t";
    } else if (i == 2) {
      size_t idx = 0;
      while (p < pos) {
        char c = line[p];
        if (std::isdigit(static_cast<unsigned char>(c))) {
          idx = idx * 10 + (c - '0');
        } else if (c == '+' || c == '-') {
          std::cout << mapping[std::to_string(idx)] << c << line[p + 1];
          idx = 0;
          ++p; // skip comma
        }
        ++p;
      }
    } else {
      std::cout << "*" << std::endl;
    }

    ++i;
    p = pos + 1;

    if (pos == line.size())
      break;
  }
}

void rewrite_W_line(std::string &line,
                    std::map<std::string, std::string> &mapping) {
  size_t p = 0;
  size_t i = 0;
  while (p <= line.size()) {
    std::size_t pos = line.find('\t', p);
    if (pos == std::string::npos)
      pos = line.size();

    if (i < 6) {
      std::cout << line.substr(p, pos - p) << "\t";
    } else if (i == 6) {
      size_t idx = 0;
      char c = line[p];
      char c0 = c;
      ++p;
      while (p < pos) {
        c = line[p];
        if (std::isdigit(static_cast<unsigned char>(c))) {
          idx = idx * 10 + (c - '0');
        } else {
          std::cout << c0 << mapping[std::to_string(idx)];
          idx = 0;
          c0 = c;
        }
        ++p;
      }
      std::cout << c0 << mapping[std::to_string(idx)] << std::endl;
    }

    ++i;
    p = pos + 1;

    if (pos == line.size())
      break;
  }

  // std::vector<size_t> path;
  // size_t s = 0;
  // size_t i = 0;
  // while (s <= line.size()) {
  //   std::size_t pos = line.find('\t', s);
  //   if (pos == std::string::npos)
  //     pos = line.size();
  //   if (i == 1) {
  //     sample = line.substr(s, pos - s);
  //   } else if (i == 2) {
  //     haplotype = line.substr(s, pos - s);
  //   } else if (i == 3) {
  //     contig = line.substr(s, pos - s);
  //   } else if (i == 4) {
  //     start = line.substr(s, pos - s);
  //   } else if (i == 5) {
  //     end = line.substr(s, pos - s);
  //   } else {
  //     bool strand = line[s] == '>';
  //     ++s;
  //     size_t ss = 0;
  //     while (s <= line.size()) {
  //       if (line[s + ss] == '>' || line[s + ss] == '<') {
  //         std::string idx_s = line.substr(s, ss - s);
  //         std::cerr << strand << " " << idx_s << std::endl;
  //       }
  //     }

  //     // seq = line.substr(start, pos - start);
  //     break;
  //   }
  //   ++i;
  //   s = pos + 1;

  //   if (pos == line.size())
  //     break;
  // }
  // return path;
}

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
        }
      });
    }
  }

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

bool flag_vertex(
    const handlegraph::handle_t &handle,
    const std::set<handlegraph::handle_t> &known_vertices,
    const std::map<handlegraph::handle_t,
                   std::pair<std::set<std::string>, std::set<std::string>>>
        &novel_vertices,
    size_t min_supp) {
  if (known_vertices.find(handle) == known_vertices.end()) {
    // if vertex is novel
    if (novel_vertices.find(handle) == novel_vertices.end() ||
        novel_vertices.at(handle).second.size() < min_supp) {
      // if vertex is not real novel OR its support is low
      return true;
    }
  }

  return false;
}

bool flag_edge(
    const handlegraph::edge_t &edge,
    const std::set<handlegraph::edge_t> &known_edges,
    const std::map<handlegraph::edge_t,
                   std::pair<std::set<std::string>, std::set<std::string>>>
        &novel_edges,
    size_t min_supp) {
  if (known_edges.find(edge) == known_edges.end()) {
    if (novel_edges.find(edge) == novel_edges.end() ||
        novel_edges.at(edge).second.size() < min_supp) {
      return true;
    }
  }
  return false;
}

int main_augment(int argc, char *argv[]) {
  // === CLI =========================================================
  bool chunk = false;
  bool force = false;
  bool delete_temp = false;
  size_t min_supp = 2;
  std::string wd = "/tmp";
  size_t maxp = 1024;
  std::string retained_gaf_fn = "";
  int nthreads = 4;

  int _c;
  while ((_c = getopt(argc, argv, "cfdn:w:g:zs:@:h")) != -1) {
    switch (_c) {
    case 'c':
      chunk = true;
      break;
    case 'f':
      force = true;
      break;
    case 'd':
      delete_temp = true;
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
  double rt = realtime();

  // =================================================================

  std::vector<std::pair<std::string, std::string>> file_pairs;

  size_t n_chunks;
  // Use existing chunks if not -f
  std::string done_fp = wd + "/DONE";
  std::ifstream f(done_fp);
  if (f.good() && !force) {
    f >> n_chunks;
    for (size_t c = 0; c < n_chunks; ++c) {
      std::string augmented_pg_fn =
          wd + "/" + std::to_string(c) + ".augmented.pg";
      std::string gaf_fn = wd + "/" + std::to_string(c) + ".gaf";
      file_pairs.push_back({augmented_pg_fn, gaf_fn});
    }
    f.close();
    fprintf(stderr, "[M::%s] Loaded %ld chunks in %.3f sec\n", __func__,
            n_chunks, realtime() - rt);
  } else {
    if (chunk) {
      // === CHUNK AND CONVERT =========================================
      std::vector<std::set<std::string>> segments;
      {
        gbwtgraph::GBZ gbz;
        double rt0 = realtime();
        sdsl::simple_sds::load_from(gbz, gbz_fn);
        fprintf(stderr, "[M::%s] Load graph in %.3f sec\n", __func__,
                realtime() - rt0);

        // Chunking required 7m:35s and 32.8GB for HPRCv2
        rt0 = realtime();
        std::pair<std::vector<gbwtgraph::GBZ>, std::vector<std::string>>
            chunks = chunk_graph(gbz, {});
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

        // --- Convert .gbz to .pg -----------------------------------
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

          sub_gbz.graph.for_each_handle(
              [&](const handlegraph::handle_t &handle) {
                segments[i].insert(sub_gbz.graph.get_segment_name(handle));
              });

          if (delete_temp)
            std::filesystem::remove(sub_gbz_fn);

#pragma omp critical(printf_lock)
          {
            fprintf(stderr,
                    "[M::%s::%d] Processed chunk %ld (%s) in %.3f sec\n",
                    __func__, omp_get_thread_num(), i, chunks.second[i].c_str(),
                    realtime() - rt0);
            fflush(stderr);
          }
        }
      }

      // --- Split GAF -----------------------------------------------
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

      // --- augment chunks ------------------------------------------
      rt = realtime();
      {
#pragma omp parallel for num_threads(2) schedule(static, 1)
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

#pragma omp critical(printf_lock)
          {
            fprintf(stderr, "[M::%s::%d] Augmented chunk %ld in %.3f sec\n",
                    __func__, omp_get_thread_num(), c, realtime() - rt);
            fflush(stderr);
          }
          if (delete_temp)
            std::filesystem::remove(original_pg_fn);
        }
      }
      fprintf(stderr, "[M::%s] Augmented %ld chunks in %.3f sec\n", __func__,
              n_chunks, realtime() - rt);

      // ---

      std::ofstream out(done_fp, std::ios::out | std::ios::trunc);
      out << n_chunks << std::endl;
      out.close();
    }

    else {
      std::string pg_fn = wd + "/0.pg";
      std::string augmented_pg_fn = wd + "/0.augmented.pg";

      std::ostringstream cmd1;
      cmd1 << "vg convert --packed-out " << gbz_fn << " > " << pg_fn;

      // XXX: do this better
      (void)std::system(cmd1.str().c_str());

      std::ostringstream cmd2;
      cmd2 << "vg augment --include-paths --min-coverage 1 --gaf " << pg_fn
           << " " << gaf_fn << " > " << augmented_pg_fn;
      // XXX: do this better
      (void)std::system(cmd2.str().c_str());

      if (delete_temp)
        std::filesystem::remove(pg_fn);

      std::ofstream out(wd + "/DONE", std::ios::out | std::ios::trunc);
      out << 1 << std::endl;
      out.close();

      file_pairs.push_back({augmented_pg_fn, gaf_fn});
      fprintf(stderr, "[M::%s] Converted and augmented graph in %.3f sec\n",
              __func__, realtime() - rt);
    }
  }
  // =================================================================

  rt = realtime();
#pragma omp parallel for num_threads(nthreads) schedule(static, 1)
  for (size_t c = 0; c < file_pairs.size(); ++c) {
    // if (c != 10)
    //   continue;
    double rt0 = realtime();
    std::string pg_fn = file_pairs[c].first;
    if (chunk)
      gaf_fn = file_pairs[c].second;

    bdsg::PackedGraph *pg = new bdsg::PackedGraph();
    pg->deserialize(pg_fn);

    // size_t total_paths = pg->get_path_count();

    std::set<handlegraph::handle_t> known_vertices;
    std::set<handlegraph::edge_t> known_edges;
    int pi = 0;
    int np = 0;

    pg->for_each_path_handle([&](const bdsg::path_handle_t p) {
      std::string pname = pg->get_path_name(p);
      std::string_view sv = pname;
      // XXX: improve this, quite inefficient (batch insert?)
      if (sv.substr(0, 5).compare("palss") != 0) {
        handlegraph::handle_t last_h;
        bool first = true;

        for (const handlegraph::handle_t &h : pg->scan_path(p)) {
          known_vertices.insert(h);
          if (!first)
            known_edges.insert(std::make_pair(last_h, h));
          first = false;
          last_h = h;
        }
      } else {
        ++np;
      }
      ++pi;
    });

#pragma omp critical(printf_lock)
    {
      fprintf(stderr, "[M::%s::%i] Subgraph %ld - Known vertices/edges: %ld/%ld\n", __func__, omp_get_thread_num(), c,
             known_vertices.size(), known_edges.size());
      fflush(stderr);
    }

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

    size_t total_consensuses = 0;
    std::map<std::string, std::vector<handlegraph::handle_t>> real_novel;
    pg->for_each_path_handle([&](const bdsg::path_handle_t path) {
      std::string pname = pg->get_path_name(path);
      // if (pname.compare("palss-37356.0") != 0)
      //   return;
      std::string_view sv = pname;
      if (sv.substr(0, 5).compare("palss") == 0) {
        ++total_consensuses;
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
          for (size_t h = 1; h < rpath.size(); ++h) {
            n_novel_vertices +=
                known_vertices.find(rpath[h]) == known_vertices.end();
            n_novel_edges +=
                known_edges.find(std::make_pair(rpath[h - 1], rpath[h])) ==
                known_edges.end();
            rseq += pg->get_sequence(rpath[h]);
          }
          assert(seq.compare(rseq) == 0);

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
        }
      }
    });

#pragma omp critical(printf_lock)
    {
      fprintf(stderr, "[M::%s::%i] Subgraph %ld - Real novel vertices: %ld\n", __func__, omp_get_thread_num(), c, real_novel.size());
      fflush(stderr);
    }

    // 2-pass cleaning. If a vertex is supported by at least one path, we want
    // to keep it even if it's not supported by others
    std::map<handlegraph::handle_t,
             std::pair<std::set<std::string>, std::set<std::string>>>
        novel_vertices;
    std::map<handlegraph::edge_t,
             std::pair<std::set<std::string>, std::set<std::string>>>
        novel_edges;

    // get novel vertices and edges from real novel paths
    for (const auto &pp : real_novel) {
      std::string cname = pp.first;
      if (known_vertices.find(pp.second[0]) == known_vertices.end()) {
        novel_vertices[pp.second[0]].first.insert(cname);
        novel_vertices[pp.second[0]].second.insert(
            cluster_support[cname].begin(), cluster_support[cname].end());
      }
      for (size_t h = 1; h < pp.second.size(); ++h) {
        if (known_vertices.find(pp.second[h]) == known_vertices.end()) {
          novel_vertices[pp.second[h]].first.insert(cname);
          novel_vertices[pp.second[h]].second.insert(
              cluster_support[cname].begin(), cluster_support[cname].end());
        }
        if (known_edges.find(std::make_pair(pp.second[h - 1], pp.second[h])) ==
            known_edges.end()) {
          novel_edges[std::make_pair(pp.second[h - 1], pp.second[h])]
              .first.insert(cname);
          novel_edges[std::make_pair(pp.second[h - 1], pp.second[h])]
              .second.insert(cluster_support[cname].begin(),
                             cluster_support[cname].end());
        }
      }
    }

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
          path.push_back(h);
        }
        // } else {
        //   path = real_novel[pname];
        // }

        handlegraph::handle_t prev_handle;
        handlegraph::handle_t handle = path[0];

	// std::cerr << pg->get_id(handle) << " " << (known_vertices.find(handle) == known_vertices.end()) << " " << ((novel_vertices.find(handle) == novel_vertices.end()) ? -1 : novel_vertices.at(handle).second.size()) << std::endl;

        if (flag_vertex(handle, known_vertices, novel_vertices, min_supp) &&
            flag_vertex(pg->flip(handle), known_vertices, novel_vertices,
                        min_supp)) {
          vertices_to_remove.insert(handle);
          vertices_to_remove.insert(pg->flip(handle));
        }

        for (size_t h = 1; h < path.size(); ++h) {
          prev_handle = path[h - 1];
          handle = path[h];

	  // std::cerr << pg->get_id(handle) << " " << (known_vertices.find(handle) == known_vertices.end()) << " " << ((novel_vertices.find(handle) == novel_vertices.end()) ? -1 : novel_vertices.at(handle).second.size()) << std::endl;

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

#pragma omp critical(printf_lock)
    {
      fprintf(stderr, "[M::%s::%i] Subgraph %ld - Vertices/Edges to remove: %ld/%ld\n", __func__, omp_get_thread_num(), c, vertices_to_remove.size(), edges_to_remove.size());
      fflush(stderr);
    }

    std::set<std::string> retained_consensuses; // store the consensuses of
                                                // novel vertices/edges we print

    std::ofstream outfile(wd + "/" + std::to_string(c) + ".final.gfa");

    // === S LINES
    // =====================================================
    pg->for_each_handle(
        [&](const handlegraph::handle_t &handle) {
          // XXX: assuming here handle is always on + strand
	  // if (novel_vertices.find(handle) != novel_vertices.end())
	  //   std::cerr << "R: " << pg->get_id(handle) << " " << (vertices_to_remove.find(handle) == vertices_to_remove.end()) << std::endl;
          if (vertices_to_remove.find(handle) == vertices_to_remove.end()) {
            outfile << "S"
                    << "\t" << pg->get_id(handle) << "\t"
                    << pg->get_sequence(handle)
                    << (novel_vertices.find(handle) != novel_vertices.end()
                            ? "\tTY:Z:new"
                            : "")
                    << std::endl;
            retained_consensuses.insert(novel_vertices[handle].first.begin(),
                                        novel_vertices[handle].first.end());
          }
        },
        false);
#pragma omp critical(printf_lock)
    {
      fprintf(stderr, "[M::%s::%i] Subgraph %ld - Retained consensuses: %ld over %ld\n", __func__, omp_get_thread_num(), c, retained_consensuses.size(), total_consensuses);
      fflush(stderr);
    }

    // === L LINES
    // =====================================================
    pg->for_each_edge(
        [&](const handlegraph::edge_t &edge) {
          if (vertices_to_remove.find(edge.first) == vertices_to_remove.end() &&
              vertices_to_remove.find(edge.second) ==
                  vertices_to_remove.end() &&
              edges_to_remove.find(edge) == edges_to_remove.end()) {
            outfile << "L\t" << pg->get_id(edge.first) << "\t"
                    << (pg->get_is_reverse(edge.first) ? '-' : '+') << "\t"
                    << pg->get_id(edge.second) << "\t"
                    << (pg->get_is_reverse(edge.second) ? '-' : '+') << "\t"
                    << "0M"
                    << (known_edges.find(edge) != known_edges.end()
                            ? "\tTY:Z:new"
                            : "")
                    << std::endl;
            retained_consensuses.insert(novel_edges[edge].first.begin(),
                                        novel_edges[edge].first.end());
          }
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
        //   std::cerr << "[gfa] warning: incorrect end offset (" <<
        //   end_offset
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

    outfile.open(wd + "/" + std::to_string(c) + ".retained.gaf");
    std::ifstream file(gaf_fn);
    std::string line;
    while (std::getline(file, line)) {
      std::istringstream iss(line);
      std::string cname;
      std::getline(iss, cname, '\t');
      if (retained_consensuses.find(cname) != retained_consensuses.end())
        outfile << line << std::endl;
    }
    outfile.close();

    delete pg;
#pragma omp critical(printf_lock)
    {
      fprintf(stderr, "[M::%s::%d] Refined chunk %ld in %.3f sec\n", __func__, omp_get_thread_num(), c,
              realtime() - rt0);
      fflush(stderr);
    }
  }
  fprintf(stderr, "[M::%s] Refined in %.3f sec\n", __func__, realtime() - rt);

  // === UNCHOP ======================================================
  rt = realtime();
#pragma omp parallel for num_threads(2) schedule(static, 1)
  for (size_t c = 0; c < file_pairs.size(); ++c) {
    double rt0 = realtime();

    std::string gfa = wd + "/" + std::to_string(c) + ".final.gfa";
    std::string ugfa = wd + "/" + std::to_string(c) + ".final.unchopped.gfa";

    std::ostringstream cmd;
    cmd << "vg mod --unchop " << gfa << " > " << ugfa;
    // XXX: do this better
    (void)std::system(cmd.str().c_str());
    if (delete_temp) {
      std::filesystem::remove(wd + "/" + std::to_string(c) + ".final.gfa");
    }

#pragma omp critical(printf_lock)
    {
      fprintf(stderr, "[M::%s::%d] Unchopped chunk %ld in %.3f sec\n", __func__,
              omp_get_thread_num(), c, realtime() - rt0);
      fflush(stderr);
    }
  }
  fprintf(stderr, "[M::%s] Unchopped in %.3f sec\n", __func__, realtime() - rt);

  // === MERGE =======================================================
  rt = realtime();
  std::ofstream out_gaf;
  if (!retained_gaf_fn.empty())
    out_gaf.open(retained_gaf_fn, std::ios::out);

  std::cout << "H\tVN:Z:1.1" << std::endl;
  size_t min_idx = 1;
  for (size_t c = 0; c < file_pairs.size(); ++c) {
    double rt0 = realtime();
    std::ifstream in;
    std::map<std::string, std::string> mapping;
    {
      in.open(wd + "/" + std::to_string(c) + ".final.unchopped.gfa");

      std::string line;
      while (std::getline(in, line)) {
        if (line[0] == 'S') {
          std::string name;
          std::string seq;
          bool is_new = false;

          size_t start = 0;
          size_t i = 0;
          while (start <= line.size()) {
            std::size_t pos = line.find('\t', start);
            if (pos == std::string::npos)
              pos = line.size();
            if (i == 1) {
              name = line.substr(start, pos - start);
            } else if (i == 2) {
              seq = line.substr(start, pos - start);
            } else if (i == 3) {
              is_new = true;
            }
            ++i;

            start = pos + 1;

            if (pos == line.size())
              break;
          }

          mapping[name] = std::to_string(min_idx);
          ++min_idx;

          std::cout << "S"
                    << "\t" << mapping[name] << "\t" << seq
                    << (is_new ? "\tTY:Z:new" : "") << std::endl;
        } else if (line[0] == 'L') {
          std::vector<std::string> fields;
          // std::string name1;
          // std::string strand1;
          // std::string name2;
          // std::string strand2;
          // std::string cigar;

          size_t start = 0;
          // size_t i = 0;
          while (start <= line.size()) {
            std::size_t pos = line.find('\t', start);
            if (pos == std::string::npos)
              pos = line.size();
            fields.push_back(line.substr(start, pos - start));
            // ++i;
            start = pos + 1;
            if (pos == line.size())
              break;
          }

          fields[1] = mapping[fields[1]];
          fields[3] = mapping[fields[3]];
          std::cout << "L"
                    << "\t" << fields[1] << "\t" << fields[2] << "\t"
                    << fields[3] << "\t" << fields[4] << "\t" << fields[5]
                    << std::endl;
        } else if (line[0] == 'P') {
          rewrite_P_line(line, mapping);
        } else if (line[0] == 'W') {
          rewrite_W_line(line, mapping);
        }
      }
      in.close();
    }

    if (!retained_gaf_fn.empty()) {
      in.open(wd + "/" + std::to_string(c) + ".retained.gaf");
      out_gaf << in.rdbuf();
      in.close();
    }
    if (delete_temp) {
      std::filesystem::remove(wd + "/" + std::to_string(c) +
                              ".final.unchopped.gfa");
      std::filesystem::remove(wd + "/" + std::to_string(c) + ".retained.gaf");
    }

    fprintf(stderr, "[M::%s] Dumped chunk %ld in %.3f sec\n", __func__, c,
            realtime() - rt0);
  }
  out_gaf.close();
  std::cout.flush();

  fprintf(stderr, "[M::%s] Merged in %.3f sec\n", __func__, realtime() - rt);

  return 0;
}
