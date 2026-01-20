// #include <algorithm>
#include <getopt.h>
#include <iostream>
#include <string>
#include <vector>

extern "C" {
#include "abpoa.h"
#include "ksw2.h"
}

// #include "usage.h"

#include "graph.hpp"
#include "kmer.hpp"
#include "misc.hpp"
#include "sfs.hpp"

float canberra(const std::vector<uint32_t> &v1,
               const std::vector<uint32_t> &v2) {
  float cd = 0.0;
  for (size_t i = 0; i < v1.size(); ++i) {
    float x1 = (float)v1[i];
    float x2 = (float)v2[i];
    cd += x1 + x2 != 0 ? std::abs(x1 - x2) / (x1 + x2) : 0;
  }
  return cd;
}

// plain sequence to 0123-kmer counts
std::vector<uint32_t> count_kmers_plain(const char *seq, int seql, int klen) {
  std::vector<uint32_t> kcounts((1 << (2 * klen)), 0);
  int c; // current char
  uint8_t *kmer = (uint8_t *)malloc(klen);
  memcpy(kmer, seq, klen);
  uint64_t kmer_d = k2d((char *)kmer, klen);
  // uint64_t rckmer_d = rc(kmer_d, klen);
  // uint64_t ckmer_d = std::min(kmer_d, rckmer_d);
  ++kcounts[kmer_d];
  for (int p = klen; p < seql; ++p) {
    c = seq[p] < 5 ? seq[p] - 1 : rand() % 4;
    kmer_d = lsappend(kmer_d, c, klen);
    // rckmer_d = rsprepend(rckmer_d, reverse_char(c), klen);
    // ckmer_d = std::min(kmer_d, rckmer_d);
    ++kcounts[kmer_d];
  }
  free(kmer);
  return kcounts;
}

// XXX: improve/remove
std::string reverseAndComplement(const std::string &input) {
  std::string reversed(input.rbegin(), input.rend());
  std::string complement;

  for (char nucleotide : reversed) {
    switch (nucleotide) {
    case 'A':
      complement += 'T';
      break;
    case 'T':
      complement += 'A';
      break;
    case 'C':
      complement += 'G';
      break;
    case 'G':
      complement += 'C';
      break;
    default:
      complement += nucleotide; // handles non-nucleotide characters
      break;
    }
  }

  return complement;
}

int main_align(int argc, char *argv[]) {
  double rt;

  int klen = 31; // TODO: get k from somewhere else
  int cklen = 5;

  // int nth = 4;    // number of threads
  uint min_w = 2; // minimum size for clusters
  // int min_as = 0; // minimum alignment score
  bool out_sam = true;

  int _c;
  while ((_c = getopt(argc, argv, "k:w:gh")) != -1) {
    switch (_c) {
    case 'k':
      klen = std::stoi(optarg);
      break;
    case 'w':
      min_w = std::stoi(optarg);
      break;
    case 'g':
      out_sam = false;
      break;
      //     // case 's':
      //     //   min_as = std::stoi(optarg);
      //     //   break;
      //     // case '@':
      //     //   nth = std::stoi(optarg);
      //     //   break;
      //     case 'p':
      //       path_prefix = optarg;
      //       break;
    case 'h':
      // fprintf(stderr, "%s", AUGMENT_USAGE_MESSAGE);
      return 0;
    }
  }
  if (argc - optind != 2) {
    // fprintf(stderr, "%s", CALL_USAGE_MESSAGE);
    return 1;
  }
  std::string gbz_fn = argv[optind++];
  std::string sfs_fn = argv[optind++];

  rt = realtime();
  Graph graph(gbz_fn, "CHM13");
  graph.load();
  graph.load_fl();
  const gbwt::FastLocate &fl = graph.get_fl();

  std::vector<sfs_t> specific_strings = load_sfs(sfs_fn);
  fprintf(stderr, "[M::%s] Loaded anchored specific strings in %.3f sec\n",
          __func__, realtime() - rt);

  // Cluster specific strings based on their anchors (anchors have been
  // reversed)
  rt = realtime();
  std::sort(specific_strings.begin(), specific_strings.end(),
            [](const sfs_t &a, const sfs_t &b) { return a.skmer < b.skmer; });
  fprintf(stderr, "[M::%s] Sorted specific strings in %.3f sec\n", __func__,
          realtime() - rt);

  rt = realtime();
  std::vector<std::vector<sfs_t>> clusters;
  if (specific_strings.size() > 0)
    clusters.push_back({specific_strings[0]});
  uint last_c = 0;
  for (uint ss = 1; ss < specific_strings.size(); ++ss) {
    const sfs_t &sfs = specific_strings[ss];

    if (sfs.skmer == clusters[last_c][0].skmer) {
      uint c = last_c;
      for (c = last_c; c < clusters.size(); ++c) {
        if (sfs.ekmer == clusters[c][0].ekmer)
          break;
      }
      if (c < clusters.size())
        clusters[c].push_back(sfs);
      else
        clusters.push_back({sfs});
    } else {
      clusters.push_back({sfs});
      last_c = clusters.size() - 1;
    }
  }
  fprintf(stderr, "[M::%s] Built %ld clusters in %.3f sec\n", __func__,
          clusters.size(), realtime() - rt);

  // Initialize abpoa
  abpoa_t *ab = abpoa_init();
  abpoa_para_t *abpt = abpoa_init_para();
  abpt->align_mode = 0; // global
  abpt->out_msa = 0;
  abpt->out_cons = 1;
  abpt->out_gfa = 0;
  abpt->disable_seeding = 1;
  abpt->progressive_poa = 0;
  abpt->amb_strand = 0;
  // abpt->wb = -1;
  abpt->max_n_cons = 2;
  abpt->min_freq = 0.25;
  abpoa_post_set_para(abpt);

  // XXX: Assuming clusters of size <= 64. TODO: reallocation
  // alternatively, we could create a struct cluster_t with the array
  uint8_t **cluster_seqs = (uint8_t **)malloc(sizeof(uint8_t *) * 64);
  int *cluster_seqs_lens = (int *)malloc(sizeof(int) * 64);

  // Initialize ksw2
  // https://github.com/lh3/minimap2/blob/69e36299168d739dded1c5549f662793af10da83/options.c#L144
  // asm5
  // int sc_mch = 1, sc_mis = -19, gapo = 39, gape = 3, gapo2 = 81, gape2 = 1;
  // asm10
  int sc_mch = 1, sc_mis = -9, gapo = 16, gape = 2, gapo2 = 41, gape2 = 1;
  int8_t a = (int8_t)sc_mch,
         b = sc_mis < 0 ? (int8_t)sc_mis : -(int8_t)sc_mis; // a>0 and b<0
  int8_t mat[25] = {a, b, b, b, 0, b, a, b, b, 0, b, b, a,
                    b, 0, b, b, b, a, 0, 0, 0, 0, 0, 0};
  ksw_extz_t ez;
  memset(&ez, 0, sizeof(ksw_extz_t));

  // TODO: parallelize

  rt = realtime();
  for (size_t cidx = 0; cidx < clusters.size(); ++cidx) {
    std::vector<sfs_t> &cluster = clusters[cidx];

    // if (cidx != 33)
    //   continue;

    if ((cidx + 1) % 1000 == 0) {
      fprintf(stderr, "[M::%s] analyzed %ld/%ld clusters in %.3f sec\n",
              __func__, cidx + 1, clusters.size(), realtime() - rt);
    }

    if (cluster.size() > 64) {
      // FIXME
      continue;
    }

    const sfs_t &s0 = cluster.front();

    fprintf(stderr, "Cluster %ld: %s>%s\n", cidx,
            graph.get_gfa_name(s0.sv >> 1).c_str(),
            graph.get_gfa_name(s0.ev >> 1).c_str());

    // is it enough to consider paths from s0 only?
    std::map<uint64_t, uint64_t> path_map;
    for (const uint64_t &p : s0.paths) {
      path_map[p >> 4] = p;
    }

    uint32_t sv = s0.sv;
    uint32_t sv_inv = sv ^ 1;
    uint32_t soff = s0.soff;
    uint32_t soff_inv = graph.get_vertex_len(sv >> 1) - soff - 1;
    uint32_t ev = s0.ev;
    uint32_t ev_inv = ev ^ 1;
    uint32_t eoff = s0.eoff;
    uint32_t eoff_inv = graph.get_vertex_len(ev >> 1) - eoff - 1;

    uint64_t p0 = s0.paths.front();
    bool p0_strand = (p0 >> 2) & 1;
    bool p0_sinv = (p0 >> 1) & 1;
    bool p0_einv = p0 & 1;

    // Find paths going through sv and ev (both versions)
    std::vector<path_t> paths1, paths2;

    if (p0_strand) {
      paths1 = graph.get_paths(p0_sinv ? sv_inv : sv, p0_einv ? ev_inv : ev, 1,
                               false);
      paths2 = graph.get_paths(p0_einv ? ev : ev_inv, p0_sinv ? sv : sv_inv, 1,
                               false);
    } else {
      paths1 = graph.get_paths(p0_einv ? ev_inv : ev, p0_sinv ? sv_inv : sv, 1,
                               false);
      paths2 = graph.get_paths(p0_sinv ? sv : sv_inv, p0_einv ? ev : ev_inv, 1,
                               false);
    }
    assert(paths1.size() > 0 || paths2.size() > 0);

    std::vector<path_t> paths;
    for (path_t &path : paths1) {
      if (path_map.find(path.id) == path_map.end())
        continue;
      paths.push_back(path);
      paths.back().strand = p0_strand;
    }
    for (path_t &path : paths2) {
      if (path_map.find(path.id) == path_map.end())
        continue;
      paths.push_back(path);
      paths.back().strand = !p0_strand;
    }

    // Keep only unique paths
    std::sort(paths.begin(), paths.end());
    paths.erase(std::unique(paths.begin(), paths.end()), paths.end());

    // TODO: keep track of all paths, even erased ones

    std::vector<std::vector<uint32_t>> path_kcounts(paths.size());
    size_t pid = 0;

    for (auto &path : paths) {
      uint64_t p = path_map[path.id];

      bool strand = (p >> 2) & 1;
      assert(strand == path.strand);

      bool sinv = (p >> 1) & 1;
      bool einv = p & 1;

      uint32_t v1 = sv;
      uint32_t off1 = soff;
      uint32_t v2 = ev;
      uint32_t off2 = eoff;

      if (sinv) {
        v1 = sv_inv;
        off1 = soff_inv;
      }
      if (einv) {
        v2 = ev_inv;
        off2 = eoff_inv;
      }

      if (strand == 0) {
        std::swap(v1, v2);
        std::swap(off1, off2);
      }

      /*
       * If on + strand, path goes from (v1, s1) to (v2, s2), on - strand, it
       * goes from (v2, ~s2) to (v1, ~s1)
       */

      assert(v1 == path.vertices.front() && v2 == path.vertices.back());

      // path start/end positions (both inclusive)
      path.ps = off1;
      path.pe = path.sequence.size() - graph.get_vertex_len(v2 >> 1) + off2;
      path.l = path.sequence.size();
      path.sequence = path.sequence.substr(path.ps, path.pe - path.ps + 1);

      // XXX: why do we need to do this only when strands differ?
      path.reversed = p0_strand != strand;

      // When needed, reverse path, count kmers (since we don't want to use
      // canonical kmers), then reverse again
      if (path.reversed)
        path.sequence = reverseAndComplement(path.sequence);
      path_kcounts[pid] =
          count_kmers_plain(path.sequence.c_str(), path.sequence.size(), cklen);
      if (path.reversed)
        path.sequence = reverseAndComplement(path.sequence);
      ++pid;
    }

    // Consensus via abpoa
    int goods = 0;
    for (sfs_t &s : cluster) {
      cluster_seqs_lens[goods] = s.l;
      cluster_seqs[goods] = s.seq;
      ++goods;
    }
    abpoa_msa(ab, abpt, goods, NULL, cluster_seqs_lens, cluster_seqs, NULL,
              NULL);

    abpoa_cons_t *abc = ab->abc;
    // total_clusters += abc->n_cons;
    for (int ci = 0; ci < abc->n_cons; ++ci) {
      int cons_l = abc->cons_len[ci];
      int abpoa_supp = abc->clu_n_seq[ci];
      char *cons_seq = (char *)abc->cons_base[ci];
      for (int i = 0; i < cons_l; ++i)
        cons_seq[i] = "ACGT"[(int)cons_seq[i]];
      cons_seq[cons_l] = '\0';
      // std::cerr << "= " << cons_seq << std::endl;

      // === select best path ===
      std::vector<uint32_t> consensus_kcounts =
          count_kmers_plain(cons_seq, cons_l, cklen);
      size_t best_p = -1U;
      float cd = 0, best_cd = INT64_MAX;
      for (size_t i = 0; i < path_kcounts.size(); ++i) {
        const auto &pkc = path_kcounts[i];
        if (pkc.empty())
          continue;
        cd = canberra(consensus_kcounts, pkc);
        if (cd < best_cd) {
          best_cd = cd;
          best_p = i;
        }
      }
      assert(best_p != -1U);

      const path_t &path = paths[best_p];
      // std::cerr << (path.reversed ? "<" : ">") << " " << path.sequence
      //           << std::endl;

      std::cerr << graph.get_path_contig(path.id >> 1) << " "
                << graph.get_gfa_name(path.vertices.front() >> 1)
                << ((path.vertices.front() & 1) ? "-" : "+") << " > "
                << graph.get_gfa_name(path.vertices.back() >> 1)
                << ((path.vertices.back() & 1) ? "-" : "+") << std::endl;

      std::string cseq_plain(cons_seq);
      if (path.reversed)
        cseq_plain = reverseAndComplement(cseq_plain);
      std::string pseq_plain(path.sequence);
      // go back to 0123
      for (int i = 0; i < cons_l; ++i)
        cons_seq[i] = to_int[(uint8_t)cseq_plain[i]] - 1;

      char *path_seq = (char *)malloc(pseq_plain.size() + 1);
      for (size_t i = 0; i < pseq_plain.size(); ++i)
        path_seq[i] = to_int[(uint8_t)pseq_plain[i]] - 1;
      path_seq[pseq_plain.size()] = '\0';

      // if (std::abs((int)pseq_plain.size() - cons_l) >= 50000) {
      //   free(path_seq);
      //   continue;
      // }

      fprintf(stderr, "Aligning %dbp against %ldbp\n", cons_l,
              path.sequence.size());
      ksw_extd2_sse(0, cons_l, (uint8_t *)cons_seq, path.sequence.size(),
                    (uint8_t *)path_seq, 5, mat, gapo, gape, gapo2, gape2, -1,
                    -1, -1, 0, &ez);

      // OUTPUT
      if ((ez.cigar[0] & 0xf) != 0 || (ez.cigar[ez.n_cigar - 1] & 0xf) != 0) {
        // clipped = 1;
        free(path_seq);
        continue; // XXX: do we want these?
      }

      // Parse CIGAR
      int opl;
      int tot_cigar_len = 0;
      int tot_res_matches = 0;

      // XXX: improve this
      std::string cigar;
      std::string cs;
      int cons_p = 0;
      int pseq_p = 0;
      for (int i = 0; i < ez.n_cigar; ++i) {
        opl = ez.cigar[i] >> 4;
        tot_cigar_len += opl;

        // cigar
        cigar += std::to_string(opl) + "MID"[ez.cigar[i] & 0xf];

        // difference string, residues, cigar length
        if ((ez.cigar[i] & 0xf) == 0) {
          // M
          int l_tmp = 0;
          for (int j = 0; j < opl; ++j) {
            if (cons_seq[cons_p + j] != path_seq[pseq_p + j]) {
              if (l_tmp > 0) {
                cs += ":" + std::to_string(l_tmp);
                l_tmp = 0;
              }
              cs += "*";
              cs += "ACGT"[(uint8_t)path_seq[pseq_p + j]];
              cs += "ACGT"[(uint8_t)cons_seq[cons_p + j]];
            } else {
              ++l_tmp;
            }
          }
          if (l_tmp > 0) {
            tot_res_matches += l_tmp;
            cs += ":" + std::to_string(l_tmp);
          }
          cons_p += opl;
          pseq_p += opl;
        } else if ((ez.cigar[i] & 0xf) == 1) {
          // I
          std::string sub(cons_seq + cons_p, opl);
          for (size_t i = 0; i < sub.size(); ++i)
            sub[i] = "ACGT"[(uint8_t)sub[i]];
          cs += "+" + sub;
          cons_p += opl;
        } else if ((ez.cigar[i] & 0xf) == 2) {
          // D
          std::string sub(path_seq + pseq_p, opl);
          for (size_t i = 0; i < sub.size(); ++i)
            sub[i] = "ACGT"[(uint8_t)sub[i]];
          cs += "-" + sub;
          pseq_p += opl;
        } else {
          fprintf(stderr,
                  "Cluster %ld --- We shouldn't be here while parsing ksw "
                  "cigar. Halting.\n",
                  cidx);
          exit(EXIT_FAILURE);
        }
      }

      for (int i = 0; i < cons_l; ++i)
        cons_seq[i] = "ACGT"[(uint8_t)cons_seq[i]];
      cons_seq[cons_l] = '\0';
      for (size_t i = 0; i < path.sequence.size(); ++i)
        path_seq[i] = "ACGT"[(uint8_t)path_seq[i]];
      path_seq[path.sequence.size()] = '\0';

      // since vertices might be split (if longer than 1024), the offset along
      // path refers to "internal" vertices, not GFA segments
      gbwtgraph::handle_t handle =
          gbwtgraph::GBWTGraph::node_to_handle(path.vertices[0]);
      size_t xx = graph.get_gbz().graph.get_segment_offset(handle);
      size_t ps = xx + path.ps;
      size_t pe = ps + path.sequence.size();
      //

      // print GAF line
      std::cout << cidx << "." << ci << "\t";
      std::cout << cons_l << "\t";
      std::cout << 0 << "\t";
      std::cout << cons_l << "\t";
      std::cout << (path.reversed ? "-" : "+") << "\t";
      //
      size_t x = 0;
      std::string segment_name = graph.get_gfa_name(path.vertices[x] >> 1);
      std::cout << ((path.vertices[x] & 1) ? "<" : ">") << segment_name;
      ++x;
      std::string last_segment_name = segment_name;
      for (; x < path.vertices.size(); ++x) {
        segment_name = graph.get_gfa_name(path.vertices[x] >> 1);
        if (segment_name.compare(last_segment_name) != 0)
          std::cout << ((path.vertices[x] & 1) ? "<" : ">") << segment_name;
        last_segment_name = segment_name;
      }
      std::cout << "\t";
      //
      std::cout << path.l << "\t";
      std::cout << ps << "\t";
      std::cout << pe << "\t";
      std::cout << tot_res_matches << "\t";
      std::cout << tot_cigar_len << "\t";
      std::cout << 60 << "\t";
      std::cout << "AS:i:" << ez.score << "\t";
      std::cout << "cg:Z:" << cigar << "\t";
      std::cout << "cs:Z:" << cs << "\t";
      // std::cout << "cl:Z:" << clipped << "\t";
      std::cout << "cw:Z:" << abpoa_supp << "\t";
      std::cout << "qs:Z:" << cseq_plain << "\t";
      std::cout << "ps:Z:" << pseq_plain;
      std::cout << std::endl;

      free(path_seq);
    }
  }

  free(ez.cigar);
  abpoa_free_para(abpt);
  abpoa_free(ab);
  free(cluster_seqs);
  free(cluster_seqs_lens);

  for (sfs_t &s : specific_strings)
    free(s.seq);

  return 0;
}
