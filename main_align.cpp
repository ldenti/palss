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

std::vector<uint32_t> count_kmers(const char *seq, int seql, int klen) {
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

  std::vector<sfs_t> specific_strings = load_anchored_sfs(sfs_fn);
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
//
https: // github.com/lh3/minimap2/blob/69e36299168d739dded1c5549f662793af10da83/options.c#L144
  // asm5
  // int sc_mch = 1, sc_mis = -19, gapo = 39, gape = 3, gapo2 = 81, gape2 =
  1;
  // asm10
  int sc_mch = 1, sc_mis = -9, gapo = 16, gape = 2, gapo2 = 41, gape2 = 1;
  int8_t a = (int8_t)sc_mch,
         b = sc_mis < 0 ? (int8_t)sc_mis : -(int8_t)sc_mis; // a>0 and b<0
  int8_t mat[25] = {a, b, b, b, 0, b, a, b, b, 0, b, b, a,
                    b, 0, b, b, b, a, 0, 0, 0, 0, 0, 0};
  ksw_extz_t ez;
  memset(&ez, 0, sizeof(ksw_extz_t));

  std::map<std::string, size_t> ref_paths;
  if (out_sam) {
    std::cout << "@HD\tVN:1.6\tSO:coordinate" << std::endl;
    ref_paths = graph.get_reference_paths();
    for (const auto &[name, l] : ref_paths)
      std::cout << "@SQ\tSN:" << name << "\tLN:" << l << std::endl;
  }

  // TODO: parallelize
  int small_clusters = 0;
  int big_clusters = 0;
  int nopaths_clusters = 0;
  int noanchors_clusters = 0;
  int strange_clusters = 0;

  rt = realtime();
  for (size_t cidx = 0; cidx < clusters.size(); ++cidx) {
    std::vector<sfs_t> &cluster = clusters[cidx];

    if (cidx != 903 && cidx != 324 && cidx != 352)
      continue;

    if ((cidx + 1) % 1000 == 0) {
      fprintf(stderr, "[M::%s] analyzed %ld/%ld clusters in %.3f sec\n",
              __func__, cidx + 1, clusters.size(), realtime() - rt);
    }

    if (cluster.size() < min_w) {
      ++small_clusters;
      continue;
    }
    if (cluster.size() > 64) {
      // FIXME
      ++big_clusters;
      continue;
    }

    int v1 = cluster.front().sv1;
    int v1b = cluster.front().sv2;
    int v2 = cluster.front().ev2;
    int v2a = cluster.front().ev1;

    fprintf(stderr, "Cluster %ld: %s/%s > %s/%s\n", cidx,
            graph.get_gfa_name(v1).c_str(), graph.get_gfa_name(v1b).c_str(),
            graph.get_gfa_name(v2a).c_str(), graph.get_gfa_name(v2).c_str());
      for (const sfs_t &s : cluster) {
        std::cerr << s.rname << std::endl;
      }


    // TODO: store path names
    std::vector<path_t> paths = graph.get_paths(v1, v2, out_sam);
    size_t npaths = paths.size();
    if (npaths == 0) {
      ++nopaths_clusters;
      continue;
    }
    if (npaths > 1) {
      std::sort(paths.begin(), paths.end());
      // paths.erase(std::unique(paths.begin(), paths.end()), paths.end());
    }
    size_t nupaths = paths.size();
    // for (path_t &p : paths) {
    //   // CHECKME: it seems we are not supposed to reverse&complement the
    //   paths, not sure why if (p.id & 1) {
    //     p_reverse(p);
    //   }
    // }

    // fprintf(stderr, "%ld/%ld paths\n", nupaths, npaths);

    // get starting/ending kmers from specific strings
    std::string first_kmer_1, last_kmer_1;
    std::string first_kmer_2, last_kmer_2;
    for (const auto &s : cluster) {
      std::string kmer(s.plain_seq, 0, klen);
      std::string last_kmer(s.plain_seq, s.l - klen, klen);
      if (first_kmer_1.empty()) {
        first_kmer_1 = kmer;
        last_kmer_1 = last_kmer;
        continue;
      }
      if (kmer.compare(first_kmer_1) == 0)
        continue;
      first_kmer_2 = kmer;
      last_kmer_2 = last_kmer;
      break;
    }
    if (first_kmer_2.empty() && last_kmer_2.empty()) {
      // we had reads on one strand only, reverse&complement the kmers anyway
      // since paths may be on other strand
      first_kmer_2 = reverseAndComplement(last_kmer_1);
      last_kmer_2 = reverseAndComplement(first_kmer_1);
    }

    std::cerr << first_kmer_1 << " " << last_kmer_1 << std::endl;
    std::cerr << first_kmer_2 << " " << last_kmer_2 << std::endl;

    bool hit = false;
    std::string first_kmer, last_kmer;
    for (const path_t &path : paths) {
      // std::cerr << path.sequence << std::endl;
      // std::cerr << (path.id & 1) << std::endl;

      size_t p1, p2;
      if ((p1 = path.sequence.find(first_kmer_1)) != std::string::npos &&
          (p2 = path.sequence.find(last_kmer_1)) != std::string::npos) {
        hit = true;
        if (p1 < p2) {
          first_kmer = first_kmer_1;
          last_kmer = last_kmer_1;
          break; // CHECKME: is checking one path enough?
        }

      } else if ((p1 = path.sequence.find(first_kmer_2)) != std::string::npos &&
                 (p2 = path.sequence.find(last_kmer_2)) != std::string::npos) {
        hit = true;
        if (p1 < p2) {
          first_kmer = first_kmer_2;
          last_kmer = last_kmer_2;
          break; // CHECKME: is checking one path enough?
        }
      }
    }

    if (first_kmer.empty() && last_kmer.empty()) {
      // XXX: anchors may not follow path sequence since in search we do not
      // filter paths with inverted anchors
      // FIXME: can we fix this while anchoring somehow?
      ++noanchors_clusters;
      /**
      fprintf(stderr,
              "Skipping cluster %ld since could not find anchors on any path, "
              "kmers: %d\n",
              cidx, hit);
      **/
      continue;
    }

    std::cerr << " -> " << first_kmer << " " << last_kmer << std::endl;

    // cut path sequences and compute kmer counts
    std::vector<std::vector<uint32_t>> path_kcounts(paths.size());
    for (size_t p = 0; p < paths.size(); ++p) {
      path_t &path = paths[p];
      size_t p1 = path.sequence.find(first_kmer);
      size_t p2 = path.sequence.find(last_kmer);
      if (p1 == std::string::npos || p2 == std::string::npos)
        // not all paths are good, since there might be bubbles
        continue;
      assert(p1 < p2);
      path.cutpfx = p1;
      path.cutsfx = path.sequence.size() - (p2 + klen);
      path.sequence = path.sequence.substr(p1, p2 + klen - p1);
      path_kcounts[p] =
          count_kmers(path.sequence.c_str(), path.sequence.size(), cklen);
    }

    for (auto &s : cluster) {
      if (s.plain_seq.substr(0, 31).compare(first_kmer) != 0) {
        int i;
        for (i = 0; i < (s.l >> 1); ++i) {
          int tmp = s.seq[s.l - 1 - i];
          s.seq[s.l - 1 - i] = 3 - s.seq[i];
          s.seq[i] = 3 - tmp;
        }
        s.seq[i] = 3 - s.seq[i];

        std::string seq;
        for (i = 0; i < s.l; ++i)
          seq += "ACGT"[s.seq[i]];
        s.plain_seq = seq;
      }
      assert(s.plain_seq.substr(0, 31).compare(first_kmer) == 0 &&
             s.plain_seq.substr(s.l - 31, 31).compare(last_kmer) == 0);
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
      // for (int i = 0; i < cons_l; ++i)
      //   cons_seq[i] = "ACGT"[(int)cons_seq[i]];

      // select best path
      std::vector<uint32_t> consensus_kcounts =
          count_kmers(cons_seq, cons_l, cklen);
      size_t best_p = -1U;
      float cd, best_cd = INT64_MAX;
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
      const path_t &path = paths[best_p];
      char *path_seq = (char *)malloc(path.sequence.size());
      for (size_t i = 0; i < path.sequence.size(); ++i)
        path_seq[i] = to_int[(uint8_t)path.sequence[i]] - 1;

      // for (int i = 0; i < cons_l; ++i)
      //   std::cerr << (int)cons_seq[i];
      // std::cerr << std::endl;
      // for (int i = 0; i < path.sequence.size(); ++i)
      //   std::cerr << (int)path_seq[i];
      // std::cerr << std::endl;

      if (std::abs((int)path.sequence.size() - cons_l) >= 50000) {
        ++strange_clusters;
        fprintf(
            stderr,
            "Skipping cluster %ld.%d since path/consensus lengths disagree (%ld/%d)\n",
            cidx, ci, path.sequence.size(), cons_l);
        continue;
      }
      // fprintf(stderr, "Aligning %dbp against %ldbp\n", cons_l,
      //         path.sequence.size());
      ksw_extd2_sse(0, cons_l, (uint8_t *)cons_seq, path.sequence.size(),
                    (uint8_t *)path_seq, 5, mat, gapo, gape, gapo2, gape2, -1,
                    -1, -1, 0, &ez);

      // OUTPUT
      if ((ez.cigar[0] & 0xf) != 0 || (ez.cigar[ez.n_cigar - 1] & 0xf) != 0) {
        // clipped = 1;
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

      if (out_sam) {
        std::string contig_name = graph.get_path_contig(path.id);
        std::cout << cidx << "." << ci << "\t";
        std::cout << 0 << "\t";
        std::cout << contig_name << "\t";
        std::cout << ref_paths[contig_name] - path.offset1 + path.cutpfx + 1
                  << "\t";
        std::cout << 60 << "\t";
        std::cout << cigar << "\t";
        std::cout << "*"
                  << "\t";
        std::cout << 0 << "\t";
        std::cout << 0 << "\t";
        std::cout << cons_seq << "\t";
        std::cout << "*"
                  << "\t";
        std::cout << "cs:Z:" << cs << "\t";
        std::cout << "cw:Z:" << abpoa_supp << "\t";
        std::cout << "AS:i:" << ez.score << "\t";
        // std::cout << "ps:Z:" << pseq;
        std::cout << std::endl;
      } else {
        // print GAF line
        std::cout << cidx << "." << ci << "\t";
        std::cout << cons_l << "\t";
        std::cout << 0 << "\t";
        std::cout << cons_l << "\t";
        std::cout << "+\t";
        std::cout << ">" << graph.get_gfa_name(path.vertices[0]);
        for (size_t i = 1; i < path.vertices.size(); ++i)
          std::cout << ">" << graph.get_gfa_name(path.vertices[i]);
        std::cout << "\t";
        std::cout << path.sequence.size() + path.cutpfx + path.cutsfx << "\t";
        std::cout << path.cutpfx << "\t";
        std::cout << path.cutpfx + path.cutpfx << "\t";
        std::cout << tot_res_matches << "\t";
        std::cout << tot_cigar_len << "\t";
        std::cout << 60 << "\t";
        std::cout << "AS:i:" << ez.score << "\t";
        std::cout << "cg:Z:" << cigar << "\t";
        std::cout << "cs:Z:" << cs << "\t";
        // std::cout << "cl:Z:" << clipped << "\t";
        std::cout << "cw:Z:" << abpoa_supp << "\t";
        // std::cout << "qs:Z:" << consensus[idx] << "\t";
        // std::cout << "ps:Z:" << pseq;
        std::cout << std::endl;
      }

      free(path_seq);
      // return 1;
    }
  }

  free(ez.cigar);
  abpoa_free_para(abpt);
  abpoa_free(ab);
  free(cluster_seqs);
  free(cluster_seqs_lens);

  for (sfs_t &s : specific_strings)
    free(s.seq);

  fprintf(stderr, "Total clusters: %ld\n", clusters.size());
  fprintf(stderr, "Small clusters: %d (%f)\n", small_clusters,
          (float)small_clusters / clusters.size());
  fprintf(stderr, "Big clusters: %d (%f)\n", big_clusters,
          (float)big_clusters / clusters.size());
  fprintf(stderr, "NoPaths clusters: %d (%f)\n", nopaths_clusters,
          (float)nopaths_clusters / clusters.size());
  fprintf(stderr, "NoAnchors clusters: %d (%f)\n", noanchors_clusters,
          (float)noanchors_clusters / clusters.size());
  fprintf(stderr, "Strange clusters: %d (%f)\n", strange_clusters,
          (float)strange_clusters / clusters.size());
  return 0;
}
