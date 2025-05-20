#include <string>
#include <vector>
#include <zlib.h>

extern "C" {
#include "abpoa.h"
#include "io.h"
#include "kseq.h"
#include "ksw2.h"
#include "sfs.h"
}
#include "graph.hpp"
// #include "usage.hpp"

KSEQ_INIT(gzFile, gzread)

std::string decode(const uint8_t *s, int l) {
  if (s == NULL)
    return "";
  char ds[l + 1];
  for (int i = 0; i < l; ++i)
    ds[i] = s[i] <= 3 ? "ACGT"[s[i]] : 'N';
  ds[l] = '\0';
  return ds;
}

int main_call(int argc, char *argv[]) {
  double rt = realtime(), rt1;

  int klen = 27; // kmer size
  int min_w = 2; // minimum support for cluster
  // float lr = 0.97; // Ratio to cluster specific strings inside same cluster

  int _c;
  while ((_c = getopt(argc, argv, "k:h")) != -1) {
    switch (_c) {
    case 'k':
      klen = std::stoi(optarg);
      break;
    case 'h':
      fprintf(stderr, "%s", "CALL_USAGE_MESSAGE");
      return 0;
    default:
      fprintf(stderr, "Error\n");
      return 1;
    }
  }
  if (argc - optind != 3) {
    fprintf(stderr, "Error!\n");
    return 1;
  }

  std::string gbz_fn = argv[optind++];
  std::string sfs_fn = argv[optind++];
  std::string fq_fn = argv[optind++];

  // Graph
  Graph graph(gbz_fn);
  graph.load();
  // graph.print_stats();
  positions_t positions = graph.get_positions();

  // cluster by starting/ending solid anchors
  int nreads = 0;
  FILE *fp;
  char *line = NULL;
  size_t len = 0;
  ssize_t read;
  std::vector<sfs_t> specifics;
  fp = fopen(sfs_fn.c_str(), "r");
  if (fp == NULL)
    exit(EXIT_FAILURE);
  while ((read = getline(&line, &len, fp)) != -1) {
    if (line[0] != 'O')
      continue;
    sfs_t ss = parse_sfs_line(line + 2);
    if (ss.qidx > nreads)
      nreads = ss.qidx;
    specifics.push_back(ss);
  }
  ++nreads; // since qidx starts from 0
  fclose(fp);

  std::sort(specifics.begin(), specifics.end(),
            [](const sfs_t &a, const sfs_t &b) { return a.a.seq < b.a.seq; });

  std::vector<std::vector<sfs_t>> clusters;
  clusters.push_back({specifics[0]});
  int last_c = 0;
  for (int ss = 1; ss < specifics.size(); ++ss) {
    sfs_t sfs = specifics[ss];
    if (sfs.a.seq == clusters[last_c][0].a.seq) {
      int c = last_c;
      for (c = last_c; c < clusters.size(); ++c) {
        if (sfs.b.seq == clusters[c][0].b.seq)
          break;
      }
      if (c < clusters.size()) {
        clusters[c].push_back(sfs);
      } else {
        clusters.push_back({sfs});
      }
    } else {
      clusters.push_back({sfs});
      last_c = clusters.size() - 1;
    }
  }

  // Split specific strings by read
  std::vector<std::vector<sfs_t *>> ss_byread(nreads);
  for (auto &cluster : clusters) {
    if (cluster.size() < min_w)
      continue;
    for (sfs_t &sfs : cluster) {
      ss_byread[sfs.qidx].push_back(&sfs);
    }
  }

  // Fill specific strings sequence (from reads)
  gzFile fq_fp = gzopen(fq_fn.c_str(), "r");
  kseq_t *seq = kseq_init(fq_fp);
  uint8_t *s;
  int l;
  int qidx = 0;
  std::vector<std::string> qnames;
  while ((l = kseq_read(seq)) >= 0) {
    qnames.push_back(seq->name.s);
    if (qidx >= ss_byread.size() || ss_byread[qidx].empty()) {
      ++qidx;
      continue;
    }

    s = (uint8_t *)seq->seq.s;
    rb3_char2nt6(seq->seq.l, s);
    if (ss_byread[qidx][0]->strand == 0) {
      // Reverse and complement the sequence
      int i;
      for (i = 0; i < (l >> 1); ++i) {
        int tmp = s[l - 1 - i];
        tmp = (tmp >= 1 && tmp <= 4) ? 5 - tmp : tmp;
        s[l - 1 - i] = (s[i] >= 1 && s[i] <= 4) ? 5 - s[i] : s[i];
        s[i] = tmp;
      }
      if (l & 1)
        s[i] = (s[i] >= 1 && s[i] <= 4) ? 5 - s[i] : s[i];
    }
    for (sfs_t *ss : ss_byread[qidx]) {
      ss->seq = (uint8_t *)malloc((ss->l + 1) * sizeof(uint8_t));
      memcpy(ss->seq, seq->seq.s + ss->s, ss->l);
      ss->seq[ss->l] = '\0';
      // XXX: convert to 0123 since it seems abpoa works better this way in
      // diploid mode
      for (int i = 0; i < ss->l; ++i)
        --ss->seq[i];
    }
    ++qidx;
  }
  kseq_destroy(seq);
  gzclose(fq_fp);

  // --- Analyze clusters

  // === INIT ABPOA
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

  // XXX: Assuming clusters of size <= 64
  uint8_t **cseqs = (uint8_t **)malloc(sizeof(uint8_t *) * 64);
  int *cseqs_lens = (int *)malloc(sizeof(int) * 64);
  int goods = 0;

  // === INIT KSW2
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

  int cidx = 0;
  rt = realtime();
  for (auto &cluster : clusters) {
    ++cidx;
    if (cidx % 1000 == 0) {
      fprintf(stderr, "[M::%s] analyzed %d/%ld clusters in %.3f sec\n",
              __func__, cidx, clusters.size(), realtime() - rt);
      rt = realtime();
    }
    if (cluster.size() < min_w)
      continue;
    if (cluster.size() > 63)
      continue;

    // std::cerr << cluster.size() << std::endl;

    // get reference locus
    sfs_t ss = cluster[0];
    int v1 = ss.a.v;
    int v2 = ss.b.v;
    std::string ref1 = positions.v2ref[v1];
    std::string ref2 = positions.v2ref[v2];
    if (ref1.compare(ref2) != 0)
      continue;
    if (positions.offsets.find(v1) == positions.offsets.end() ||
        positions.offsets.find(v2) == positions.offsets.end()) {
      // TODO: find best path by intersecting paths passing trough the two
      // vertices and report over it
      continue;
    }
    int pos1 = positions.offsets.at(v1) + ss.a.offset;
    int pos2 = positions.offsets.at(v2) + ss.b.offset;
    // std::cerr << ref1 << ":" << pos1 + 1 << "-" << pos2 + klen << std::endl;

    // === Get paths in between v1 and v2
    std::map<gbwt::size_type, gbwt::vector_type> paths =
        graph.get_subpaths(v1, v2);
    // std::cerr << "subpath extracted" << std::endl;

    gbwt::vector_type path = paths.begin()->second;
    // - Collapse same path
    // std::map<gbwt::size_type, std::string> collapsed_paths;
    // for (auto &[p, path] : paths) {
    //   int ii = 0;
    //   for (auto &[p2, _] : collapsed_paths) {
    //     if (path.size() != paths[p2].size())
    //       continue;
    //     int j = 0;
    //     for (j = 0; j < path.size(); ++j) {
    //       if (path[j] != paths[p2][j])
    //         break;
    //     }
    //     if (j == path.size())
    //       break;
    //     ++ii;
    //   }
    //   if (ii == collapsed_paths.size()) {
    //     collapsed_paths[p] = std::to_string(p);
    //   } else {
    //     collapsed_paths[ii] += "," + std::to_string(p);
    //   }
    // }

    // - Get string from path
    std::string path_seq = "";
    int tot_plen = 0;
    int lastprefix = -1;
    for (const auto &v : path) {
      std::string seq = graph.get_sequence(v >> 1);
      tot_plen += seq.size();
      if ((v >> 1) == v1) {
        path_seq += seq.substr(ss.a.offset, seq.size());
      } else if ((v >> 1) == v2) {
        path_seq += seq.substr(0, ss.b.offset + klen);
        lastprefix = seq.size() - ss.b.offset - klen;
      } else {
        path_seq += seq;
      }
    }
    assert(lastprefix != -1);
    // - convert to 0123
    int pseq_l = path_seq.size();
    uint8_t *pseq = (uint8_t *)malloc(path_seq.size() + 1);
    strncpy((char *)pseq, path_seq.c_str(), pseq_l);
    pseq[pseq_l] = '\0';
    rb3_char2nt6(pseq_l, pseq);
    for (int i = 0; i < pseq_l; ++i)
      --pseq[i];

    // std::cerr << "path extracted" << std::endl;

    // === Consensus via abpoa
    goods = 0;
    for (sfs_t &s : cluster) {
      cseqs_lens[goods] = s.l;
      cseqs[goods] = s.seq;
      ++goods;
    }
    assert(goods < 64);
    abpoa_msa(ab, abpt, goods, NULL, cseqs_lens, cseqs, NULL, NULL);

    // std::cerr << "abpoa" << std::endl;

    // === Realignment via KSW2
    abpoa_cons_t *abc = ab->abc;
    for (int ci = 0; ci < abc->n_cons; ++ci) {
      int cons_l = abc->cons_len[ci];
      uint8_t *cons_seq = abc->cons_base[ci];
      int abpoa_supp = abc->clu_n_seq[ci];
      // if (abpoa_supp < min_w)
      //   continue;

      ksw_extd2_sse(0, cons_l, cons_seq, pseq_l, pseq, 5, mat, gapo, gape,
                    gapo2, gape2, -1, -1, -1, 0, &ez);

      int clipped = 0;
      if ((ez.cigar[0] & 0xf) != 0 || (ez.cigar[ez.n_cigar - 1] & 0xf) != 0) {
        clipped = 1;
        // continue; // XXX: do we want these?
      }

      // - Parse CIGAR
      int opl;
      int tot_cigar_len = 0;
      int tot_res_matches = 0;

      // XXX: avoid reallocation at each iteration
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
            if (cons_seq[cons_p + j] != pseq[pseq_p + j]) {
              if (l_tmp > 0) {
                cs += ":" + std::to_string(l_tmp);
                l_tmp = 0;
              }
              cs += "*";
              cs += pseq[pseq_p + j] <= 3 ? "ACGT"[pseq[pseq_p + j]] : 'N';
              cs += cons_seq[cons_p + j] <= 3 ? "ACGT"[cons_seq[cons_p + j]]
                                              : 'N';
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
          cs += "+" + decode(cons_seq + cons_p, opl);
          cons_p += opl;
        } else if ((ez.cigar[i] & 0xf) == 2) {
          // D
          cs += "-" + decode(pseq + pseq_p, opl);
          pseq_p += opl;
        } else {
          fprintf(stderr,
                  "Cluster %d --- We shouldn't be here while parsing ksw "
                  "cigar. Halting.\n",
                  cidx);
          exit(EXIT_FAILURE);
        }
      }

      // print GAF line
      std::cout << cidx << "." << ci << "\t";
      std::cout << cons_l << "\t";
      std::cout << 0 << "\t";
      std::cout << cons_l << "\t";
      std::cout << "+\t";
      std::cout << ">" << graph.get_gfa_idx(path[0] >> 1);
      for (int i = 1; i < path.size(); ++i)
        std::cout << ">" << graph.get_gfa_idx(path[i] >> 1).c_str();
      std::cout << "\t";
      std::cout << tot_plen << "\t";
      std::cout << ss.a.offset << "\t";
      std::cout << ss.a.offset + pseq_l << "\t";
      std::cout << tot_res_matches << "\t";
      std::cout << tot_cigar_len << "\t";
      std::cout << 60 << "\t";
      std::cout << "AS:i:" << ez.score << "\t";
      std::cout << "cg:Z:" << cigar << "\t";
      std::cout << "cs:Z:" << cs << "\t";
      std::cout << "cl:Z:" << clipped << "\t";
      std::cout << "cw:Z:" << abpoa_supp << "\t";
      std::cout << "rp:Z:" << ref1 << ":" << pos1 + 1 << "-" << pos2 + klen
                << "\t";
      std::cout << "qs:Z:" << decode(cons_seq, cons_l) << "\t";
      std::cout << "ps:Z:" << decode(pseq, pseq_l);
      std::cout << std::endl;
    }
    // === === ===

    for (sfs_t &s : cluster)
      free(s.seq);
  }

  abpoa_free_para(abpt);
  abpoa_free(ab);
  free(cseqs);
  free(cseqs_lens);

  return 0;
}
