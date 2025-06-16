#include <fstream>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

extern "C" {
#include "ksw2.h"
}

#include "graph.hpp"
#include "usage.h"

static const uint8_t to_int[128] = {0, 0, 1, 2, 3, 0, 0, 0, 0, 0, // 0
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 10
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 20
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 30
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 40
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 50
                                    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, // 60
                                    0, 2, 0, 0, 0, 0, 0, 0, 0, 0, // 70
                                    0, 0, 0, 0, 3, 0, 0, 0, 0, 0, // 80
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, // 90
                                    0, 0, 0, 2, 0, 0, 0, 0, 0, 0, // 100
                                    0, 0, 0, 0, 0, 0, 3, 0, 0, 0, // 110
                                    0, 0, 0, 0, 0, 0, 0, 0};      // 120

int main_realign(int argc, char *argv[]) {
  // TODO: we should directly take GraphAligner .gaf (no need for python script)

  // double rt = realtime(), rt1;

  int klen = 27;  // kmer size
  int min_as = 0; // minimumm alignment score

  int _c;
  while ((_c = getopt(argc, argv, "k:s:h")) != -1) {
    switch (_c) {
    case 'k':
      klen = std::stoi(optarg);
      break;
    case 's':
      min_as = std::stoi(optarg);
      break;
    case 'h':
      fprintf(stderr, "%s", "CALL_USAGE_MESSAGE");
      return 0;
    default:
      fprintf(stderr, "Error\n");
      return 1;
    }
  }
  if (argc - optind != 2) {
    fprintf(stderr, "Error!\n");
    return 1;
  }

  std::string gfa_fn = argv[optind++];
  std::string txt_fn = argv[optind++];

  // Graph
  Graph graph(gfa_fn);
  graph.load_vertices();

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

  std::ifstream fp;
  fp.open(txt_fn);
  for (std::string line; std::getline(fp, line);) {
    std::istringstream iss(line);
    std::string idx, path_s, consensus;
    iss >> idx;
    iss >> path_s;
    iss >> consensus;

    int cidx, subcidx, /*length,*/ support;
    std::string a1, a2;
    int /*v1,*/ off1, v2, off2;

    std::stringstream idx_ss(idx);
    std::string token;

    getline(idx_ss, token, '.');
    cidx = std::stoi(token);

    getline(idx_ss, token, '.');
    subcidx = std::stoi(token);

    getline(idx_ss, token, '.');
    // _length = std::stoi(token);

    getline(idx_ss, token, '.');
    support = std::stoi(token);

    getline(idx_ss, a1, '.');
    std::stringstream a1_ss(a1);
    getline(a1_ss, token, ':');
    // v1 = std::stoi(token);
    getline(a1_ss, token, ':');
    off1 = std::stoi(token);

    getline(idx_ss, a2, '.');
    std::stringstream a2_ss(a2);
    getline(a2_ss, token, ':');
    v2 = std::stoi(token);
    getline(a2_ss, token, ':');
    off2 = std::stoi(token);

    // XXX: pseq can be built directly without substr
    std::vector<int> path;
    std::string full_pseq;
    std::stringstream path_ss(path_s);
    int ps = off1, pe = 0;
    while (getline(path_ss, token, ',')) {
      path.push_back(std::stoi(token));
      if (path.back() == v2)
        pe = full_pseq.size() + off2 + klen;
      full_pseq += graph.get_sequence(graph.get_iidx(path.back()));
    }
    std::string pseq = full_pseq.substr(ps, pe - ps);

    uint8_t *pseq_c = (uint8_t *)malloc(pseq.size() + 1);
    for (uint i = 0; i < pseq.size(); ++i)
      pseq_c[i] = to_int[pseq[i]];
    pseq_c[pseq.size()] = '\0';

    uint8_t *consensus_c = (uint8_t *)malloc(consensus.size() + 1);
    for (uint i = 0; i < consensus.size(); ++i)
      consensus_c[i] = to_int[consensus[i]];
    consensus_c[consensus.size()] = '\0';

    ksw_extd2_sse(0, consensus.size(), consensus_c, pseq.size(), pseq_c, 5, mat,
                  gapo, gape, gapo2, gape2, -1, -1, -1, 0, &ez);

    // int clipped = 0;
    if ((ez.cigar[0] & 0xf) != 0 || (ez.cigar[ez.n_cigar - 1] & 0xf) != 0) {
      // clipped = 1;
      continue; // XXX: do we want these?
    }

    if (ez.score < min_as)
      continue;

    // Parse CIGAR
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
          if (consensus[cons_p + j] != pseq[pseq_p + j]) {
            if (l_tmp > 0) {
              cs += ":" + std::to_string(l_tmp);
              l_tmp = 0;
            }
            cs += "*";
            cs += pseq[pseq_p + j];
            cs += consensus[cons_p + j];
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
        cs += "+" + consensus.substr(cons_p, opl);
        cons_p += opl;
      } else if ((ez.cigar[i] & 0xf) == 2) {
        // D
        cs += "-" + pseq.substr(pseq_p, opl);
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
    std::cout << cidx << "." << subcidx << "\t";
    std::cout << consensus.size() << "\t";
    std::cout << 0 << "\t";
    std::cout << consensus.size() << "\t";
    std::cout << "+\t";
    std::cout << ">" << path[0];
    for (uint i = 1; i < path.size(); ++i)
      std::cout << ">" << path[i];
    std::cout << "\t";
    std::cout << full_pseq.size() << "\t";
    std::cout << off1 << "\t";
    std::cout << off1 + pseq.size() << "\t";
    std::cout << tot_res_matches << "\t";
    std::cout << tot_cigar_len << "\t";
    std::cout << 60 << "\t";
    std::cout << "AS:i:" << ez.score << "\t";
    std::cout << "cg:Z:" << cigar << "\t";
    std::cout << "cs:Z:" << cs << "\t";
    // std::cout << "cl:Z:" << clipped << "\t";
    std::cout << "cw:Z:" << support << "\t";
    // std::cout << "rp:Z:" << ref1 << ":" << pos1 + 1 << "-" << pos2 + klen
    //           << "\t";
    std::cout << "qs:Z:" << consensus << "\t";
    std::cout << "ps:Z:" << pseq;
    std::cout << std::endl;

    free(pseq_c);
    free(consensus_c);
  }
  fp.close();

  return 0;
}
