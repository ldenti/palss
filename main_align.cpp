// #include <algorithm>
#include <getopt.h>
#include <iostream>
#include <omp.h>
#include <string>
#include <vector>

extern "C" {
#include "abpoa.h"
#include "ksw2.h"
}

#include "gaf.hpp"
#include "graph.hpp"
#include "kmer.hpp"
#include "misc.hpp"
#include "sfs.hpp"
#include "usage.hpp"

typedef struct {
  std::vector<gbwt::size_type> vertices; // with strand
  std::string sequence;
  std::vector<uint32_t> kcounts;
  //
  uint32_t total_length;
  uint32_t skipped_prefix;
  uint32_t skipped_suffix;

} ph_t;
// bool operator==(const path_t &path1, const path_t &path2);

// TODO: clusters are pair of integers defining intervals in specific_strings
std::vector<std::vector<sfs_t>>
cluster_sfs(const std::vector<sfs_t> &specific_strings) {
  std::vector<std::vector<sfs_t>> clusters;
  // XXX: we are sure to have at least one specific string here
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
  return clusters;
}

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

// 0123-sequence to 0123-kmer counts
std::vector<uint32_t> count_kmers(const char *seq, int seql, int klen) {
  std::vector<uint32_t> kcounts((1 << (2 * klen)), 0);
  int c; // current char

  uint64_t kmer_d = 0;
  for (uint8_t i = 0; i < klen; ++i) {
    kmer_d = (kmer_d << 2) | (seq[i] < 4 ? seq[i] : rand() % 4);
  }
  ++kcounts[kmer_d];
  for (int p = klen; p < seql; ++p) {
    c = seq[p];
    kmer_d = lsappend(kmer_d, c, klen);
    ++kcounts[kmer_d];
  }
  return kcounts;
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
    // ACGT -> 0123
    c = to_int[(uint8_t)seq[p]] - 1;
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

  int cklen = 5;
  int nth = 4;
  size_t max_plen = 100000;

  int _c;
  while ((_c = getopt(argc, argv, "@:m:h")) != -1) {
    switch (_c) {
    case '@':
      nth = std::stoi(optarg);
      break;
    case 'm':
      max_plen = std::stoi(optarg);
      break;
    case 'h':
      fprintf(stderr, "%s", ALIGN_USAGE_MESSAGE);
      return 0;
    }
  }
  if (argc - optind != 2) {
    fprintf(stderr, "%s", ALIGN_USAGE_MESSAGE);
    return 1;
  }

  std::string gbz_fn = argv[optind++];
  std::string sfs_fn = argv[optind++];

  /** Load graph ****************************************************/
  Graph graph(gbz_fn, "CHM13");
  graph.load();
  graph.load_fl();

  const gbwtgraph::GBZ &gbz = graph.get_gbz();
  const gbwt::FastLocate &fl = graph.get_fl();

  /** Load specific strings *****************************************/
  rt = realtime();
  std::vector<sfs_t> specific_strings = load_sfs(sfs_fn);
  if (specific_strings.empty()) {
    fprintf(stderr, "[W::%s] No specific strings. Halting..\n", __func__);
    return 1;
  }
  fprintf(stderr, "[M::%s] Loaded %ld specific strings in %.3f sec\n", __func__,
          specific_strings.size(), realtime() - rt);

  /** Cluster specific strings based on their anchors ***************/
  rt = realtime();
  std::sort(specific_strings.begin(), specific_strings.end(),
            [](const sfs_t &a, const sfs_t &b) { return a.skmer < b.skmer; });
  fprintf(stderr, "[M::%s] Sorted specific strings in %.3f sec\n", __func__,
          realtime() - rt);

  rt = realtime();
  std::vector<std::vector<sfs_t>> clusters = cluster_sfs(specific_strings);
  fprintf(stderr, "[M::%s] Built %ld clusters in %.3f sec\n", __func__,
          clusters.size(), realtime() - rt);

  // Initialize what we need for main loop
  std::vector<abpoa_t *> abpoa_ts(nth);
  std::vector<abpoa_para_t *> abpoa_para_ts(nth);
  std::vector<ksw_extz_t> ksws(nth);

  // ksw2 parameters
  // https://github.com/lh3/minimap2/blob/69e36299168d739dded1c5549f662793af10da83/options.c#L144
  // asm5
  // int sc_mch = 1, sc_mis = -19, gapo = 39, gape = 3, gapo2 = 81, gape2 = 1;
  // asm10
  int sc_mch = 1, sc_mis = -9, gapo = 16, gape = 2, gapo2 = 41, gape2 = 1;
  int8_t a = (int8_t)sc_mch,
         b = sc_mis < 0 ? (int8_t)sc_mis : -(int8_t)sc_mis; // a>0 and b<0
  int8_t mat[25] = {a, b, b, b, 0, b, a, b, b, 0, b, b, a,
                    b, 0, b, b, b, a, 0, 0, 0, 0, 0, 0};

  std::vector<uint8_t **> cluster_seqs(nth);
  std::vector<int *> cluster_seqs_lens(nth);

  for (int t = 0; t < nth; ++t) {
    // init abpoa
    abpoa_ts[t] = abpoa_init();
    abpoa_para_ts[t] = abpoa_init_para();
    abpoa_para_t *abpt = abpoa_para_ts[t];
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

    // init ksw2
    memset(&ksws[t], 0, sizeof(ksw_extz_t));

    // Array to store specific strings sequences
    // XXX: Assuming clusters of size <= 64
    // TODO: reallocation or we could create a struct cluster_t with the array
    cluster_seqs[t] = (uint8_t **)malloc(sizeof(uint8_t *) * 64);
    cluster_seqs_lens[t] = (int *)malloc(sizeof(int) * 64);
  }

  rt = realtime();
#pragma omp parallel for num_threads(nth) schedule(static, 1)
  for (size_t cidx = 0; cidx < clusters.size(); ++cidx) {
    int tt = omp_get_thread_num();

    std::vector<sfs_t> &cluster = clusters[cidx];

    // if (cidx > 1000)
    //   continue;
    // if (cidx != 30337)
    //   continue;

    // std::cerr << "=== " << cidx << " ===" << std::endl;

    //   if ((cidx + 1) % 10000 == 0) {
    //     fprintf(stderr, "[M::%s] analyzed %ld/%ld clusters in %.3f sec\n",
    //             __func__, cidx + 1, clusters.size(), realtime() - rt);
    //   }

    if (cluster.size() > 64) {
      // FIXME
      fprintf(stderr, "[W::%s] skipping cluster %ld since >64\n", __func__,
              cidx);
      continue;
    }

    const sfs_t &s0 = cluster.front();
    // std::cerr << s0.sv << " " << s0.ev << std::endl;
    // std::cerr << s0.skmer << " " << s0.ekmer << std::endl;
    // std::cerr << graph.get_gfa_name(s0.sv >> 1) << ("+-"[s0.ev & 1]) << " "
    //           << graph.get_gfa_name(s0.ev >> 1) << ("+-"[s0.ev & 1])
    //           << std::endl;

    /** Extract paths from first specific string *********************/
    // XXX: is it enough to consider paths from s0 only?
    std::vector<ph_t> paths;
    for (const uint64_t &p : s0.paths) {

      // Get anchor information from first specific string
      uint32_t sv = s0.sv;
      uint32_t soff = s0.soff;
      if ((p >> 1) & 1) {
        sv = s0.sv ^ 1;
        soff = graph.get_vertex_len(sv >> 1) - s0.soff - 1;
      }

      uint32_t ev = s0.ev;
      uint32_t eoff = s0.eoff;
      if (p & 1) {
        ev = s0.ev ^ 1;
        eoff = graph.get_vertex_len(ev >> 1) - s0.eoff - 1;
      }

      gbwt::size_type seqid = p >> 4;

      std::vector<gbwt::size_type> intervals1 = fl.decompressSA(sv);
      std::vector<gbwt::size_type> intervals2 = fl.decompressSA(ev);

      for (size_t i = 0; i < intervals1.size(); ++i) {
        gbwt::FastLocate::size_type int1 = intervals1[i];
        gbwt::size_type seqid1 = fl.seqId(int1);
        if (seqid1 != seqid)
          continue;

        gbwt::size_type seqoff1 = fl.seqOffset(int1);
        for (size_t j = 0; j < intervals2.size(); ++j) {
          gbwt::FastLocate::size_type int2 = intervals2[j];
          gbwt::size_type seqid2 = fl.seqId(int2);
          if (seqid1 != seqid2)
            continue;

          gbwt::size_type seqoff2 = fl.seqOffset(int2);

          gbwt::edge_type position = std::make_pair(sv, i);
          gbwt::edge_type position2 = std::make_pair(ev, j);

          if (sv == ev && seqoff2 != seqoff1) {
            // if we have cycles (e.g., same vertex more times along the same
            // path), we might end up with the same vertex but different offset
            // along the path and this will break everything
            continue;
          }
          // XXX: what about cycles when sv != ev?
          // We should be able to consider only good pairs using bits in the
          // path identifier

          uint32_t local_sv = sv, local_ev = ev;
          uint32_t local_soff = soff, local_eoff = eoff;
          if (seqoff2 > seqoff1) {
            std::swap(position, position2);
            std::swap(local_sv, local_ev);
            std::swap(local_soff, local_eoff);
          } else if (seqoff2 == seqoff1 && local_soff > local_eoff) {
            std::swap(local_soff, local_eoff);
          }

          assert(fl.index->contains(position));
          assert(fl.index->contains(position2));

          // std::cerr << gbz.index.metadata.fullPath(seqid1 >> 1).sample_name
          //           << " "
          //           << gbz.index.metadata.fullPath(seqid1 >> 1).contig_name
          //           << " " << graph.get_gfa_name(position.first >> 1)
          //           << "+-"[position.first & 1] << ":" << local_soff << " "
          //           << graph.get_gfa_name(position2.first >> 1)
          //           << "+-"[position2.first & 1] << ":" << local_eoff
          //           << std::endl;

          ph_t path;
          path.total_length = 0;
          while (position.first != position2.first) {
            // assert(position.first != 0);

            // position.first is encoded with strand
            path.vertices.push_back(position.first);
            gbwtgraph::handle_t handle =
                gbwtgraph::GBWTGraph::node_to_handle(position.first);
            // in-place view of the sequence: (start, length)
            std::string_view view = gbz.graph.get_sequence_view(handle);

            path.total_length += view.size();

            if (position.first == local_sv) {
              // cut prefix
              path.skipped_prefix = local_soff;
              path.sequence.append(view.data() + local_soff,
                                   view.size() - local_soff);
            } else
              path.sequence.append(view.data(), view.size());
            position = gbz.index.LF(position);
          }

          assert(position.first == local_ev);
          path.vertices.push_back(position.first);
          gbwtgraph::handle_t handle =
              gbwtgraph::GBWTGraph::node_to_handle(position.first);
          std::string_view view = gbz.graph.get_sequence_view(handle);
          path.total_length += view.size();
          if (position.first == local_sv) {
            path.skipped_prefix = local_soff;
            path.sequence.append(view.data() + local_soff,
                                 local_eoff + 1 - local_soff);
            path.skipped_suffix = view.size() - local_eoff - 1;
          } else {
            path.sequence.append(view.data(), local_eoff + 1);
            path.skipped_suffix = view.size() - local_eoff - 1;
          }

          // std::cerr << gbwt::Path::is_reverse(seqid) << " " << path.sequence
          //           << std::endl;

          bool has_n = false;
          for (size_t i = 0; i < path.sequence.size(); ++i) {
            if (path.sequence[i] == 'N') {
              has_n = true;
              break;
            }
            // ACGT -> 0123
            path.sequence[i] = to_int[(uint8_t)path.sequence[i]] - 1;
          }
          if (has_n)
            continue;

          path.kcounts =
              count_kmers(path.sequence.c_str(), path.sequence.size(), cklen);

          if (path.sequence.size() <= max_plen)
            // XXX: do we want this?
            paths.push_back(path);
        }
      }
    }

    if (paths.empty())
      continue;

    // XXX: path sequences should all start/end in the same way (same "strand")

    /** **************************************************************/

    // Consensus via abpoa
    int goods = 0;
    for (sfs_t &s : cluster) {
      cluster_seqs_lens[tt][goods] = s.l;
      cluster_seqs[tt][goods] = s.seq;
      ++goods;
    }
    abpoa_t *ab = abpoa_ts[tt];
    abpoa_msa(ab, abpoa_para_ts[tt], goods, NULL, cluster_seqs_lens[tt],
              cluster_seqs[tt], NULL, NULL);

    abpoa_cons_t *abc = ab->abc;
    for (int sub_cidx = 0; sub_cidx < abc->n_cons; ++sub_cidx) {
      int cons_l = abc->cons_len[sub_cidx];
      // int abpoa_supp = abc->clu_n_seq[sub_cidx];
      char *cons_seq = (char *)abc->cons_base[sub_cidx];

      // XXX: seems we need ACGT just to count kmers from plain... we might
      // avoid it and keep the string 0123-encoded

      // === select best path ===
      std::vector<uint32_t> consensus_kcounts =
          count_kmers(cons_seq, cons_l, cklen);
      size_t best_p = -1U;
      float cd = 0, best_cd = INT64_MAX;
      for (size_t i = 0; i < paths.size(); ++i) {
        const auto &pkc = paths[i].kcounts;
        if (pkc.empty())
          continue;
        cd = canberra(consensus_kcounts, pkc);
        if (cd < best_cd) {
          best_cd = cd;
          best_p = i;
        }
      }
      assert(best_p != -1U);

      ph_t &path = paths[best_p];

      //   if (path.reversed)
      //     cseq_plain = reverseAndComplement(cseq_plain);

      // XXX: path longer are actually filtered out when computing them above
      // if (path.sequence.size() > max_plen)
      //   continue;

      ksw_extz_t ez = ksws[tt];
      ksw_extd2_sse(0, cons_l, (uint8_t *)cons_seq, path.sequence.size(),
                    (uint8_t *)path.sequence.c_str(), 5, mat, gapo, gape, gapo2,
                    gape2, -1, -1, -1, 0, &ez);

      /** OUTPUT *****************************************************/

      // go back to ACGT
      for (int i = 0; i < cons_l; ++i)
        cons_seq[i] = "ACGT"[(int)cons_seq[i]];
      cons_seq[cons_l] = '\0';

      // std::cerr << "\n" << cons_seq << std::endl;

      std::string pseq(path.sequence);
      for (size_t i = 0; i < path.sequence.size(); ++i)
        pseq[i] = "ACGT"[(int)pseq[i]];

      // bool clipped = false;
      if ((ez.cigar[0] & 0xf) != 0 || (ez.cigar[ez.n_cigar - 1] & 0xf) != 0) {
        // clipped = true;
        // free(path_seq);
        continue; // XXX: do we want these?
      }

      // std::cerr << path.sequence << std::endl;
      // std::cerr << cons_seq << std::endl;

      // Build CIGAR and difference string
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
            if (cons_seq[cons_p + j] != pseq[pseq_p + j]) {
              if (l_tmp > 0) {
                cs += ":" + std::to_string(l_tmp);
                l_tmp = 0;
              }
              cs += "*";
              cs += pseq[pseq_p + j];
              cs += cons_seq[cons_p + j];
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
          cs += "+" + std::string(cons_seq + cons_p, opl);
          cons_p += opl;
        } else if ((ez.cigar[i] & 0xf) == 2) {
          // D
          cs += "-" + pseq.substr(pseq_p, opl);
          pseq_p += opl;
        } else {
          std::cerr << "Cluster " << cidx << ": error in ksw2 CIGAR"
                    << std::endl;
          exit(EXIT_FAILURE);
        }
      }

      // Get offsets along path
      // since vertices might be split (if longer than 1024), the offset along
      // path refers to "internal" vertices, not GFA segments
      gbwtgraph::handle_t handle =
          gbwtgraph::GBWTGraph::node_to_handle(path.vertices[0]);
      size_t xx = graph.get_gbz().graph.get_segment_offset(handle);
      size_t ps = xx + path.skipped_prefix;

      // Build and write GAF record
      GAFREC gr;
      gr.qname =
          "palss-" + std::to_string(cidx) + "." + std::to_string(sub_cidx);
      gr.qlen = cons_l;
      gr.qs = 0;
      gr.qe = cons_l;
      gr.strand = true; // TODO
      std::string segment_name = graph.get_gfa_name(path.vertices[0] >> 1);
      gr.path.push_back(((path.vertices[0] & 1) ? "<" : ">") + segment_name);
      std::string last_segment_name = segment_name;
      for (size_t x = 1; x < path.vertices.size(); ++x) {
        segment_name = graph.get_gfa_name(path.vertices[x] >> 1);
        if (segment_name.compare(last_segment_name) != 0) {
          gr.path.push_back(((path.vertices[x] & 1) ? "<" : ">") +
                            segment_name);
          last_segment_name = segment_name;
        }
      }
      gr.plen = path.total_length;
      gr.ps = ps;
      gr.pe = ps + pseq.size();
      gr.tot_res_matches = tot_res_matches;
      gr.tot_cigar_len = tot_cigar_len;
      gr.mapq = 60;
      gr.as = ez.score;
      gr.cigar = cigar;
      gr.cs = cs;
      // gr.clipped = clipped;
      for (int x = 0; x < abc->clu_n_seq[sub_cidx]; ++x)
        gr.reads.push_back(cluster[abc->clu_read_ids[sub_cidx][x]].rname);
      gr.qseq = cons_seq;
      gr.pseq = pseq;

      gr.write();
    }
  }

  for (int t = 0; t < nth; ++t) {
    free(ksws[t].cigar);
    abpoa_free_para(abpoa_para_ts[t]);
    abpoa_free(abpoa_ts[t]);
    free(cluster_seqs[t]);
    free(cluster_seqs_lens[t]);
  }

  for (sfs_t &s : specific_strings)
    free(s.seq);

  return 0;
}
