#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <map>
#include <stdio.h>
#include <string>
#include <zlib.h>

extern "C" {
#include "abpoa.h"
#include "fm-index.h"
#include "kseq.h"
}
#include "graph.hpp"
#include "kmer.hpp"
#include "misc.hpp"
#include "sfs.hpp"
#include "sketch.hpp"
#include "usage.h"
#include <iostream>

KSEQ_INIT(gzFile, gzread)

std::string decode(const uint8_t *s, int l) {
  if (s == NULL)
    return "";
  char ds[l + 1];
  for (int i = 0; i < l; ++i)
    ds[i] = "ACGT"[s[i]];
  ds[l] = '\0';
  return ds;
}

typedef struct {
  std::string idx;
  std::string seq;
} read_t;

// Compute SFS strings from P and store them into solutions
std::vector<sfs_t> ping_pong_search(const rb3_fmi_t *index, uint8_t *P, int l,
                                    const std::string qname) {
  std::vector<sfs_t> out;
  rb3_sai_t ik;
  int begin = l - 1;
  while (begin >= 0) {
    // Backward search. Stop at first mismatch.
    int bmatches = 0;
    rb3_fmd_set_intv(index, P[begin], &ik);
    while (ik.size != 0 && begin > 0) {
      --begin;
      ++bmatches;
      rb3_sai_t ok[RB3_ASIZE]; // output SA intervals (one for each symbol
                               // between 0 and 5)
      rb3_fmd_extend(index, &ik, ok, 1);
      ik = ok[P[begin]];
    }
    // last checked char (i.e., first of the query) was a match
    if (begin == 0 && ik.size != 0) {
      break;
    }

    // Forward search
    int end = begin;
    int fmatches = 0;
    rb3_fmd_set_intv(index, P[end], &ik);
    while (ik.size != 0) {
      ++end;
      ++fmatches;
      rb3_sai_t ok[RB3_ASIZE];
      rb3_fmd_extend(index, &ik, ok, 0);
      ik = ok[P[end] >= 1 && P[end] <= 4 ? 5 - P[end] : P[end]];
    }
    sfs_t s;
    s.qname = qname;
    s.s = begin;
    s.l = end - begin + 1;
    out.push_back(s);

    if (begin == 0)
      break;
    // if (config->overlap == 0) // Relaxed
    //   begin -= 1;
    // else
    begin = end - 1;
  }
  return out;
}

// Merge specifics strings that are too close (d-bp apart) on read
void assemble(std::vector<sfs_t> &S, int d) {
  // Reverse the vector
  std::reverse(S.begin(), S.end());

  size_t i = 0;
  while (i < S.size()) {
    size_t j;
    for (j = i + 1; j <= S.size(); ++j) {
      if (j == S.size() || S[j - 1].s + S[j - 1].l <= S[j].s - d) {
        // non-overlapping: update first, clean others
        S[i].l = S[j - 1].s + S[j - 1].l - S[i].s;
        for (size_t j2 = i + 1; j2 < j; ++j2)
          S[j2].l = 0;
        break;
      }
    }
    i = j;
  }

  // Remove gaps by shifting left
  int new_n = 0;
  i = 0;
  while (i < S.size()) {
    if (S[i].l > 0) {
      S[new_n].qname = S[i].qname;
      S[new_n].s = S[i].s;
      S[new_n].l = S[i].l;
      ++new_n;
    }
    ++i;
  }

  // Resise to new size
  S.resize(new_n);
}

bool gcolinear(const sfs_t &s1, const sfs_t &s2) {
  assert(s1.a.v < s1.b.v || (s1.a.v == s1.b.v && s1.a.offset <= s1.b.offset));
  assert(s2.a.v < s2.b.v || (s2.a.v == s2.b.v && s2.a.offset <= s2.b.offset));
  int v11 = s1.a.v, v12 = s1.b.v;
  int v21 = s2.a.v, v22 = s2.b.v;
  int o11 = s1.a.offset, o12 = s1.b.offset;
  int o21 = s2.a.offset, o22 = s2.b.offset;
  bool ret = (v11 < v21 || (v11 == v21 && o11 <= o21)) &&
             (v12 < v22 || (v12 == v22 && o12 <= o22));
  // std::cerr << ret << std::endl;
  return ret;
}

// Merge anchored specifics strings overlapping on same read
void assemble2(std::vector<sfs_t> &S, int strand) {
  // std::cerr << S.size() << std::endl;
  if (S.size() > 0 && strand == 0) {
    // Reverse the vector
    std::reverse(S.begin(), S.end());
  }

  // std::cerr << "Strand: " << strand << std::endl;
  // int ii = -1;
  // for (const auto &s : S) {
  //   ++ii;
  //   // if (s.qidx == -1 || s.strand != strand)
  //   //   continue;
  //   std::cerr << ii << "\t" << s.qidx << " " << (s.strand == strand) << "\t"
  //             << s.s << ":" << s.l << " " << s.a.v << ":" << s.a.offset << "
  //             - "
  //             << s.b.v << ":" << s.b.offset << std::endl;
  // }

  size_t i = 0;
  while (i < S.size()) {
    // skip unanchored or on wrong strand
    if (!S[i].good) {
      ++i;
      continue;
    }

    size_t j;
    size_t last_j = i;
    for (j = i + 1; j <= S.size(); ++j) {
      // skip unanchored or on wrong strand
      if (!S[j].good)
        continue;

      if (j == S.size() || S[last_j].s + S[last_j].l <= S[j].s ||
          !gcolinear(S[last_j], S[j])) {
        // non-overlapping: update first, clean others
        // if (i != last_j)
        //   std::cerr << "XXX" << std::endl;
        S[i].l = S[last_j].s + S[last_j].l - S[i].s;
        // if (S[i].strand == 1)
        S[i].b = S[last_j].b;
        // else
        //   S[i].a = S[last_j].a;
        for (size_t j2 = i + 1; j2 < j; ++j2)
          S[j2].l = 0;
        break;
      }
      last_j = j;
    }
    i = j;
  }

  // Remove gaps by shifting left
  int new_n = 0;
  i = 0;
  while (i < S.size()) {
    if (S[i].l > 0) {
      S[new_n].qname = S[i].qname;
      S[new_n].s = S[i].s;
      S[new_n].l = S[i].l;
      S[new_n].a = S[i].a;
      S[new_n].b = S[i].b;
      ++new_n;
    }
    ++i;
  }

  // Resise to new size
  S.resize(new_n);
  // std::cerr << "--- " << new_n << std::endl;
  // std::cerr << std::endl;
  // for (const auto &s : S) {
  //   if (s.qidx == -1 || s.strand != strand)
  //     continue;
  //   std::cerr << s.s << ":" << s.l << " " << s.a.v << ":" << s.a.offset << "
  //   - "
  //             << s.b.v << ":" << s.b.offset << std::endl;
  // }

  if (S.size() > 0 && strand == 0) {
    // Reverse the vector
    std::reverse(S.begin(), S.end());
  }
}

// Anchor specific strings on graph using graph sketch
void anchor(Graph &graph, sketch_t *sketch, std::vector<sfs_t> &SS,
            const char *Q, int ql, int klen, int NA) {

  int beg, end;
  char *kmer = (char *)malloc(klen);
  uint64_t kmer_d;       // kmer
  uint64_t rckmer_d = 0; // reverse and complemented kmer
  uint64_t ckmer_d = 0;  // canonical kmer
  int c;                 // current char
  hit_t vx;              // hit from sketch

  for (uint sidx = 0; sidx < SS.size(); ++sidx) {
    sfs_t &s = SS[sidx];

    // Finding anchors in flanking regions
    // XXX: Do we want anchors overlapping the string?
    beg = s.s - klen;
    beg = beg < 0 ? 0 : beg;
    end = s.s + s.l;
    end = end > ql - klen ? ql - klen : end;
    memcpy(kmer, Q + beg, klen);

    kmer_d = k2d(kmer, klen);
    rckmer_d = rc(kmer_d, klen);
    ckmer_d = std::min(kmer_d, rckmer_d);

    anchor_t sanchors[NA];
    int nsa = 0;
    while (beg > 0 && nsa < NA) {
      if ((vx = sk_get(sketch, ckmer_d)).first != -1) {
        sanchors[nsa].v = vx.first;
        sanchors[nsa].offset = vx.second;
        sanchors[nsa].p = beg;
        sanchors[nsa].seq = ckmer_d;
        ++nsa;
      }
      --beg;
      c = Q[beg] < 5 ? Q[beg] - 1 : rand() % 4;
      kmer_d = rsprepend(kmer_d, c, klen);
      rckmer_d = lsappend(rckmer_d, reverse_char(c), klen);
      ckmer_d = std::min(kmer_d, rckmer_d);
    }

    memcpy(kmer, Q + end, klen);
    kmer_d = k2d(kmer, klen);
    rckmer_d = rc(kmer_d, klen);
    ckmer_d = std::min(kmer_d, rckmer_d);
    anchor_t eanchors[NA];
    int nea = 0;
    while (end < ql - klen + 1 && nea < NA) {
      if ((vx = sk_get(sketch, ckmer_d)).first != -1) {
        eanchors[nea].v = vx.first;
        eanchors[nea].offset = vx.second;
        eanchors[nea].p = end;
        eanchors[nea].seq = ckmer_d;
        ++nea;
      }
      ++end;
      c = Q[end + klen - 1] < 5 ? Q[end + klen - 1] - 1 : rand() % 4;
      kmer_d = lsappend(kmer_d, c, klen);
      rckmer_d = rsprepend(rckmer_d, reverse_char(c), klen);
      ckmer_d = std::min(kmer_d, rckmer_d);
    }

    if (nsa == 0 || nea == 0) {
      s.good = false; // tag as invalid
      continue;
    }

    // keep only unique anchors on read
    for (int i1 = 0; i1 < nsa; ++i1) {
      // Starting anchors vs starting anchors
      if (sanchors[i1].v == -1)
        continue;
      for (int i2 = i1 + 1; i2 < nsa; ++i2) {
        if (sanchors[i2].v == -1)
          continue;
        if (sanchors[i1].seq == sanchors[i2].seq) {
          sanchors[i1].v = -1;
          sanchors[i2].v = -1;
        }
      }
      // Starting anchors vs ending anchors
      for (int i2 = 0; i2 < nea; ++i2) {
        if (eanchors[i2].v == -1)
          continue;
        if (sanchors[i1].seq == eanchors[i2].seq) {
          sanchors[i1].v = -1;
          eanchors[i2].v = -1;
        }
      }
    }
    // Ending anchors vs ending anchors
    for (int i1 = 0; i1 < nea; ++i1) {
      if (eanchors[i1].v == -1)
        continue;
      for (int i2 = i1 + 1; i2 < nea; ++i2) {
        if (eanchors[i2].v == -1)
          continue;
        if (eanchors[i1].seq == eanchors[i2].seq) {
          eanchors[i1].v = -1;
          eanchors[i2].v = -1;
        }
      }
    }

    // Finding best pair of anchors
    int mind = 100, dist;
    int sax = -1, eax = -1; // index for selected anchors

    std::map<std::pair<int, int>, int>::iterator hhit;
    int x = 0, y = 0;
    int xoff = 0, yoff = 0;
    for (int i = 0; i < nsa; ++i) {
      x = sanchors[i].v;
      if (x == -1)
        // anchor has been filtered out since it was repeated in the read
        continue;
      xoff = sanchors[i].offset;
      for (int j = 0; j < nea; ++j) {
        y = eanchors[j].v;
        if (y == -1)
          // anchor has been filtered out since it was repeated in the read
          continue;
        yoff = eanchors[j].offset;
        dist = graph.distance(graph.get_iidx(sanchors[i].v),
                              graph.get_iidx(eanchors[j].v));
        if (dist < 0)
          continue;
        if (x == y && (xoff == yoff || (xoff < yoff && xoff + klen >= yoff) ||
                       (xoff > yoff && yoff + klen >= xoff)))
          continue;
        if (dist < mind) {
          sax = i;
          eax = j;
          mind = dist;
        }
      }
    }
    if (sax == -1 || eax == -1) {
      s.good = false; // tag as invalid
      continue;
    }

    // Assigning the anchors
    anchor_t sa = sanchors[sax];
    anchor_t ea = eanchors[eax];
    int b = sa.p;
    int l = ea.p + klen - sa.p;
    int strand = 1;
    if (sa.v > ea.v || (sa.v == ea.v && sa.offset > ea.offset)) {
      anchor_t tmp = sa;
      sa = ea;
      ea = tmp;
      strand = 0;
    }

    s.s = b;
    s.l = l;
    s.a = sa;
    s.a.closest = sax;
    s.b = ea;
    s.b.closest = eax;
    s.strand = strand;

    assert(sa.v <= ea.v);
  }
  free(kmer);
}

int load_batch(kseq_t *seq, std::vector<read_t> &entries, int nb) {
  int l = 0;
  int n = 0;
  while (n < nb && (l = kseq_read(seq)) >= 0) {
    entries[n++] = {seq->name.s, seq->seq.s};
  }
  return n;
}

int main_search(int argc, char *argv[]) {
  double rt, rt1;

  int klen = 27;
  int d = 0;         // merge specific strings this close
  int bsize = 10000; // batch size
  // int hd = 0;  // hamming distance for fixing anchors
  int NA = 20;    // number of kmers to check for anchoring
  int nth = 4;    // number of threads
  uint min_w = 2; // minimum size for clusters
  int _c;
  while ((_c = getopt(argc, argv, "k:w:a:b:@:h")) != -1) {
    switch (_c) {
    case 'k':
      klen = std::stoi(optarg);
      break;
    case 'w':
      min_w = std::stoi(optarg);
      break;
    case 'b':
      bsize = std::stoi(optarg);
      break;
    case 'a':
      NA = std::stoi(optarg);
      break;
    case '@':
      nth = std::stoi(optarg);
      break;
    case 'h':
      fprintf(stderr, "%s", SEARCH_USAGE_MESSAGE);
      return 0;
    }
  }
  if (argc - optind != 4) {
    fprintf(stderr, "%s", SEARCH_USAGE_MESSAGE);
    return 1;
  }
  std::string gfa_fn = argv[optind++];
  std::string skt_fn = argv[optind++];
  std::string fmd_fn = argv[optind++];
  std::string fq_fn = argv[optind++];

  // Graph and sketch loading
  rt = realtime();
  sketch_t *sketch = sk_load(skt_fn.c_str());
  fprintf(stderr, "[M::%s] loaded %ld sketches in %.3f sec\n", __func__,
          sketch->n, realtime() - rt);

  rt = realtime();
  Graph graph(gfa_fn);
  graph.load_vertices();
  graph.load_edges();
  graph.load_paths("CHM13"); // XXX: hardcoded
  fprintf(stderr,
          "[M::%s] loaded %ld vertices, %d edges, and %ld paths in %.3f sec\n",
          __func__, graph.vertices.size(), graph.ne, graph.paths.size(),
          realtime() - rt);

  // rt = realtime();
  // graph.build_distance_index();
  // fprintf(stderr, "[M::%s] build distance index in %.3f sec\n", __func__,
  //         realtime() - rt);

  // FMD-index loading
  rt = realtime();
  rb3_fmi_t fmd;
  rb3_fmi_restore(&fmd, fmd_fn.c_str(), 0);
  if (fmd.e == 0 && fmd.r == 0) {
    fprintf(stderr, "[E::%s] FMD-index cannot be restored", __func__);
    exit(EXIT_FAILURE);
  }
  fprintf(stderr, "[M::%s] restored FMD index in %.3f sec\n", __func__,
          realtime() - rt);
  rt = realtime();

  // Specific strings computation and anchoring

  std::vector<read_t> entries(bsize);
  std::vector<std::vector<sfs_t>> output(bsize);
  fprintf(stderr, "[M::%s] pre-allocation in %.3f sec\n", __func__,
          realtime() - rt);
  rt = realtime();

  gzFile fp = gzopen(fq_fn.c_str(), "r");
  kseq_t *seq = kseq_init(fp);

  // Some statistics
  std::vector<sfs_t> specific_strings;

  int anchored_n = 0;
  int unanchored_n = 0;
  int assembled_n = 0;

  int x;
  int analyzed = 0;
  rt1 = realtime();
  while ((x = load_batch(seq, entries, bsize)) > 0) {
    fprintf(stderr, "[M::%s] loaded %d reads in %.3f sec\n", __func__, x,
            realtime() - rt1);
    rt1 = realtime();
#pragma omp parallel for num_threads(nth) schedule(static, 1)
    for (int qq = 0; qq < x; ++qq) {
      uint8_t *seq = (uint8_t *)entries[qq].seq.data();
      int l = entries[qq].seq.size();
      rb3_char2nt6(l, seq);

      output[qq] = ping_pong_search(&fmd, seq, l, entries[qq].idx);
      assemble(output[qq], d);
      anchor(graph, sketch, output[qq], entries[qq].seq.c_str(), l, klen, NA);

      // strand stuff
      int strands[2] = {0, 0};
      for (uint ss = 0; ss < output[qq].size(); ++ss) {
        if (output[qq][ss].good)
          ++strands[output[qq][ss].strand];
      }
      // + strand if tie
      int strand = strands[0] > strands[1] ? 0 : 1;
      for (uint ss = 0; ss < output[qq].size(); ++ss) {
        if (!output[qq][ss].good)
          continue;
        if (output[qq][ss].strand == 0)
          // reverse
          output[qq][ss].s = l - (output[qq][ss].s + output[qq][ss].l);
        if (output[qq][ss].strand != strand)
          output[qq][ss].good = false;
      }

      // Finally reassemble good ones
      assemble2(output[qq], strand);

      if (strand == 0) {
        // Reverse and complement the sequence
        int i;
        for (i = 0; i < (l >> 1); ++i) {
          int tmp = seq[l - 1 - i];
          tmp = (tmp >= 1 && tmp <= 4) ? 5 - tmp : tmp;
          seq[l - 1 - i] = (seq[i] >= 1 && seq[i] <= 4) ? 5 - seq[i] : seq[i];
          seq[i] = tmp;
        }
        if (l & 1)
          seq[i] = (seq[i] >= 1 && seq[i] <= 4) ? 5 - seq[i] : seq[i];
      }

      // Fill sequence
      for (sfs_t &s : output[qq]) {
        s.seq = (uint8_t *)malloc(s.l + 1);
        memcpy(s.seq, seq + s.s, s.l);
        s.seq[s.l] = '\0';
        // -1 since abpoa works on 0123
        for (int i = 0; i < s.l; ++i)
          --s.seq[i];
      }
    }
    fprintf(stderr, "[M::%s] searched in %.3f sec\n", __func__,
            realtime() - rt1);
    rt1 = realtime();

    // output
    for (int qq = 0; qq < x; ++qq) {
      assembled_n += output[qq].size();
      for (uint j = 0; j < output[qq].size(); ++j) {
        sfs_t s = output[qq][j];
        if (!s.good) {
          ++unanchored_n;
          // printf("X %d %s %d %d . . .:.:.:. .:.:.:.\n", qq,
          //        entries[qq].idx.c_str(), s.s, s.l);
        } else {
          ++anchored_n;

          specific_strings.push_back(s);
          // printf("O %s %d %d %d %s %ld:%d:%ld %ld:%d:%ld %d %d\n",
          //        entries[qq].idx.c_str(), s.s, s.l, s.strand, ".", s.a.v,
          //        s.a.offset, s.a.seq, s.b.v, s.b.offset, s.b.seq,
          //        s.a.closest, s.b.closest);
        }
      }
    }
    analyzed += x;
  }

  // At this point, specific strings are anchored. Anchors follow + strand
  // on graph. If read was on -, we have reversed the specific strings so
  // that everything is on + strand

  fprintf(
      stderr,
      "[M::%s] Computed %d specific strings (%.3f%% anchored, %d unanchored) "
      "in %.3f sec\n",
      __func__, assembled_n, anchored_n / (float)assembled_n * 100,
      unanchored_n, realtime() - rt);

  kseq_destroy(seq);
  gzclose(fp);
  sk_destroy(sketch);
  rb3_fmi_free(&fmd);

  // Cluster specific strings based on their anchors
  rt = realtime();
  std::sort(specific_strings.begin(), specific_strings.end(),
            [](const sfs_t &a, const sfs_t &b) { return a.a.seq < b.a.seq; });

  std::vector<std::vector<sfs_t>> clusters;
  clusters.push_back({specific_strings[0]});
  uint last_c = 0;
  for (uint ss = 1; ss < specific_strings.size(); ++ss) {
    sfs_t sfs = specific_strings[ss];
    // printf("%s\n", decode(sfs.seq, sfs.l).c_str());
    if (sfs.a.seq == clusters[last_c][0].a.seq) {
      uint c = last_c;
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

  fprintf(stderr, "[M::%s] Built %ld clusters in %.3f sec\n", __func__,
          clusters.size(), realtime() - rt);

  // === Consensus via abpoa
  rt = realtime();

  // Init abpoa
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

  int cidx = 0;
  rt = realtime();
  // TODO: parallelize
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

    sfs_t ss = cluster[0];

    // std::ofstream ofs(std::to_string(cidx) + ".txt", std::ofstream::out);
    // === Consensus via abpoa
    goods = 0;
    for (sfs_t &s : cluster) {
      // ofs << s.qname << ":" + s.s << "-" << s.s + s.l << "/" << s.l << "/"
      //     << s.strand << "/" << s.a.v << ">" << s.b.v << "/" << s.a.offset
      //     << ":" << s.b.offset << "/" << s.a.seq << ":" << s.b.seq << "\n";
      // ofs << decode(s.seq, s.l) << "\n";
      cseqs_lens[goods] = s.l;
      cseqs[goods] = s.seq;
      ++goods;
    }
    // ofs.close();
    assert(goods < 64);
    abpoa_msa(ab, abpt, goods, NULL, cseqs_lens, cseqs, NULL, NULL);

    abpoa_cons_t *abc = ab->abc;
    for (int ci = 0; ci < abc->n_cons; ++ci) {
      int cons_l = abc->cons_len[ci];
      // char *cons_seq = (char *)abc->cons_base[ci];
      uint8_t *cons_seq = abc->cons_base[ci];
      for (int i = 0; i < cons_l; ++i)
        cons_seq[i] = "ACGT"[cons_seq[i]];
      cons_seq[cons_l] = '\0';
      int abpoa_supp = abc->clu_n_seq[ci];
      // if (abpoa_supp < min_w)
      //   continue;
      std::cout << ">" << cidx << "." << ci << "." << cons_l << "."
                << abpoa_supp << "." << ss.a.v << ":" << ss.a.offset << "."
                << ss.b.v << ":" << ss.b.offset << std::endl;
      std::cout << cons_seq << std::endl;
    }
    // === === ===

    for (sfs_t &s : cluster)
      free(s.seq);
  }

  abpoa_free_para(abpt);
  abpoa_free(ab);
  free(cseqs);
  free(cseqs_lens);

  fprintf(stderr, "[M::%s] Built consensus in %.3f sec\n", __func__,
          realtime() - rt);

  return 0;
}
