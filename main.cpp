#include <algorithm>
#include <cstdint>
#include <iostream>
#include <vector>
#include <zlib.h>

#include "fm-index.h"
// #include "kmc_api/kmc_file.h"
#include "kseq.h"

#include "gsketch.hpp"

// KSEQ_INIT(gzFile, gzread) // we already init kstream in gsketch
__KSEQ_TYPE(gzFile)
__KSEQ_BASIC(static, gzFile)
__KSEQ_READ(static)

using namespace std;

int main_call(int argc, char *argv[]);

struct sfs_t {
  int s;       // start on query
  int l;       // length
  int sn = -1; // starting node with "unique" kmer
  int en = -1; // ending node with "unique" kmer
  int strand = 1;
};

vector<sfs_t> assemble_lr(const vector<sfs_t> &sfs) {
  vector<sfs_t> assembled_sfs;
  int i = 0;
  while (i < sfs.size()) {
    int j;
    for (j = i + 1; j < sfs.size(); ++j) {
      if (sfs[j - 1].s + sfs[j - 1].l <= sfs[j].s - 200) { // FIXME hardcoded
        // non-overlapping
        int l = sfs[j - 1].s + sfs[j - 1].l - sfs[i].s;
        assembled_sfs.push_back({sfs[i].s, l});
        i = j;
        break;
      }
    }
    if (j == sfs.size()) {
      int l = sfs[j - 1].s + sfs[j - 1].l - sfs[i].s;
      assembled_sfs.push_back({sfs[i].s, l});
      i = j;
    }
  }
  return assembled_sfs;
}

vector<sfs_t> assemble_rl(const vector<sfs_t> &sfs) {
  vector<sfs_t> assembled_sfs;
  int i = sfs.size() - 1;
  while (i >= 0) {
    int j;
    for (j = i - 1; j >= 0; --j) {
      if (sfs[j + 1].s + sfs[j + 1].l <= sfs[j].s - 200) { // FIXME hardcoded
        // non-overlapping
        int l = sfs[j + 1].s + sfs[j + 1].l - sfs[i].s;
        assembled_sfs.push_back({sfs[i].s, l});
        i = j;
        break;
      }
    }
    if (j < 0) {
      int l = sfs[j + 1].s + sfs[j + 1].l - sfs[i].s;
      assembled_sfs.push_back({sfs[i].s, l});
      i = j;
    }
  }
  return assembled_sfs;
}

/* Compute SFS strings from P and store them into solutions*/
vector<sfs_t> ping_pong_search(const rb3_fmi_t *index, uint8_t *P, int l) {
  vector<sfs_t> S;
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
    S.push_back({begin, end - begin + 1});

    if (begin == 0)
      break;

    //   if (config->overlap == 0) // Relaxed
    //     begin -= 1;
    //   else
    begin = end - 1;
  }
  return S;
}

vector<sfs_t> extend(const vector<sfs_t> &sfs, uint8_t *P, int l, GSK &gsk
                     /* , const string &seq_name */) {
  vector<sfs_t> extended_sfs;
  int k = gsk.k;
  int b, e;
  int n = -1, sn = -1, en = -1;
  char *kmer = (char *)malloc((k + 1) * sizeof(char));
  kmer[k] = '\0';
  uint64_t kmer_d;       // kmer
  uint64_t rckmer_d = 0; // reverse and complemented kmer
  uint64_t ckmer_d = 0;  // canonical kmer
  int c;                 // current char

  // FIXME: hardcoded
  int N = 20;
  for (const sfs_t &s : sfs) {
    // cerr << s.s << " " << s.l << endl;
    b = s.s - k;
    b = b < 0 ? 0 : b;
    e = s.s + s.l;
    e = e > l - k + 1 ? l - k + 1 : e;

    memcpy(kmer, P + b, k);
    kmer_d = k2d(kmer, k);
    rckmer_d = rc(kmer_d, k);
    ckmer_d = std::min(kmer_d, rckmer_d);
    vector<pair<int, int>> snodes;
    while (b > 0 && snodes.size() < N) {
      if ((sn = gsk.get(ckmer_d)) != -1)
        snodes.push_back(make_pair(sn, b));
      --b;
      c = P[b] < 5 ? P[b] - 1 : rand() % 4;
      kmer_d = rsprepend(kmer_d, c, k);
      rckmer_d = lsappend(rckmer_d, reverse_char(c), k);
      ckmer_d = std::min(kmer_d, rckmer_d);
    }

    memcpy(kmer, P + e, k);
    kmer_d = k2d(kmer, k);
    rckmer_d = rc(kmer_d, k);
    ckmer_d = std::min(kmer_d, rckmer_d);
    vector<pair<int, int>> enodes;
    while (e < l - k + 1 && enodes.size() < N) {
      if ((en = gsk.get(ckmer_d)) != -1)
        enodes.push_back(make_pair(en, e));
      ++e;
      c = P[e] < 5 ? P[e] - 1 : rand() % 4;
      kmer_d = lsappend(kmer_d, c, k);
      rckmer_d = rsprepend(rckmer_d, reverse_char(c), k);
      ckmer_d = std::min(kmer_d, rckmer_d);
    }

    if (snodes.size() == 0 || enodes.size() == 0) {
      extended_sfs.push_back({b, e + k + 1 - b, -1, -1, 0});
    } else {
      // for (int i = 0; i < snodes.size(); ++i)
      //   cout << snodes[i].first << " ";
      // cout << " >>>  ";
      // for (int i = 0; i < enodes.size(); ++i)
      //   cout << enodes[i].first << " ";
      // cout << endl;

      int mind = 100;
      int d;
      for (int i = 0; i < snodes.size(); ++i) {
        for (int j = 0; j < enodes.size(); ++j) {
          d = abs(snodes[i].first - enodes[j].first);
          if (d < mind) {

            sn = snodes[i].first;
            b = snodes[i].second;
            en = enodes[j].first;
            e = enodes[j].second;
            mind = d;
          }
        }
      }
      if (sn <= en)
        extended_sfs.push_back({b, e + k + 1 - b, sn, en, 1});
      else
        extended_sfs.push_back({b, e + k + 1 - b, en, sn, 0});

      // if (sn != -1 && en != -1) {
      //   std::cout << (sn <= en ? sn : en) << ">" << (sn <= en ? en : sn) <<
      //   endl; CXXGraph::Node<int> n1(to_string(sn <= en ? sn : en), sn <= en
      //   ? sn : en); CXXGraph::Node<int> n2(to_string(sn <= en ? en : sn), sn
      //   <= en ? en : sn); std::cout << "d(" << n1 << "," << n2 << ") = " <<
      //   std::flush; auto res = gsk.graph.dijkstra(n1, n2); cout << res.result
      //   << "\n";
      // }
    }
  }
  free(kmer);
  return extended_sfs;
}

int main(int argc, char *argv[]) {
  if (strcmp(argv[1], "call") == 0)
    return main_call(argc - 1, argv + 1);

  char *gfa_fn = argv[1];
  char *fmd_fn = argv[2];
  char *fq_fn = argv[3];
  int k = stoi(argv[4]);

  cerr << "Building graph sketch..." << endl;
  GSK gsk(gfa_fn);
  gsk.build_sketch(k);
  gsk.build_graph();
  cerr << "Done." << endl;
  return 1;

  rb3_fmi_t f;
  rb3_fmi_restore(&f, fmd_fn, 0);
  if (f.e == 0 && f.r == 0) {
    cerr << "Error restoring index" << endl;
    return 1;
  }

  cerr << "Computing specific strings..." << endl;
  gzFile fp = gzopen(fq_fn, "r");
  kseq_t *seq = kseq_init(fp);
  int l;
  // char *ss = (char *)malloc(100000); // FIXME
  uint8_t *s;
  vector<sfs_t> S;
  while ((l = kseq_read(seq)) >= 0) {
    //   strncpy(ss, seq->seq.s, seq->seq.l + 1);
    s = (uint8_t *)seq->seq.s;
    rb3_char2nt6(seq->seq.l, s);
    S = ping_pong_search(&f, s, seq->seq.l);
    // for (const auto &s : S)
    //   cerr << s.s << "\t" << s.l << endl;
    // cerr << "---" << endl;
    S = assemble_rl(S);
    //   // for (const auto &s : S)
    //   //   cerr << s.s << "\t" << s.l << endl;

    S = extend(S, s, l, gsk);
    //   // for (const auto &s : S)
    //   //   cerr << s.s << "\t" << s.l << endl;
    // S = assemble_lr(S);
    for (const auto &s : S)
      cout << seq->name.s << "\t" << s.s << "\t" << s.l << "\t" << s.sn << ">"
           << s.en << "\t" << s.strand << "\t" << l << endl;
  }
  kseq_destroy(seq);
  gzclose(fp);
  rb3_fmi_free(&f);

  cerr << "END" << endl;
  return 0;
}
