#include <algorithm>
#include <bit>
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

struct anchor_t {
  int v = -1;        // vertex on graph
  int p = -1;        // position on query
  uint64_t seq = -1; // kmer
};

struct sfs_t {
  int qidx; // read index
  int s;    // start on query
  int l;    // length
  anchor_t a = {};
  anchor_t b = {};
  int strand = 1;
  uint64_t esk = -1,
           eek = -1; // expected starting and ending kmers (from cluster)
  int good = true;
  string seq = "";
};

struct cluster_t {
  vector<sfs_t> specifics;
  int va = -1, vb = -1;      // starting and ending vertices
  uint64_t ka = -1, kb = -1; // starting and ending kmers
};

vector<sfs_t> assemble(const vector<sfs_t> &sfs) {
  vector<sfs_t> assembled_sfs;
  int i = sfs.size() - 1;
  while (i >= 0) {
    int j;
    for (j = i - 1; j >= 0; --j) {
      if (sfs[j + 1].s + sfs[j + 1].l <= sfs[j].s - 200) { // FIXME hardcoded
        // non-overlapping
        int l = sfs[j + 1].s + sfs[j + 1].l - sfs[i].s;
        assembled_sfs.push_back({sfs[i].qidx, sfs[i].s, l});
        i = j;
        break;
      }
    }
    if (j < 0) {
      int l = sfs[j + 1].s + sfs[j + 1].l - sfs[i].s;
      assembled_sfs.push_back({sfs[i].qidx, sfs[i].s, l});
      i = j;
    }
  }
  return assembled_sfs;
}

/* Compute SFS strings from P and store them into solutions */
vector<sfs_t> ping_pong_search(const rb3_fmi_t *index, uint8_t *P, int qidx,
                               int l) {
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
    S.push_back({qidx, begin, end - begin + 1});

    if (begin == 0)
      break;
    //   if (config->overlap == 0) // Relaxed
    //     begin -= 1;
    //   else
    begin = end - 1;
  }
  return S;
}

vector<sfs_t> anchor(const vector<sfs_t> &sfs, uint8_t *P, int l, GSK &gsk) {
  vector<sfs_t> anchored_sfs;
  int k = gsk.k;
  int b, e;
  int vx = -1;
  char *kmer = (char *)malloc((k + 1) * sizeof(char));
  kmer[k] = '\0';
  uint64_t kmer_d;       // kmer
  uint64_t rckmer_d = 0; // reverse and complemented kmer
  uint64_t ckmer_d = 0;  // canonical kmer
  int c;                 // current char

  int N = 20; // FIXME: hardcoded
  for (const sfs_t &s : sfs) {
    b = s.s - k;
    b = b < 0 ? 0 : b;
    e = s.s + s.l;
    e = e > l - k + 1 ? l - k + 1 : e;
    /*cerr << "Extending from " << b << " " << e << endl;*/

    memcpy(kmer, P + b, k);
    kmer_d = k2d(kmer, k);
    rckmer_d = rc(kmer_d, k);
    ckmer_d = std::min(kmer_d, rckmer_d);
    vector<anchor_t> sanchors;
    /*cerr << "s " << b << " " << ckmer_d << endl;*/
    while (b > 0 && sanchors.size() < N) {
      if ((vx = gsk.get(ckmer_d)) != -1)
        sanchors.push_back({vx, b, ckmer_d});
      --b;
      c = P[b] < 5 ? P[b] - 1 : rand() % 4;
      kmer_d = rsprepend(kmer_d, c, k);
      rckmer_d = lsappend(rckmer_d, reverse_char(c), k);
      ckmer_d = std::min(kmer_d, rckmer_d);
      /*cerr << "s " << b << " " << ckmer_d << endl;*/
    }
    /*cerr << "-" << endl;*/

    memcpy(kmer, P + e, k);
    kmer_d = k2d(kmer, k);
    rckmer_d = rc(kmer_d, k);
    ckmer_d = std::min(kmer_d, rckmer_d);
    /*cerr << "e " << e << " " << ckmer_d << endl;*/
    vector<anchor_t> eanchors;
    while (e < l - k + 1 && eanchors.size() < N) {
      if ((vx = gsk.get(ckmer_d)) != -1)
        eanchors.push_back({vx, e, ckmer_d});
      ++e;
      c = P[e + k - 1] < 5 ? P[e + k - 1] - 1 : rand() % 4;
      kmer_d = lsappend(kmer_d, c, k);
      rckmer_d = rsprepend(rckmer_d, reverse_char(c), k);
      ckmer_d = std::min(kmer_d, rckmer_d);
      /*cerr << "e " << e << " " << ckmer_d << endl;*/
    }

    if (sanchors.size() == 0 || eanchors.size() == 0)
      // extended_sfs.push_back({b, e + k + 1 - b, -1, -1, 0});
      continue;
    // for (int i = 0; i < snodes.size(); ++i)
    // cout << snodes[i].first << " ";
    // cout << " >>>  ";
    // for (int i = 0; i < enodes.size(); ++i)
    //   cout << enodes[i].first << " ";
    // cout << endl;
    /*for (int i = 0; i < sanchors.size(); ++i)*/
    /*  cerr << "S " << sanchors[i].v << " " << sanchors[i].p << " " <<
     * sanchors[i].seq << endl;*/
    /*for (int i = 0; i < eanchors.size(); ++i)*/
    /*  cerr << "E " << eanchors[i].v << " " << eanchors[i].p << " " <<
     * eanchors[i].seq << endl;*/
    /*cerr << endl;*/
    int mind = 100;
    int d;
    int sax = -1, eax = -1; // index for selected anchors
    for (int i = 0; i < sanchors.size(); ++i) {
      for (int j = 0; j < eanchors.size(); ++j) {
        d = abs(sanchors[i].v - eanchors[j].v);
        if (d < mind) {
          sax = i;
          eax = j;
          mind = d;
        }
      }
    }
    if (sax == -1 || eax == -1)
      continue;

    anchor_t sa = sanchors[sax];
    anchor_t ea = eanchors[eax];
    int b = sa.p;
    int l = ea.p + k + 1 - sa.p;
    int strand = 1;
    if (sa.v > ea.v) {
      anchor_t tmp = sa;
      sa = ea;
      ea = tmp;
      strand = 0;
    }
    anchored_sfs.push_back({s.qidx, b, l, sa, ea, strand});
    assert(sa.v <= ea.v);
    // if (sn != -1 && en != -1) {
    //   std::cout << (sn <= en ? sn : en) << ">" << (sn <= en ? en : sn) <<
    //   endl; CXXGraph::Node<int> n1(to_string(sn <= en ? sn : en), sn <= en
    //   ? sn : en); CXXGraph::Node<int> n2(to_string(sn <= en ? en : sn), sn
    //   <= en ? en : sn); std::cout << "d(" << n1 << "," << n2 << ") = " <<
    //   std::flush; auto res = gsk.graph.dijkstra(n1, n2); cout << res.result
    //   << "\n";
    // }
  }
  free(kmer);

  // CHECKME: specifics strings are already sorted by position on query
  // std::sort(extended_sfs.begin(), extended_sfs.end(),
  //           [](const sfs_t &a, const sfs_t &b) {
  //             return a.s < b.s;
  //           });
  return anchored_sfs;
}

void add(cluster_t &c, const sfs_t &s) {
  c.specifics.push_back(s);
  if (c.va == -1 || s.a.v < c.va) {
    c.va = s.a.v;
    c.ka = s.a.seq;
  }
  if (c.vb == -1 || s.b.v > c.vb) {
    c.vb = s.b.v;
    c.kb = s.b.seq;
  }
}

// sweep line clustering
vector<cluster_t> cluster(const vector<sfs_t> SS, const GSK &gsk) {
  vector<cluster_t> clusters(1);
  add(clusters.back(), SS[0]);
  for (int i = 1; i < SS.size(); ++i) {
    if (SS[i].a.v > clusters.back().vb) {
      // no overlap, so new cluster
      clusters.push_back({});
    }
    add(clusters.back(), SS[i]);
  }
  return clusters;
}

void merge(cluster_t &C) {
  vector<sfs_t> newC;
  map<int, vector<sfs_t>> byread;
  for (const sfs_t &s : C.specifics)
    byread[s.qidx].push_back(s);

  for (const auto &c : byread) {
    int mins = 100000, maxe = 0; // FIXME: assuming HiFi
    int first, last;
    for (int i = 0; i < c.second.size(); ++i) {
      if (c.second[i].s < mins) {
        mins = c.second[i].s;
        first = i;
      }
      if (c.second[i].s + c.second[i].l > maxe) {
        maxe = c.second[i].s + c.second[i].l;
        last = i;
      }
    }
    /*cerr << c.first << endl;*/
    /*cerr << c.second[first].strand << " " << c.second[last].strand << endl;*/
    /*cerr << c.second[first].a.p << " " << c.second[first].a.v << " : " <<
     * c.second[first].b.p << " " << c.second[first].b.v << endl;*/
    /*cerr << c.second[last].a.p << " " << c.second[last].a.v << " : " <<
     * c.second[last].b.p << " " << c.second[last].b.v << endl;*/
    /*cerr << endl;*/
    newC.push_back({c.first, c.second[first].s,
                    c.second[last].s + c.second[last].l - c.second[first].s,
                    c.second[0].strand ? c.second[first].a : c.second[last].a,
                    c.second[0].strand ? c.second[last].b : c.second[first].b,
                    c.second[first].strand});
    assert(newC.back().l > 0);
    assert(newC.back().a.v <= newC.back().b.v);
  }
  C.specifics = newC;
}

string d2s(uint64_t kmer, int k) {
  char kk[k + 1];
  for (int i = 1; i <= k; ++i)
    kk[i - 1] = "ACGT"[(kmer >> (k - i) * 2) & 3];
  kk[k] = '\0';
  return kk;
}

int main(int argc, char *argv[]) {
  if (strcmp(argv[1], "call") == 0)
    return main_call(argc - 1, argv + 1);

  char *gfa_fn = argv[1];
  char *fmd_fn = argv[2];
  char *fq_fn = argv[3];
  int k = stoi(argv[4]);
  int w = 2;  // FIXME: hardcoded, minimum weight for clusters
  int hd = 0; // FIXME: hardcoded, hamming distance

  cerr << "Building graph sketch..." << endl;
  GSK gsk(gfa_fn);
  gsk.build_sketch(k);
  // gsk.build_graph();
  cerr << "Done." << endl;

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
  vector<sfs_t> SS;
  uint qidx = 0;
  vector<string> qnames;
  vector<int> strands(2);
  int strand;
  while ((l = kseq_read(seq)) >= 0) {
    //   strncpy(ss, seq->seq.s, seq->seq.l + 1);
    s = (uint8_t *)seq->seq.s;
    rb3_char2nt6(seq->seq.l, s);

    S = ping_pong_search(&f, s, qidx, seq->seq.l);
    /*for (const auto &s : S)*/
    /*  cout << seq->name.s << " " << s.s << " " << s.l << endl;*/
    /*cout << "_" << endl;*/
    S = assemble(S);
    /*for (const auto &s : S)*/
    /*  cout << seq->name.s << " " << s.s << " " << s.l << endl;*/
    /*cout << "_" << endl;*/
    S = anchor(S, s, l, gsk);

    strands[0] = 0;
    strands[1] = 0;
    for (const auto &s : S)
      ++strands[s.strand];
    // FIXME: plus strand if tie
    strand = 1;
    if (strands[0] > strands[1])
      strand = 0;
    for (const auto &s : S)
      if (s.strand == strand)
        SS.push_back(s);
    // SS.reserve(SS.size() + strands[strand];
    // SS.insert(SS.end(), S.begin(), S.end());

    qnames.push_back(seq->name.s);
    ++qidx;
  }
  for (const sfs_t &s : SS)
    assert(s.l > 0);
  kseq_destroy(seq);
  gzclose(fp);
  rb3_fmi_free(&f);

  cerr << "Clustering " << SS.size() << " specific strings" << endl;
  std::sort(SS.begin(), SS.end(),
            [](const sfs_t &a, const sfs_t &b) { return a.a.v < b.b.v; });
  vector<cluster_t> Cs = cluster(SS, gsk);

  cerr << "Merging " << Cs.size() << " clusters" << endl;
  for (auto &c : Cs) {
    if (c.specifics.size() < w) {
      c.specifics.clear();
    } else {
      merge(c);
      if (c.specifics.size() < w)
        c.specifics.clear();
    }
  }

  // checking if specifics start/end with "correct" kmers
  vector<vector<sfs_t *>> pspecifics(qnames.size());
  for (auto &c : Cs) {
    if (c.specifics.empty())
      continue;
    for (const sfs_t &s : c.specifics)
      assert(s.l > 0);
    cout << c.specifics.size() << " " << c.va << ">" << c.vb << " "
         << (c.va - 1) * 32 << " " << c.ka << " " << c.kb << endl;
    for (auto &s : c.specifics) {
      if (s.a.seq == c.ka && s.b.seq == c.kb)
        continue;
      if (s.strand) {
        s.esk = c.ka;
        s.eek = c.kb;
      } else {
        s.esk = c.kb;
        s.eek = c.ka;
      }
      pspecifics[s.qidx].push_back(&s);
      cout << qnames[s.qidx] << ":" << s.s << "-" << s.s + s.l << " " << s.l
           << " " << s.strand << " " << s.a.v << ">" << s.b.v << " " << s.a.seq
           << " " << s.b.seq << endl;
    }
  }
  cout << "---" << endl;
  fp = gzopen(fq_fn, "r");
  seq = kseq_init(fp);
  qidx = 0;

  char *kmer = (char *)malloc((k + 1) * sizeof(char));
  kmer[k] = '\0';
  uint64_t kmer_d;       // kmer
  uint64_t rckmer_d = 0; // reverse and complemented kmer
  uint64_t ckmer_d = 0;  // canonical kmer
  int c;                 // current char
  int p;
  int pc;
  while ((l = kseq_read(seq)) >= 0) {
    s = (uint8_t *)seq->seq.s;
    rb3_char2nt6(seq->seq.l, s);

    for (sfs_t *s : pspecifics[qidx]) {
      if (s->a.seq == s->esk && s->b.seq == s->eek)
        continue;
      /*cerr << "Extending " << s << " on " << seq->name.s << " " << s->s << " "
       * << s->l << " " << s->strand << endl;*/

      if (s->a.seq != s->esk) {
        p = s->s - 1;
        memcpy(kmer, seq->seq.s + p, k);
        kmer_d = k2d(kmer, k);
        rckmer_d = rc(kmer_d, k);
        ckmer_d = std::min(kmer_d, rckmer_d);
        // cerr << "B " << p << " : " << s->esk << " " << ckmer_d << " " <<
        // d2s(ckmer_d, k) << endl;
        while (p > 0 && (pc = popcount(ckmer_d ^ s->esk)) > hd) {
          --p;
          c = seq->seq.s[p] < 5 ? seq->seq.s[p] - 1 : rand() % 4;
          // cerr << "ACGT"[c] << endl;
          kmer_d = rsprepend(kmer_d, c, k);
          rckmer_d = lsappend(rckmer_d, reverse_char(c), k);
          ckmer_d = std::min(kmer_d, rckmer_d);
          // cerr << "B " << p << " : " << s->esk << " " << ckmer_d << " " <<
          // d2s(ckmer_d, k) << endl;
        }
        if (pc > hd) {
          s->good = 0;
        } else {
          s->a.seq = ckmer_d; // FIXME: change also the vertex or don't change
                              // anything and mark the anchors as not updated
          s->l += s->s - p;
          s->s = p;
        }
      }
      if (s->b.seq != s->eek) {
        p = s->s + s->l;
        memcpy(kmer, seq->seq.s + p, k);
        kmer_d = k2d(kmer, k);
        rckmer_d = rc(kmer_d, k);
        ckmer_d = std::min(kmer_d, rckmer_d);
        while (p < l - k + 1 && (pc = popcount(ckmer_d ^ s->eek)) > hd) {
          ++p;
          c = seq->seq.s[p + k - 1] < 5 ? seq->seq.s[p + k - 1] - 1
                                        : rand() % 4;
          kmer_d = lsappend(kmer_d, c, k);
          rckmer_d = rsprepend(rckmer_d, reverse_char(c), k);
          ckmer_d = std::min(kmer_d, rckmer_d);
        }
        if (pc > hd) {
          s->good = 0;
        } else {
          s->b.seq = ckmer_d; // FIXME: change also the vertex or don't change
                              // anything and mark the anchors as not updated
          s->l += p + k - (s->s + s->l) + 1; // CHECKME this +1
        }
      }
    }
    ++qidx;
  }
  free(kmer);
  kseq_destroy(seq);
  gzclose(fp);

  // checking if specifics start/end with "correct" kmers
  for (auto &c : Cs) {
    if (c.specifics.empty())
      continue;

    cout << c.specifics.size() << " " << c.va << ">" << c.vb << " "
         << (c.va - 1) * 32 << "-" << (c.vb) * 32 << " " << c.ka << " " << c.kb
         << endl;
    for (auto &s : c.specifics) {
      cout << s.good << " " << qnames[s.qidx] << ":" << s.s << "-" << s.s + s.l
           << " " << s.l << " " << s.strand << " " << s.a.v << ">" << s.b.v
           << " " << s.a.seq << " " << s.b.seq << endl;
    }
  }

  cerr << "END" << endl;
  return 0;
}
