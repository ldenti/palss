#include <algorithm>
#include <bit>
#include <cstdint>
#include <iostream>
#include <vector>
#include <zlib.h>

#include "abpoa.h"
#include "fm-index.h"
#include "kseq.h"
#include "ksw2.h"

#include "gsketch.hpp"
#include "utils.h"

/*int main_index(int argc, char *argv[]);*/

// KSEQ_INIT(gzFile, gzread) // we already init kstream in gsketch
__KSEQ_TYPE(gzFile)
__KSEQ_BASIC(static, gzFile)
__KSEQ_READ(static)

using namespace std;

struct anchor_t {
  int64_t v = -1;    // vertex on graph
  int offset = -1;   // offset on vertex
  int p = -1;        // position on query
  uint64_t seq = -1; // kmer
};

struct sfs_t {
  int qidx;        // read index
  int s;           // start on query
  int l;           // length
  anchor_t a = {}; // left anchor
  anchor_t b = {}; // right anchor
  int strand = 1;  // inferred strand
  uint64_t esk = -1,
           eek = -1; // expected starting and ending kmers (from cluster)
  int good = 1;      // is it good for calling step?
  char *seq = NULL;  // sequence
};

struct cluster_t {
  vector<sfs_t> specifics;
  int va = -1, vb = -1;      // starting and ending vertices
  int offa = -1, offb = -1;  // offsets on the two vertices
  uint64_t ka = -1, kb = -1; // starting and ending kmers
};

string decode(const char *s, int l, int shift) {
  if (s == NULL)
    return "";
  char ds[l + 1];
  for (int i = 0; i < l; ++i)
    ds[i] = "NACGTN"[s[i] + shift];
  ds[l] = '\0';
  return ds;
}

int build_consensus(abpoa_t *ab, abpoa_para_t *abpt,
                    const vector<sfs_t *> &specifics, char **cons,
                    int *cons_c) {
  int *seq_lens = (int *)malloc(sizeof(int) * specifics.size());
  uint8_t **bseqs = (uint8_t **)malloc(sizeof(uint8_t *) * specifics.size());
  int goods = 0, i = 0;
  for (const sfs_t *s : specifics) {
    assert(s->good);
    /*if (!s.good)*/
    /*  continue;*/
    seq_lens[goods] = s->l;
    bseqs[goods] = (uint8_t *)malloc(sizeof(uint8_t) * (s->l + 1));
    for (i = 0; i < s->l; ++i)
      bseqs[goods][i] = s->seq[i] - 1;
    bseqs[goods][s->l] = '\0';

    if (!s->strand) {
      // rc
      for (i = 0; i < (s->l >> 1); ++i) {
        int tmp = bseqs[goods][s->l - 1 - i];
        tmp = (tmp >= 0 && tmp <= 3) ? 3 - tmp : tmp;
        bseqs[goods][s->l - 1 - i] =
            (bseqs[goods][i] >= 0 && bseqs[goods][i] <= 3) ? 3 - bseqs[goods][i]
                                                           : bseqs[goods][i];
        bseqs[goods][i] = tmp;
      }
      if (s->l & 1)
        bseqs[goods][i] = (bseqs[goods][i] >= 0 && bseqs[goods][i] <= 3)
                              ? 3 - bseqs[goods][i]
                              : bseqs[goods][i];
    }
    ++goods;
  }
  assert(goods > 0);
  /*if (!goods) {*/
  /*  (*cons)[0] = '\0';*/
  /*  return 0;*/
  /*}*/
  /*for (int i = 0; i < goods; ++i) {*/
  /*  cout << seq_lens[i] << " " << decode((char *)bseqs[i], seq_lens[i], 1)*/
  /*       << endl;*/
  /*}*/
  abpoa_msa(ab, abpt, goods, NULL, seq_lens, bseqs, NULL, NULL);
  abpoa_cons_t *abc = ab->abc;
  int cons_l = 0;
  if (abc->n_cons > 0) {
    cons_l = abc->cons_len[0];
    if (cons_l + 1 > *cons_c) {
      fprintf(stderr, "--- Reallocating from %d to %d\n", *cons_c,
              (cons_l + 1) * 2);
      char *temp = (char *)realloc(*cons, (cons_l + 1) * 2 * sizeof(char));
      if (temp == NULL) {
        free(cons);
        fprintf(stderr,
                "Error while reallocating memory for consensus string\n");
        exit(2);
      } else {
        *cons = temp;
      }
      *cons_c = (cons_l + 1) * 2;
    }
    for (i = 0; i < cons_l; ++i)
      (*cons)[i] = abc->cons_base[0][i]; // "ACGTN"[abc->cons_base[0][i]];
    (*cons)[i] = '\0';
  }

  for (i = 0; i < goods; ++i)
    free(bseqs[i]);
  free(bseqs);
  free(seq_lens);

  return cons_l;
}

/* Merge specifics strings that are too close on the same read */
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

/* Anchor specific strings on graph using graph sketch */
vector<sfs_t> anchor(const vector<sfs_t> &sfs, uint8_t *P, int l, int N,
                     GSK &gsk) {
  vector<sfs_t> anchored_sfs;
  int k = gsk.klen;
  int b, e;
  pair<int64_t, uint16_t> vx = make_pair(-1, -1);
  char *kmer = (char *)malloc((k + 1) * sizeof(char));
  kmer[k] = '\0';
  uint64_t kmer_d;       // kmer
  uint64_t rckmer_d = 0; // reverse and complemented kmer
  uint64_t ckmer_d = 0;  // canonical kmer
  int c;                 // current char

  for (const sfs_t &s : sfs) {
    b = s.s - k;
    b = b < 0 ? 0 : b;
    e = s.s + s.l;
    e = e > l - k ? l - k : e;

    memcpy(kmer, P + b, k);
    kmer_d = k2d(kmer, k);
    rckmer_d = rc(kmer_d, k);
    ckmer_d = std::min(kmer_d, rckmer_d);
    vector<anchor_t> sanchors;
    while (b > 0 && sanchors.size() < N) {
      if ((vx = gsk.get(ckmer_d)).first != -1)
        sanchors.push_back({vx.first, vx.second, b, ckmer_d});
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
    vector<anchor_t> eanchors;
    while (e < l - k + 1 && eanchors.size() < N) {
      if ((vx = gsk.get(ckmer_d)).first != -1)
        eanchors.push_back({vx.first, vx.second, e, ckmer_d});
      ++e;
      c = P[e + k - 1] < 5 ? P[e + k - 1] - 1 : rand() % 4;
      kmer_d = lsappend(kmer_d, c, k);
      rckmer_d = rsprepend(rckmer_d, reverse_char(c), k);
      ckmer_d = std::min(kmer_d, rckmer_d);
    }

    if (sanchors.size() == 0 || eanchors.size() == 0)
      // extended_sfs.push_back({b, e + k + 1 - b, -1, -1, 0});
      continue;

    int mind = 100;
    int d;
    int sax = -1, eax = -1; // index for selected anchors
    map<pair<int, int>, int> memo;
    map<pair<int, int>, int>::iterator hhit;
    int comp;
    int x = 0, xoff = 0, y = 0, yoff = 0;
    pair<int, int> xy = {x, y};
    for (uint i = 0; i < sanchors.size(); ++i) {
      x = sanchors[i].v;
      xoff = sanchors[i].offset;
      xy.first = x;
      for (int j = 0; j < eanchors.size(); ++j) {
        y = eanchors[j].v;
        yoff = eanchors[j].offset;
        xy.second = y;
        if ((hhit = memo.find(xy)) == memo.end()) {
          memo[xy] = gsk.compatible(sanchors[i].v, eanchors[j].v);
          /*memo[make_pair(y, x)] = memo[xy];*/
        }
        comp = memo[xy];
        if (!comp)
          continue;
        if (x == y && (xoff == yoff || (xoff < yoff && xoff + k >= yoff) ||
                       (xoff > yoff && yoff + k >= xoff)))
          continue;
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
    int l = ea.p + k - sa.p;
    int strand = 1;
    if (sa.v > ea.v || (sa.v == ea.v && sa.offset > ea.offset)) {
      anchor_t tmp = sa;
      sa = ea;
      ea = tmp;
      strand = 0;
    }
    anchored_sfs.push_back({s.qidx, b, l, sa, ea, strand});
    assert(sa.v <= ea.v);
  }
  free(kmer);

  // CHECKME: specifics strings are already sorted by position on query
  // std::sort(extended_sfs.begin(), extended_sfs.end(),
  //           [](const sfs_t &a, const sfs_t &b) {
  //             return a.s < b.s;
  //           });
  return anchored_sfs;
}

/* Add specific string s to cluster c and update its anchors */
void add(cluster_t &c, const sfs_t &s) {
  c.specifics.push_back(s);
  if (c.va == -1 || s.a.v < c.va) {
    c.va = s.a.v;
    c.offa = s.a.offset;
    c.ka = s.a.seq;
  }
  if (c.vb == -1 || s.b.v > c.vb) {
    c.vb = s.b.v;
    c.offb = s.b.offset;
    c.kb = s.b.seq;
  }
}

/* Sweep line clustering of specific strings based on their subgraphs - assuming
 * topological sorted DAG */
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
    int first = -1, last = -1;
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
    if (first == -1 || last == -1) {
      assert(false);
      continue;
    }
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

int align(char *tseq, char *qseq) {
  int sc_mch = 2, sc_mis = -4, gapo = 6, gape = 1;
  int i;
  int8_t a = sc_mch, b = sc_mis < 0 ? sc_mis : -sc_mis; // a>0 and b<0
  int8_t mat[25] = {a, b, b, b, 0, b, a, b, b, 0, b, b, a,
                    b, 0, b, b, b, a, 0, 0, 0, 0, 0, 0};
  int tl = strlen(tseq), ql = strlen(qseq);
  uint8_t *ts, *qs, c[256];
  ksw_extz_t ez;

  memset(&ez, 0, sizeof(ksw_extz_t));
  memset(c, 4, 256);
  c['A'] = c['a'] = 0;
  c['C'] = c['c'] = 1;
  c['G'] = c['g'] = 2;
  c['T'] = c['t'] = 3; // build the encoding table
  ts = (uint8_t *)malloc(tl);
  qs = (uint8_t *)malloc(ql);
  for (i = 0; i < tl; ++i)
    ts[i] = c[(uint8_t)tseq[i]]; // encode to 0/1/2/3
  for (i = 0; i < ql; ++i)
    qs[i] = c[(uint8_t)qseq[i]];
  // ksw_extz(0, ql, qs, tl, ts, 5, mat, gapo, gape, -1, -1, 0, &ez);

  ksw_extz2_sse(0, ql, qs, tl, ts, 5, mat, gapo, gape, -1, -1, 0, 0, &ez);

  for (i = 0; i < ez.n_cigar; ++i) // print CIGAR
    printf("%d%c", ez.cigar[i] >> 4, "MID"[ez.cigar[i] & 0xf]);
  putchar('\n');
  free(ez.cigar);
  free(ts);
  free(qs);
  return 0;
}

int main(int argc, char *argv[]) {
  if (strcmp(argv[1], "align") == 0)
    return 1; // main_index(argc-1, argv+1);
  else if (strcmp(argv[1], "align") == 0)
    return align(argv[2], argv[3]);
  char *gfa_fn = argv[1];
  char *fmd_fn = argv[2];
  char *fq_fn = argv[3];
  int k = stoi(argv[4]);
  int w = 2;     // FIXME: hardcoded, minimum weight for clusters
  int hd = 0;    // FIXME: hardcoded, hamming distance for fixing anchors
  int minl = 50; // FIXME: hardcoded, minimum SV length
  int N = 20;    // FIXME: hardcoded, number of kmers to check for anchoring
  double rt0, rt, rt1;
  rt0 = realtime();
  rt = rt0;

  // Graph sketching and path extraction
  GSK gsk(gfa_fn, k);
  gsk.build_sketch();
  fprintf(stderr, "[M::%s] sketched graph with %d vertices in %.3f sec\n",
          __func__, gsk.nvertices, realtime() - rt);
  rt = realtime();

  gsk.build_graph();
  fprintf(stderr, "[M::%s] loaded %ld paths in %.3f sec\n", __func__,
          gsk.paths.size(), realtime() - rt);
  rt = realtime();
  // ---

  // FMD-index loading
  rb3_fmi_t f;
  rb3_fmi_restore(&f, fmd_fn, 0);
  if (f.e == 0 && f.r == 0) {
    fprintf(stderr, "Error restoring index");
    return 1;
  }
  fprintf(stderr, "[M::%s] restored FMD index in %.3f sec\n", __func__,
          realtime() - rt);
  rt = realtime();
  rt1 = rt;
  // ---

  // Specific strings computation and anchoring
  gzFile fp = gzopen(fq_fn, "r");
  kseq_t *seq = kseq_init(fp);
  int l;
  uint8_t *s;
  vector<sfs_t> S;
  vector<sfs_t> SS;
  uint qidx = 0;
  vector<string> qnames;
  vector<int> strands(2);
  int strand;
  while ((l = kseq_read(seq)) >= 0) {
    s = (uint8_t *)seq->seq.s;
    rb3_char2nt6(seq->seq.l, s);

    S = ping_pong_search(&f, s, qidx, seq->seq.l);
    S = assemble(S);

    for (int i = 0; i < (int)S.size() - 1; ++i)
      assert(S[i].s < S[i + 1].s);

    S = anchor(S, s, l, N, gsk);
    for (int i = 0; i < (int)S.size() - 1; ++i)
      assert(S[i].s < S[i + 1].s);

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

    qnames.push_back(seq->name.s);
    ++qidx;
    if (qidx % 10000 == 0) {
      fprintf(stderr, "[M::%s] parsed %d reads %.3f sec\n", __func__, qidx,
              realtime() - rt1);
      rt1 = realtime();
    }
  }
  for (const sfs_t &s : SS)
    assert(s.l > 0);
  kseq_destroy(seq);
  gzclose(fp);
  rb3_fmi_free(&f);

  fprintf(stderr, "[M::%s] computed %ld specific strings in %.3f sec\n",
          __func__, SS.size(), realtime() - rt);
  rt = realtime();

  std::sort(SS.begin(), SS.end(), [](const sfs_t &a, const sfs_t &b) {
    return a.a.v < b.a.v || (a.a.v == b.a.v && a.a.offset < b.a.offset);
  });
  fprintf(stderr, "[M::%s] sorted specific strings in %.3f sec\n", __func__,
          realtime() - rt);
  rt = realtime();
  // ---

  // Clustering
  vector<cluster_t> Cs = cluster(SS, gsk);
  fprintf(stderr, "[M::%s] created %ld clusters in %.3f sec\n", __func__,
          Cs.size(), realtime() - rt);
  rt = realtime();

  int lowsc_n = 0;
  int lowsc_am_n = 0;
  int sc_n = 0;
  for (auto &c : Cs) {
    if (c.specifics.size() < w) {
      ++lowsc_n;
      c.specifics.clear();
    } else {
      merge(c);
      if (c.specifics.size() < w) {
        c.specifics.clear();
        ++lowsc_am_n;
      } else {
        ++sc_n;
      }
    }
  }
  fprintf(stderr,
          "[M::%s] merged clusters in %.3f sec (%d+%d=%d clusters filtered)\n",
          __func__, realtime() - rt, lowsc_n, lowsc_am_n, lowsc_n + lowsc_am_n);
  rt = realtime();
  // ---

  // Specific strings cleaning (make sure they start/end with "correct" kmers)
  vector<vector<sfs_t *>> pspecifics(qnames.size());
  for (auto &c : Cs) {
    if (c.specifics.empty())
      continue;
    for (const sfs_t &s : c.specifics)
      assert(s.l > 0);
    for (auto &s : c.specifics) {
      if (s.a.seq == c.ka && s.b.seq == c.kb) {
        s.good = 1;
        s.esk = c.ka;
        s.eek = c.kb;
      } else {
        s.good = 0;
        if (s.strand) {
          s.esk = c.ka;
          s.eek = c.kb;
        } else {
          s.esk = c.kb;
          s.eek = c.ka;
        }
      }
      pspecifics[s.qidx].push_back(&s);
    }
  }

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
      if (s->a.seq != s->esk && s->s > 0) {
        p = s->s - 1;
        memcpy(kmer, seq->seq.s + p, k);
        kmer_d = k2d(kmer, k);
        rckmer_d = rc(kmer_d, k);
        ckmer_d = std::min(kmer_d, rckmer_d);
        pc = hd + 1;
        while (p > 0 && (pc = popcount(ckmer_d ^ s->esk)) > hd) {
          --p;
          c = seq->seq.s[p] < 5 ? seq->seq.s[p] - 1 : rand() % 4;
          kmer_d = rsprepend(kmer_d, c, k);
          rckmer_d = lsappend(rckmer_d, reverse_char(c), k);
          ckmer_d = std::min(kmer_d, rckmer_d);
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
      if (s->b.seq != s->eek && s->s + s->l <= l - k) {
        p = s->s + s->l;
        memcpy(kmer, seq->seq.s + p, k);
        kmer_d = k2d(kmer, k);
        rckmer_d = rc(kmer_d, k);
        ckmer_d = std::min(kmer_d, rckmer_d);
        pc = hd + 1;
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
          s->l += p + k - (s->s + s->l);
        }
      }
      s->good = s->a.seq == s->esk && s->b.seq == s->eek;
      assert(s->good >= 0 && s->good <= 2);
      if (s->good > 0) {
        s->seq = (char *)malloc((s->l + 1) * sizeof(char));
        memcpy(s->seq, seq->seq.s + s->s, s->l);
        s->seq[s->l] = '\0';
      }
    }
    ++qidx;
  }
  free(kmer);
  kseq_destroy(seq);
  gzclose(fp);
  fprintf(stderr, "[M::%s] cleaned %d clusters in %.3f sec\n", __func__, sc_n,
          realtime() - rt);
  rt = realtime();
  rt1 = rt;
  // ---

  // Calling SVs
  char *pseq = (char *)malloc(8192 * sizeof(char));
  int pseq_c = 8192;
  int pseq_l = 0;
  char *cons = (char *)malloc(16384 * sizeof(char));
  int cons_c = 16384;
  int cons_l = 0;

  int sc_mch = 2, sc_mis = -4, gapo = 6, gape = 1;
  int8_t a = (int8_t)sc_mch,
         b = sc_mis < 0 ? (int8_t)sc_mis : -(int8_t)sc_mis; // a>0 and b<0
  int8_t mat[25] = {a, b, b, b, 0, b, a, b, b, 0, b, b, a,
                    b, 0, b, b, b, a, 0, 0, 0, 0, 0, 0};

  // INIT ABPOA
  abpoa_t *ab = abpoa_init();
  abpoa_para_t *abpt = abpoa_init_para();
  abpt->disable_seeding = 1;
  abpt->align_mode = 0; // global
  abpt->out_msa = 0;
  abpt->out_cons = 1;
  abpt->out_gfa = 0;
  // abpt->is_diploid = 1; // TODO: maybe this works now
  abpt->progressive_poa = 1;
  abpt->amb_strand = 0;
  abpt->wb = -1;
  abpt->max_n_cons = 1; // to generate 1 consensus sequences
  // abpt->match = 2;      // match score
  // abpt->mismatch = 4;   // mismatch penalty
  // abpt->gap_mode = ABPOA_CONVEX_GAP; // gap penalty mode
  // abpt->gap_open1 = 4;  // gap open penalty #1
  // abpt->gap_ext1 = 2;   // gap extension penalty #1
  // abpt->gap_open2 = 24; // gap open penalty #2
  // abpt->gap_ext2 = 1;   // gap extension penalty #2
  // gap_penalty = min{gap_open1 + gap_len*gap_ext1, gap_open2+gap_len*gap_ext2}
  abpoa_post_set_para(abpt);

  int vuidx = 0;
  char *cigar = (char *)malloc(16384 * sizeof(char)); // FIXME: hardcoded
  char *ins = (char *)malloc(100000 * sizeof(char));  // FIXME: hardcoded
  char *del = (char *)malloc(100000 * sizeof(char));  // FIXME: hardcoded
  int cp, cc = 0;
  for (auto &c : Cs) {
    if (c.specifics.empty())
      continue;
    ++cc;
    if (cc % 5000 == 0) {
      fprintf(stderr, "[M::%s] analyzed %d clusters (%d left) in %.3f sec\n",
              __func__, cc, sc_n - cc, realtime() - rt1);
      rt1 = realtime();
    }
    /*cout << c.specifics.size() << " " << c.va << ">" << c.vb << " "*/
    /*     << (c.va - 1) * 512 << "-" << (c.vb) * 512 << endl;*/
    /*for (const auto &s : c.specifics)*/
    /*  if (s.good)*/
    /*    cout << qnames[s.qidx] << " " << s.s << " " << s.l << " "*/
    /*         << decode(s.seq, s.l, 0) << endl;*/

    // split cluster based on strings length
    vector<int> subclusters_l(1);
    vector<vector<sfs_t *>> subclusters(1);
    int i = 0;
    while (i < c.specifics.size() && !c.specifics[i].good)
      ++i;
    if (i == c.specifics.size())
      continue;

    subclusters_l.back() = c.specifics[i].l;
    subclusters.back().push_back(&c.specifics[i]);
    for (i = i + 1; i < c.specifics.size(); ++i) {
      if (!c.specifics[i].good)
        continue;
      int j = 0;
      for (j = 0; j < subclusters.size(); ++j) {
        if (min(c.specifics[i].l, subclusters_l[j]) /
                (float)max(c.specifics[i].l, subclusters_l[j]) >=
            0.9) { // FIXME: hardcoded
          break;
        }
      }
      if (j == subclusters.size()) {
        subclusters.push_back({&c.specifics[i]});
        subclusters_l.push_back(c.specifics[i].l);
      } else {
        subclusters_l[j] =
            (subclusters_l[j] * subclusters.size() + c.specifics[i].l) /
            (subclusters.size() + 1);
        subclusters[j].push_back(&c.specifics[i]);
      }
    }

    sort(subclusters.begin(), subclusters.end(),
         [](const vector<sfs_t *> &a, const vector<sfs_t *> &b) {
           return a.size() > b.size();
         });

    vector<path_t *> subpaths = gsk.get_subpaths(c.va, c.vb);
    vector<pair<path_t *, string>> collapsed_subpaths;
    for (path_t *p : subpaths) {
      if (p == NULL)
        continue;
      int i = 0;
      for (const auto &cp : collapsed_subpaths) {
        if (p->l != cp.first->l)
          continue;
        int j = 0;
        for (j = 0; i < p->l; ++j) {
          if (p->vertices[j] != cp.first->vertices[j])
            break;
        }
        if (j == p->l)
          break;
        ++i;
      }
      if (i == collapsed_subpaths.size())
        collapsed_subpaths.push_back(make_pair(p, p->idx));
      else
        collapsed_subpaths[i].second += "," + string(p->idx);
    }
    for (int i = 0; i < (subclusters.size() == 1 ? 1 : 2); ++i) {
      cons_l = build_consensus(ab, abpt, subclusters[i], &cons, &cons_c);
      if (cons_l == 0)
        continue;
      for (pair<path_t *, string> collp : collapsed_subpaths) {
        path_t *p = collp.first;
        // if (p == NULL)
        //   continue;
        pseq_l = gsk.get_sequence(p, &pseq, &pseq_c);
        /*cout << cons_l << " " << cons_l << " " << decode(cons, cons_l, 1)*/
        /*     << endl;*/
        /*cout << pseq_l << " "*/
        /*     << pseq_l - c.offa - (gsk.get_vl(c.vb) - c.offb - k) << " "*/
        /*     << decode(pseq + c.offa,*/
        /*               pseq_l - c.offa - (gsk.get_vl(c.vb) - c.offb - k), 1)*/
        /*     << endl;*/
        ksw_extz_t ez;
        memset(&ez, 0, sizeof(ksw_extz_t));
        ksw_extz2_sse(0, cons_l, (uint8_t *)cons,
                      pseq_l - c.offa - (gsk.get_vl(c.vb) - c.offb - k),
                      (uint8_t *)(pseq + c.offa), 5, mat, gapo, gape, -1, -1, 0,
                      0, &ez);

        // OUTPUT
        int opl;
        int op;
        int qp = 0;
        int tp = 0;
        cp = 0;
        for (int i = 0; i < ez.n_cigar; ++i) {
          opl = ez.cigar[i] >> 4;
          cp += sprintf(cigar + cp, "%d%c", opl, "MID"[ez.cigar[i] & 0xf]);
        }

        for (int i = 0; i < ez.n_cigar; ++i) {
          l = ez.cigar[i] >> 4;
          op = ez.cigar[i] & 0xf; // 0:M, 1:I, 2:D
          if (op == 0) {
            // M
            qp += l;
            tp += l;
          } else if (op == 1) {
            // I
            if (l >= minl) {
              cout << vuidx << "\t"
                   << "INS"
                   << "\t" << l << "\t" << collp.second << "\t";
              cout << p->vertices[0];
              for (int i = 1; i < p->l; ++i)
                cout << ">" << p->vertices[i];
              cout << "\t" << 100 << "\t"
                   << "ACGT"[pseq[tp - 1]] << "\t";
              memcpy(ins, cons + qp, l);
              ins[l] = '\0';
              cout << decode(ins, l, 1) << "\t" << tp << "\t" << cigar << endl;
              ++vuidx;
            }
            qp += l;
          } else if (op == 2) {
            // D
            if (l >= minl) {
              cout << vuidx << "\t"
                   << "DEL"
                   << "\t" << l << "\t" << collp.second << "\t";
              cout << p->vertices[0];
              for (int i = 1; i < p->l; ++i)
                cout << ">" << p->vertices[i];
              memcpy(del, pseq + tp, l);
              del[l] = '\0';
              cout << "\t" << 100 << "\t" << decode(del, l, 1) << "\t";
              cout << "ACGT"[cons[qp - 1]] << "\t" << tp << "\t" << cigar
                   << endl;
              ++vuidx;
            }
            tp += l;
          }
        }
        free(ez.cigar);
      }
    }
    for (path_t *p : subpaths) {
      if (p != NULL)
        destroy_path(p);
    }
  }

  fprintf(stderr, "[M::%s] called %d variations in %.3f sec\n", __func__, vuidx,
          realtime() - rt);
  rt = realtime();
  // ---

  // Cleaning up
  for (auto &c : Cs) {
    for (const sfs_t s : c.specifics) {
      if (s.seq != NULL)
        free(s.seq);
    }
  }
  abpoa_free(ab);
  abpoa_free_para(abpt);
  free(ins);
  free(del);
  free(cigar);
  free(cons);
  free(pseq);

  gsk.destroy_graph();
  fprintf(stderr, "[M::%s] completed in %.3f sec\n", __func__,
          realtime() - rt0);
  // ---

  return 0;
}
