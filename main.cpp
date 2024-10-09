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

// KSEQ_INIT(gzFile, gzread) // we already init kstream in gsketch
__KSEQ_TYPE(gzFile)
__KSEQ_BASIC(static, gzFile)
__KSEQ_READ(static)

using namespace std;

struct anchor_t {
  int v = -1;        // vertex on graph
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
  uint64_t ka = -1, kb = -1; // starting and ending kmers
};

int build_consensus(const vector<sfs_t> &specifics, char **cons, int *cons_c) {
  // TODO: move this outside and init just once
  // INIT ABPOA
  abpoa_t *ab = abpoa_init();
  abpoa_para_t *abpt = abpoa_init_para();
  abpt->disable_seeding = 0;
  abpt->align_mode = 0; // global
  abpt->out_msa = 0;
  abpt->out_cons = 1;
  abpt->out_gfa = 0;
  // abpt->is_diploid = 1; // TODO: maybe this works now
  abpt->progressive_poa = 1;
  abpt->amb_strand = 0;
  abpoa_post_set_para(abpt);
  // abpt->match = 2;      // match score
  // abpt->mismatch = 4;   // mismatch penalty
  // abpt->gap_mode = ABPOA_CONVEX_GAP; // gap penalty mode
  // abpt->gap_open1 = 4;  // gap open penalty #1
  // abpt->gap_ext1 = 2;   // gap extension penalty #1
  // abpt->gap_open2 = 24; // gap open penalty #2
  // abpt->gap_ext2 = 1;   // gap extension penalty #2
  // gap_penalty = min{gap_open1 + gap_len*gap_ext1, gap_open2+gap_len*gap_ext2}

  int *seq_lens = (int *)malloc(sizeof(int) * specifics.size());
  uint8_t **bseqs = (uint8_t **)malloc(sizeof(uint8_t *) * specifics.size());
  int goods = 0, i = 0;
  for (const sfs_t &s : specifics) {
    if (!s.good)
      continue;
    seq_lens[goods] = s.l;
    bseqs[goods] = (uint8_t *)malloc(sizeof(uint8_t) * (s.l + 1));
    for (i = 0; i < s.l; ++i)
      bseqs[goods][i] = s.seq[i] - 1;
    bseqs[goods][s.l] = '\0';

    if (!s.strand) {
      // rc
      for (i = 0; i < (s.l >> 1); ++i) {
        int tmp = bseqs[goods][s.l - 1 - i];
        tmp = (tmp >= 0 && tmp <= 3) ? 3 - tmp : tmp;
        bseqs[goods][s.l - 1 - i] =
            (bseqs[goods][i] >= 0 && bseqs[goods][i] <= 3) ? 3 - bseqs[goods][i]
                                                           : bseqs[goods][i];
        bseqs[goods][i] = tmp;
      }
      if (s.l & 1)
        bseqs[goods][i] = (bseqs[goods][i] >= 0 && bseqs[goods][i] <= 3)
                              ? 3 - bseqs[goods][i]
                              : bseqs[goods][i];
    }
    ++goods;
  }

  abpoa_msa(ab, abpt, goods, NULL, seq_lens, bseqs, NULL, NULL);
  abpoa_cons_t *abc = ab->abc;
  int cons_l = 0;
  if (abc->n_cons > 0) {
    cons_l = abc->cons_len[0];
    if (cons_l + 1 > *cons_c) {
      cerr << "--- Reallocating from " << *cons_c << " to " << cons_l + 1
           << endl;
      char *temp = (char *)realloc(*cons, (cons_l + 1) * sizeof(char));
      if (temp == NULL) {
        free(cons);
        cerr << "Error while reallocating memory for consensus string" << endl;
        exit(2);
      } else {
        *cons = temp;
      }
      *cons_c = cons_l + 1;
    }
    for (i = 0; i < cons_l; ++i)
      (*cons)[i] = abc->cons_base[0][i]; // "ACGTN"[abc->cons_base[0][i]];
    (*cons)[i] = '\0';
  }

  for (i = 0; i < goods; ++i)
    free(bseqs[i]);
  free(bseqs);
  free(seq_lens);
  abpoa_free(ab);
  abpoa_free_para(abpt);

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

    memcpy(kmer, P + b, k);
    kmer_d = k2d(kmer, k);
    rckmer_d = rc(kmer_d, k);
    ckmer_d = std::min(kmer_d, rckmer_d);
    vector<anchor_t> sanchors;
    while (b > 0 && sanchors.size() < N) {
      if ((vx = gsk.get(ckmer_d)) != -1)
        sanchors.push_back({vx, b, ckmer_d});
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
      if ((vx = gsk.get(ckmer_d)) != -1)
        eanchors.push_back({vx, e, ckmer_d});
      ++e;
      c = P[e + k - 1] < 5 ? P[e + k - 1] - 1 : rand() % 4;
      kmer_d = lsappend(kmer_d, c, k);
      rckmer_d = rsprepend(rckmer_d, reverse_char(c), k);
      ckmer_d = std::min(kmer_d, rckmer_d);
    }

    if (sanchors.size() == 0 || eanchors.size() == 0)
      // extended_sfs.push_back({b, e + k + 1 - b, -1, -1, 0});
      continue;

    /*for (const auto &sa:sanchors) {*/
    /*  cerr << sa.v << "|";*/
    /*}*/
    /*cerr << " ";*/
    /*for (const auto &ea:eanchors) {*/
    /*  cerr << ea.v << "|";*/
    /*}*/
    /*cerr << endl;*/

    int mind = 100;
    int d;
    int sax = -1, eax = -1; // index for selected anchors
    map<pair<int, int>, int> memo;
    map<pair<int, int>, int>::iterator hhit;
    int comp;
    int x = 0, y = 0;
    pair<int, int> xy = {x, y};
    for (int i = 0; i < sanchors.size(); ++i) {
      x = sanchors[i].v;
      xy.first = x;
      for (int j = 0; j < eanchors.size(); ++j) {
        y = eanchors[j].v;
        xy.second = y;
        cerr << "Checking " << x << ">" << y << endl;
        if ((hhit = memo.find(xy)) == memo.end()) {
          memo[xy] = gsk.compatible(sanchors[i].v, eanchors[j].v);
          cerr << "Computing from graph: " << memo[xy] << endl;
          /*memo[make_pair(y, x)] = memo[xy];*/
        }
        comp = memo[xy];
        cerr << "Comp: " << comp << endl;
        if (!comp)
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
    if (sa.v > ea.v) {
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
    c.ka = s.a.seq;
  }
  if (c.vb == -1 || s.b.v > c.vb) {
    c.vb = s.b.v;
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

string decode(const char *s, int l, int shift) {
  if (s == NULL)
    return "";
  char ds[l + 1];
  for (int i = 0; i < l; ++i)
    ds[i] = "NACGTN"[s[i] + shift];
  ds[l] = '\0';
  return ds;
}

int main(int argc, char *argv[]) {
  char *gfa_fn = argv[1];
  char *fmd_fn = argv[2];
  char *fq_fn = argv[3];
  int k = stoi(argv[4]);
  int w = 2;  // FIXME: hardcoded, minimum weight for clusters
  int hd = 0; // FIXME: hardcoded, hamming distance for fixing anchors

  cerr << "Building graph sketch..." << endl;
  GSK gsk(gfa_fn);
  gsk.build_sketch(k);
  gsk.build_graph();

  cerr << "Restoring FMD index..." << endl;
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
  uint8_t *s;
  vector<sfs_t> S;
  vector<sfs_t> SS;
  uint qidx = 0;
  vector<string> qnames;
  vector<int> strands(2);
  int strand;
  int n = 0;
  while ((l = kseq_read(seq)) >= 0) {
    cout << seq->name.s << endl;
    s = (uint8_t *)seq->seq.s;
    rb3_char2nt6(seq->seq.l, s);

    S = ping_pong_search(&f, s, qidx, seq->seq.l);
    S = assemble(S);
    cout << S.size() << " specific strings" << endl;
    S = anchor(S, s, l, gsk);
    cout << S.size() << " anchored specific strings" << endl;

    strands[0] = 0;
    strands[1] = 0;
    for (const auto &s : S)
      ++strands[s.strand];
    // FIXME: plus strand if tie
    cout << strands[0] << "/" << strands[1] << endl;
    strand = 1;
    if (strands[0] > strands[1])
      strand = 0;
    for (const auto &s : S)
      if (s.strand == strand)
        SS.push_back(s);

    qnames.push_back(seq->name.s);
    ++qidx;
  }
  for (const sfs_t &s : SS)
    assert(s.l > 0);
  kseq_destroy(seq);
  gzclose(fp);
  rb3_fmi_free(&f);

  cerr << "Sorting " << SS.size() << " specific strings..." << endl;
  std::sort(SS.begin(), SS.end(),
            [](const sfs_t &a, const sfs_t &b) { return a.a.v < b.b.v; });
  cerr << "Clustering..." << endl;
  vector<cluster_t> Cs = cluster(SS, gsk);

  cerr << "Merging " << Cs.size() << " clusters" << endl;
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
  cerr << lowsc_n << " clusters filtered" << endl;
  cerr << lowsc_am_n << " clusters filtered after merging" << endl;
  cerr << sc_n << " clusters will be analyzed" << endl;

  // checking if specifics start/end with "correct" kmers
  vector<vector<sfs_t *>> pspecifics(qnames.size());
  for (auto &c : Cs) {
    if (c.specifics.empty())
      continue;
    for (const sfs_t &s : c.specifics)
      assert(s.l > 0);
    for (auto &s : c.specifics) {
      if (s.a.seq == c.ka && s.b.seq == c.kb) {
        s.esk = c.ka;
        s.eek = c.kb;
      } else {
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
      if (s->a.seq != s->esk) {
        s->good = 2;
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
      if (s->b.seq != s->eek) {
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
          s->l += p + k - (s->s + s->l) + 1;
        }
      }

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

  char *pseq = (char *)malloc(4096 * sizeof(char));
  int pseq_c = 4096;
  int pseq_l = 0;
  char *cons = (char *)malloc(4096 * sizeof(char));
  int cons_c = 4096;
  int cons_l = 0;

  // From minimap2 asm5
  // https://github.com/lh3/minimap2/blob/69e36299168d739dded1c5549f662793af10da83/options.c#L141
  // int a, b, q, e, q2, e2; // matching score, mismatch, gap-open and gap-ext
  // penalties mo->a = 1, mo->b = 19, mo->q = 39, mo->q2 = 81, mo->e = 3, mo->e2
  // = 1, mo->zdrop = mo->zdrop_inv = 200;
  int sc_mch = 1, sc_mis = -9, gapo = 81, gape = 1;
  int8_t a = (int8_t)sc_mch,
         b = sc_mis < 0 ? (int8_t)sc_mis : -(int8_t)sc_mis; // a>0 and b<0
  int8_t mat[25] = {a, b, b, b, 0, b, a, b, b, 0, b, b, a,
                    b, 0, b, b, b, a, 0, 0, 0, 0, 0, 0};
  ksw_extz_t ez;

  for (auto &c : Cs) {
    if (c.specifics.empty())
      continue;

    cout << c.specifics.size() << " " << c.va << ">" << c.vb << " "
         << (c.va - 1) * 32 << "-" << (c.vb) * 32 << " " << d2s(c.ka, k) << " "
         << d2s(c.kb, k) << endl;
    vector<path_t *> subpaths = gsk.get_subpaths(c.va, c.vb);
    for (auto &s : c.specifics) {
      /*  assert((s.good > 0 && s.seq != NULL) || (s.good == 0 && s.seq ==
       * NULL));*/
      cout << s.good << " " << qnames[s.qidx] << ":" << s.s << "-" << s.s + s.l
           << " (" << s.l << ") " << s.strand << " " << s.a.v << ">" << s.b.v
           << " " << d2s(s.a.seq, k) << " " << d2s(s.b.seq, k) << " "
           << decode(s.seq, s.l, 0) << endl;
    }

    cons_l = build_consensus(c.specifics, &cons, &cons_c);
    cout << decode(cons, cons_l, 1) << endl;
    cout << "Analyzing " << subpaths.size() << " subpaths" << endl;
    for (const path_t *p : subpaths) {
      pseq_l = gsk.get_sequence(p, &pseq, &pseq_c);
      cout << decode(pseq, pseq_l, 1) << endl;

      memset(&ez, 0, sizeof(ksw_extz_t));
      ksw_extz2_sse(0, cons_l, (uint8_t *)cons, pseq_l, (uint8_t *)pseq, 5, mat,
                    4, 2, -1, 200, 0, 0, &ez);
      for (int i = 0; i < ez.n_cigar; ++i)
        printf("%d%c", ez.cigar[i] >> 4, "MID"[ez.cigar[i] & 0xf]);
      putchar('\n');
      free(ez.cigar);
    }
    cout << endl;
  }
  free(cons);
  free(pseq);

  gsk.destroy_graph();
  cerr << "END" << endl;
  return 0;
}
