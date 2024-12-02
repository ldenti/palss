#include <assert.h>
#include <vector>
#include <zlib.h>

#include "fm-index.h"
#include "ketopt.h"
#include "kseq.h"

extern "C" {
#include "graph.h"
#include "sfs.h"
}
#include "sketch.hpp"
#include "utils.h"

// KSEQ_INIT(gzFile, gzread) // we already init kstream in graph.h
// XXX: there should be a better way to do this
__KSEQ_TYPE(gzFile)
__KSEQ_BASIC(static, gzFile)
__KSEQ_READ(static)

using namespace std;

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

/* Merge specifics strings that are too close (d-bp apart) on the same read */
vector<sfs_t> assemble(const vector<sfs_t> &sfs, int d) {
  vector<sfs_t> assembled_sfs;
  int i = sfs.size() - 1;
  while (i >= 0) {
    int j;
    for (j = i - 1; j >= 0; --j) {
      if (sfs[j + 1].s + sfs[j + 1].l <= sfs[j].s - d) {
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

/* Anchor specific strings on graph using graph sketch */
vector<sfs_t> anchor(const sketch_t &sketch, graph_t *graph,
                     const vector<sfs_t> &sfs, uint8_t *P, int l, int klen,
                     int N) {
  vector<sfs_t> anchored_sfs;
  int beg, end, ext = 0;
  pair<int64_t, uint16_t> vx = make_pair(-1, -1);
  char *kmer = (char *)malloc((klen + 1) * sizeof(char));
  kmer[klen] = '\0';
  uint64_t kmer_d;       // kmer
  uint64_t rckmer_d = 0; // reverse and complemented kmer
  uint64_t ckmer_d = 0;  // canonical kmer
  int c;                 // current char

  for (const sfs_t &s : sfs) {
    // Finding anchors in flanking regions
    // XXX: Do we want anchors overlapping the string?
    beg = s.s - klen;
    beg = beg < 0 ? 0 : beg;
    end = s.s + s.l;
    end = end > l - klen ? l - klen : end;

    memcpy(kmer, P + beg, klen);
    kmer_d = k2d(kmer, klen);
    rckmer_d = rc(kmer_d, klen);
    ckmer_d = std::min(kmer_d, rckmer_d);
    vector<anchor_t> sanchors;
    // XXX: do we want at least one anchor?
    while (beg > 0 && (sanchors.size() < N /*|| sanchors.empty()*/)) {
      if ((vx = sk_get(sketch, ckmer_d)).first != -1)
        sanchors.push_back({vx.first, vx.second, beg, ckmer_d});
      --beg;
      ++ext;
      c = P[beg] < 5 ? P[beg] - 1 : rand() % 4;
      kmer_d = rsprepend(kmer_d, c, klen);
      rckmer_d = lsappend(rckmer_d, reverse_char(c), klen);
      ckmer_d = std::min(kmer_d, rckmer_d);
    }
    ext = 0;
    memcpy(kmer, P + end, klen);
    kmer_d = k2d(kmer, klen);
    rckmer_d = rc(kmer_d, klen);
    ckmer_d = std::min(kmer_d, rckmer_d);
    vector<anchor_t> eanchors;
    // XXX: do we want at least one anchor?
    while (end < l - klen + 1 &&
           (eanchors.size() < N /*|| eanchors.empty()*/)) {
      if ((vx = sk_get(sketch, ckmer_d)).first != -1)
        eanchors.push_back({vx.first, vx.second, end, ckmer_d});
      ++end;
      ++ext;
      c = P[end + klen - 1] < 5 ? P[end + klen - 1] - 1 : rand() % 4;
      kmer_d = lsappend(kmer_d, c, klen);
      rckmer_d = rsprepend(rckmer_d, reverse_char(c), klen);
      ckmer_d = std::min(kmer_d, rckmer_d);
    }

    if (sanchors.size() == 0 || eanchors.size() == 0) {
      anchored_sfs.push_back({s.qidx, s.s, s.l, -1, -1, -1});
      // fprintf(stderr, "Lost (I) %s:%d-%d\n", qname, s.s, s.s+s.l);
      continue;
    }

    // Finding best pair of anchors
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
          memo[xy] = compatible(graph, sanchors[i].v, eanchors[j].v);
          // memo[make_pair(y, x)] = memo[xy];
        }
        comp = memo[xy];
        if (!comp)
          continue;
        if (x == y && (xoff == yoff || (xoff < yoff && xoff + klen >= yoff) ||
                       (xoff > yoff && yoff + klen >= xoff)))
          continue;
        d = abs(sanchors[i].v - eanchors[j].v);
        if (d < mind) {
          sax = i;
          eax = j;
          mind = d;
        }
      }
    }
    if (sax == -1 || eax == -1) {
      anchored_sfs.push_back({s.qidx, s.s, s.l, -1, -1, -1});
      // fprintf(stderr, "Lost (II) %s:%d-%d\n", qname, s.s, s.s+s.l);
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
    anchored_sfs.push_back({s.qidx, b, l, sa, ea, strand});
    assert(sa.v <= ea.v);
  }
  free(kmer);

  return anchored_sfs;
}

int main_search(int argc, char *argv[]) {
  double rt0 = realtime();
  double rt = rt0, rt1;

  int klen = 27; // kmer size
  int d = 100;   // merge specific strings this close
  int hd = 0;    // hamming distance for fixing anchors
  int N = 20;    // number of kmers to check for anchoring
  static ko_longopt_t longopts[] = {};
  ketopt_t opt = KETOPT_INIT;
  int _c;
  while ((_c = ketopt(&opt, argc, argv, 1, "k:d:w:d:l:a:r:", longopts)) >= 0) {
    if (_c == 'k')
      klen = atoi(opt.arg);
    else if (_c == 'd')
      d = atoi(opt.arg);
    else if (_c == 'd')
      hd = atoi(opt.arg);
    else if (_c == 'a')
      N = atoi(opt.arg);
  }
  if (argc - opt.ind != 4) {
    fprintf(stderr, "Argh");
    return 1;
  }
  char *gfa_fn = argv[opt.ind++];
  char *skt_fn = argv[opt.ind++];
  char *fmd_fn = argv[opt.ind++];
  char *fq_fn = argv[opt.ind++];

  // FMD-index loading
  rb3_fmi_t fmd;
  rb3_fmi_restore(&fmd, fmd_fn, 0);
  if (fmd.e == 0 && fmd.r == 0) {
    fprintf(stderr, "Error restoring index");
    return 1;
  }
  fprintf(stderr, "[M::%s] restored FMD index in %.3f sec\n", __func__,
          realtime() - rt);
  rt = realtime();
  // ---

  // Graph sketching and path extraction
  sketch_t sketch;
  sk_load(sketch, skt_fn);
  fprintf(stderr, "[M::%s] loaded %ld sketches in %.3f sec\n", __func__,
          sketch.size(), realtime() - rt);
  rt = realtime();

  graph_t *graph = init_graph(gfa_fn);
  load_paths(graph);

  fprintf(stderr, "[M::%s] loaded %d paths in %.3f sec\n", __func__, graph->np,
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
  uint qidx = 0;
  vector<int> strands(2);
  int strand;

  int specifics_n = 0;
  int anchored_n = 0;
  while ((l = kseq_read(seq)) >= 0) {
    s = (uint8_t *)seq->seq.s;
    rb3_char2nt6(seq->seq.l, s);

    S = ping_pong_search(&fmd, s, qidx, seq->seq.l);
    // strings are sorted wrt their position on read (reverse)
    S = assemble(S, d);
    // strings are now sorted wrt their position on read
    for (int i = 0; i < (int)S.size() - 1; ++i)
      assert(S[i].s < S[i + 1].s);
    specifics_n += S.size();

    S = anchor(sketch, graph, S, s, l, klen, N);

    strands[0] = 0;
    strands[1] = 0;
    for (const auto &s : S) {
      if (s.a.v != -1 && s.b.v != -1)
        ++strands[s.strand];
    }
    // + strand if tie
    strand = strands[0] > strands[1] ? 0 : 1;

    for (sfs_t &s : S) {
      // a.v always precedes b.v (even for sfs on - strand)
      printf("%s:%d-%d %d %d %d %d %ld:%d:%ld>%ld:%d:%ld\n", seq->name.s, s.s,
             s.s + s.l, s.strand ? s.s : l - (s.s + s.l), s.l, s.strand,
             s.strand == strand, s.a.v, s.a.offset, s.a.seq, s.b.v, s.b.offset,
             s.b.seq);
      if (s.strand == -1 || s.a.v == -1 || s.b.v == -1)
        continue;
      ++anchored_n;

      // we may want to "reverse" the sfs on the read if - strand
      // if (!s.strand)
      //   s.s = l - (s.s + s.l);
    }
    ++qidx;
    if (qidx % 10000 == 0) {
      fprintf(stderr, "[M::%s] parsed %d reads %.3f sec\n", __func__, qidx,
              realtime() - rt1);
      rt1 = realtime();
    }
  }
  fprintf(stderr, "%d specific strings, %d anchored strings (lost: %d)\n",
          specifics_n, anchored_n, specifics_n - anchored_n);

  kseq_destroy(seq);
  gzclose(fp);
  rb3_fmi_free(&fmd);
  destroy_graph(graph);
  fprintf(stderr, "[M::%s] done in %.3f sec\n", __func__, realtime() - rt0);

  return 0;
}
