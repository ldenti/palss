#include <algorithm>
#include <cstdlib>
#include <getopt.h>
#include <map>
#include <stdio.h>
#include <string>
#include <zlib.h>

#include "kseq.h"

extern "C" {
#include "fm-index.h"
#include "kmer.h"
#include "sfs.h"
#include "sketch.h"
}
#include "graph.hpp"
#include "misc.hpp"
#include "usage.hpp"

KSEQ_INIT(gzFile, gzread)

typedef struct {
  std::string idx;
  std::string seq;
} read_t;

// Compute SFS strings from P and store them into solutions
std::vector<sfs_t> ping_pong_search(const rb3_fmi_t *index, uint8_t *P, int l,
                                    int qidx) {
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

    out.push_back({qidx, begin, end - begin + 1});

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
      S[new_n].qidx = S[i].qidx;
      S[new_n].s = S[i].s;
      S[new_n].l = S[i].l;
      ++new_n;
    }
    ++i;
  }

  // Resise to new size
  S.resize(new_n);
}

// Anchor specific strings on graph using graph sketch
void anchor(sketch_t *sketch, Graph &graph, std::vector<sfs_t> &SS, char *Q,
            int ql, int klen, int NA) {

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
    ckmer_d = MIN(kmer_d, rckmer_d);

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
      ckmer_d = MIN(kmer_d, rckmer_d);
    }

    memcpy(kmer, Q + end, klen);
    kmer_d = k2d(kmer, klen);
    rckmer_d = rc(kmer_d, klen);
    ckmer_d = MIN(kmer_d, rckmer_d);
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
      ckmer_d = MIN(kmer_d, rckmer_d);
    }

    if (nsa == 0 || nea == 0) {
      s.qidx = -1; // tag as invalid
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

    std::map<std::pair<int, int>, int> memo;
    std::map<std::pair<int, int>, int>::iterator hhit;
    int x = 0, y = 0;
    int xoff = 0, yoff = 0;
    std::pair<int, int> xy = {x, y};
    for (int i = 0; i < nsa; ++i) {
      x = sanchors[i].v;
      if (x == -1)
        // anchor has been filtered out since it was repeated in the read
        continue;
      xoff = sanchors[i].offset;
      xy.first = x;
      for (int j = 0; j < nea; ++j) {
        y = eanchors[j].v;
        if (y == -1)
          // anchor has been filtered out since it was repeated in the read
          continue;
        yoff = eanchors[j].offset;
        xy.second = y;
        if ((hhit = memo.find(xy)) == memo.end()) {
          memo[xy] = graph.distance(sanchors[i].v, eanchors[j].v);
          memo[std::make_pair(y, x)] = memo[xy];
        }
        dist = memo[xy];
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
      s.qidx = -1; // tag as invalid
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
  int d = 0;        // merge specific strings this close
  int bsize = 1000; // batch size
  // int hd = 0;  // hamming distance for fixing anchors
  int NA = 20; // number of kmers to check for anchoring
  int nth = 4; // number of threads
  int _c;
  while ((_c = getopt(argc, argv, "k:a:b:@:h")) != -1) {
    switch (_c) {
    case 'k':
      klen = std::stoi(optarg);
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
  if (argc - optind != 3) {
    fprintf(stderr, "%s", SEARCH_USAGE_MESSAGE);
    return 1;
  }
  std::string gbz_fn = argv[optind++];
  std::string fmd_fn = argv[optind++];
  std::string fq_fn = argv[optind++];

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

  // Graph sketching and extraction
  sketch_t *sketch = sk_load((gbz_fn + ".skt").c_str());
  fprintf(stderr, "[M::%s] loaded %ld sketches in %.3f sec\n", __func__,
          sketch->n, realtime() - rt);
  // Graph
  Graph graph(gbz_fn);
  graph.load();
  // graph.print_stats();

  // Specific strings computation and anchoring

  std::vector<read_t> entries(bsize);
  std::vector<std::vector<sfs_t>> output(bsize);
  fprintf(stderr, "[M::%s] pre-allocation in %.3f sec\n", __func__,
          realtime() - rt);
  rt = realtime();

  gzFile fp = gzopen(fq_fn.c_str(), "r");
  kseq_t *seq = kseq_init(fp);

  // Some statistics
  int anchored_n = 0;
  int unanchored_n = 0;
  int assembled_n = 0;

  int x;
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

      output[qq] =
          ping_pong_search(&fmd, seq, l, qq); // query idx is local to the batch

      // if (n == 0)
      //   continue;
      assemble(output[qq], d);
      anchor(sketch, graph, output[qq], entries[qq].seq.data(), l, klen, NA);

      // strand stuff
      int strands[2] = {0, 0};
      for (uint ss = 0; ss < output[qq].size(); ++ss) {
        if (output[qq][ss].qidx != -1)
          ++strands[output[qq][ss].strand];
      }
      // + strand if tie
      int strand = strands[0] > strands[1] ? 0 : 1;
      for (uint ss = 0; ss < output[qq].size(); ++ss) {
        if (output[qq][ss].qidx == -1)
          continue;
        if (output[qq][ss].strand == 0)
          // reverse
          output[qq][ss].s = l - (output[qq][ss].s + output[qq][ss].l);
        if (output[qq][ss].strand != strand)
          output[qq][ss].strand = 2;
      }
    }
    fprintf(stderr, "[M::%s] searched in %.3f sec\n", __func__,
            realtime() - rt1);
    rt1 = realtime();

    // output
    // TODO: use a thread to do this
    for (int qq = 0; qq < x; ++qq) {
      assembled_n += output[qq].size();
      for (uint j = 0; j < output[qq].size(); ++j) {
        sfs_t s = output[qq][j];
        if (s.qidx == -1) {
          ++unanchored_n;
          printf("X %d %s %d %d . . .:.:.:. .:.:.:.\n", qq,
                 entries[qq].idx.c_str(), s.s, s.l);
        } else {
          ++anchored_n;
          char t = s.strand == 2 ? 'S' : 'O';
          printf("%c %d %s %d %d %d %s %ld:%d:%ld %ld:%d:%ld %d %d\n", t,
                 s.qidx, entries[qq].idx.c_str(), s.s, s.l, s.strand, ".",
                 s.a.v, s.a.offset, s.a.seq, s.b.v, s.b.offset, s.b.seq,
                 s.a.closest, s.b.closest);
        }
      }
    }

    fprintf(stderr, "[M::%s] output in %.3f sec\n", __func__, realtime() - rt1);
  }

  // At this point, specific strings are anchored. Anchors follow + strand
  // on graph. If read was on -, we have reversed the specific strings so
  // that everything is on + strand

  fprintf(stderr,
          "[M::%s] Computed %d specific strings (%d anchored, %d unanchored) "
          "in %.3f sec\n",
          __func__, assembled_n, anchored_n, unanchored_n, realtime() - rt);

  kseq_destroy(seq);
  gzclose(fp);
  sk_destroy(sketch);
  rb3_fmi_free(&fmd);

  return 0;
}
