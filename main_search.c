#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "fm-index.h"
#include "ketopt.h"
#include "kseq.h"

#include "graph.h"
#include "kmer.h"
#include "misc.h"
#include "sfs.h"
#include "sketch.h"
#include "usage.h"

// KSEQ_INIT(gzFile, gzread) // we already init kstream in graph.h
// XXX: there should be a better way to do this
__KSEQ_TYPE(gzFile)
__KSEQ_BASIC(static, gzFile)
__KSEQ_READ(static)

/* Compute SFS strings from P and store them into solutions */
int ping_pong_search(const rb3_fmi_t *index, uint8_t *P, int l, int qidx,
                     sfs_t **output) {
  int n = 0;
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
    if (n == 8192) {
      fprintf(stderr, "We have to many sfs. Please reallocate\n");
      exit(EXIT_FAILURE);
    }
    output[n]->qidx = qidx;
    output[n]->s = begin;
    output[n]->l = end - begin + 1;
    ++n;

    if (begin == 0)
      break;
    //   if (config->overlap == 0) // Relaxed
    //     begin -= 1;
    //   else
    begin = end - 1;
  }
  return n;
}

/* Merge specifics strings that are too close (d-bp apart) on read
 */
int assemble(sfs_t **S, int n, int d) {
  /* Reverse the vector */
  sfs_t *s;
  int i = 0;
  for (i = 0; i < n >> 1; ++i) {
    s = S[n - 1 - i];
    S[n - 1 - i] = S[i];
    S[i] = s;
  }

  i = 0;
  while (i < n) {
    int j;
    for (j = i + 1; j <= n; ++j) {
      if (j == n || S[j - 1]->s + S[j - 1]->l <= S[j]->s - d) {
        /* non-overlapping: update first, clean others */
        S[i]->l = S[j - 1]->s + S[j - 1]->l - S[i]->s;
        for (int j2 = i + 1; j2 < j; ++j2)
          S[j2]->l = 0;
        break;
      }
    }
    i = j;
  }

  /* Remove gaps by shifting left */
  int new_n = 0;
  i = 0;
  while (i < n) {
    if (S[i]->l > 0) {
      S[new_n]->qidx = S[i]->qidx;
      S[new_n]->s = S[i]->s;
      S[new_n]->l = S[i]->l;
      ++new_n;
    }
    ++i;
  }

  return new_n;
}

/* Anchor specific strings on graph using graph sketch */
void anchor(sketch_t *sketch, graph_t *graph, sfs_t **SS, int nSS, char *Q,
            int ql, int klen, int NA) {

  int beg, end;
  /* char *kmer_s = malloc(klen + 1); */
  char *kmer = malloc(klen);
  uint64_t kmer_d;       // kmer
  uint64_t rckmer_d = 0; // reverse and complemented kmer
  uint64_t ckmer_d = 0;  // canonical kmer
  int c;                 // current char
  hit_t vx;              // hit from sketch

  sfs_t *s;
  for (int sidx = 0; sidx < nSS; ++sidx) {
    s = SS[sidx];

    /* Finding anchors in flanking regions */
    /* XXX: Do we want anchors overlapping the string? */
    beg = s->s - klen;
    beg = beg < 0 ? 0 : beg;
    end = s->s + s->l;
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
      s->qidx = -1; // tag as invalid
      continue;
    }

    /* keep only unique anchors on read */
    for (int i1 = 0; i1 < nsa; ++i1) {
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
    }

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

    /* Finding best pair of anchors */
    int mind = 100, d;
    int sax = -1, eax = -1; // index for selected anchors
    /* map<pair<int, int>, int> memo; */
    /* map<pair<int, int>, int>::iterator hhit; */
    /* int comp; */
    int x = 0, y = 0;
    int xoff = 0, yoff = 0;
    /* pair<int, int> xy = {x, y}; */
    for (uint i = 0; i < nsa; ++i) {
      x = sanchors[i].v;
      if (x == -1)
        // anchor has been filtered out since it was repeated in the read
        continue;
      xoff = sanchors[i].offset;
      /* xy.first = x; */
      for (int j = 0; j < nea; ++j) {
        y = eanchors[j].v;
        if (y == -1)
          // anchor has been filtered out since it was repeated in the read
          continue;
        yoff = eanchors[j].offset;
        /* xy.second = y; */
        /* if ((hhit = memo.find(xy)) == memo.end()) { */
        /*   memo[xy] = compatible(graph, sanchors[i].v, eanchors[j].v); */
        /*   // memo[make_pair(y, x)] = memo[xy]; */
        /* } */
        /* comp = memo[xy]; */
        /* if (!comp) */
        /*   continue; */
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
      /* fprintf(stderr, "F2\n"); */
      s->qidx = -1; // tag as invalid
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
    s->s = b;
    s->l = l;
    {
      // memcpy(s->a, &sa, sizeof(anchor_t));
      s->a->v = sa.v;
      s->a->offset = sa.offset;
      s->a->p = sa.p;
      s->a->seq = sa.seq;
    }
    {
      // memcpy(s->v, &ea, sizeof(anchor_t));
      s->b->v = ea.v;
      s->b->offset = ea.offset;
      s->b->p = ea.p;
      s->b->seq = ea.seq;
    }
    assert(sa.v == s->a->v);
    assert(ea.v == s->b->v);

    s->strand = strand;
    assert(sa.v <= ea.v);
  }
  free(kmer);
}

typedef struct {
  char *idx;
  char *seq;
  int l; // seq length
  int m; // allocated space for seq
} read_t;

// XXX: hardcoded
// FIXME: no overflow checks
read_t *init_read() {
  read_t *r = malloc(sizeof(read_t));
  r->idx = malloc(1024);
  r->seq = malloc(32768);
  r->m = 32768;
  return r;
}

void destroy_read(read_t *r) {
  free(r->idx);
  free(r->seq);
  free(r);
}

int load_batch(kseq_t *seq, read_t **entries, int nb) {
  int i = 0;
  int l = 0;
  read_t *r;
  while (i < nb && (l = kseq_read(seq)) >= 0) {
    r = entries[i];
    strncpy(r->idx, seq->name.s, seq->name.l);
    r->idx[seq->name.l] = '\0';
    strncpy(r->seq, seq->seq.s, seq->seq.l);
    r->seq[seq->seq.l] = '\0';
    r->l = seq->seq.l;
    ++i;
  }
  return i;
}

int main_search(int argc, char *argv[]) {
  double rt0 = realtime();
  double rt = rt0, rt1;

  int klen = 27;
  int d = 0;        // merge specific strings this close
  int bsize = 1000; // batch size
  // int hd = 0;  // hamming distance for fixing anchors
  int NA = 20; // number of kmers to check for anchoring
  int nth = 4; // number of threads
  static ko_longopt_t longopts[] = {{NULL, 0, 0}};
  ketopt_t opt = KETOPT_INIT;
  int _c;
  while ((_c = ketopt(&opt, argc, argv, 1, "k:a:b:@:h", longopts)) >= 0) {
    if (_c == 'k')
      klen = atoi(opt.arg);
    else if (_c == 'b')
      bsize = atoi(opt.arg);
    else if (_c == 'a')
      NA = atoi(opt.arg);
    else if (_c == '@')
      nth = atoi(opt.arg);
    else if (_c == 'h') {
      fprintf(stderr, "%s", SEARCH_USAGE_MESSAGE);
      return 0;
    }
  }
  if (argc - opt.ind != 4) {
    fprintf(stderr, "%s", SEARCH_USAGE_MESSAGE);
    return 1;
  }
  char *gfa_fn = argv[opt.ind++];
  char *skt_fn = argv[opt.ind++];
  char *fmd_fn = argv[opt.ind++];
  char *fq_fn = argv[opt.ind++];

  /* FMD-index loading */
  rb3_fmi_t fmd;
  rb3_fmi_restore(&fmd, fmd_fn, 0);
  if (fmd.e == 0 && fmd.r == 0) {
    fprintf(stderr, "[E::%s] FMD-index cannot be restored", __func__);
    exit(EXIT_FAILURE);
  }
  fprintf(stderr, "[M::%s] restored FMD index in %.3f sec\n", __func__,
          realtime() - rt);
  rt = realtime();

  /* Graph sketching and path extraction */
  sketch_t *sketch = sk_load(skt_fn);
  fprintf(stderr, "[M::%s] loaded %ld sketches in %.3f sec\n", __func__,
          sketch->n, realtime() - rt);
  rt = realtime();

  graph_t *graph = init_graph(gfa_fn);
  load_vertices(graph);
  load_paths(graph);

  fprintf(stderr, "[M::%s] loaded %d vertices and %d paths in %.3f sec\n",
          __func__, graph->nv, graph->np, realtime() - rt);
  rt = realtime();
  rt1 = rt;

  /* Specific strings computation and anchoring */

  read_t **entries = malloc(bsize * sizeof(read_t *));
  for (int i = 0; i < bsize; ++i)
    entries[i] = init_read();
  fprintf(stderr, "[M::%s] input allocation in %.3f sec\n", __func__,
          realtime() - rt);
  rt = realtime();

  int *output_n = malloc(bsize * sizeof(int));
  for (int i = 0; i < bsize; ++i)
    output_n[i] = 0;
  sfs_t ***output = malloc(bsize * sizeof(sfs_t **));
  for (int i = 0; i < bsize; ++i) {
    output[i] = malloc(8192 * sizeof(sfs_t *)); // XXX: hardcoded
    for (int j = 0; j < 8192; ++j)
      output[i][j] = init_sfs();
  }
  fprintf(stderr, "[M::%s] output allocation in %.3f sec\n", __func__,
          realtime() - rt);
  rt = realtime();

  gzFile fp = gzopen(fq_fn, "r");
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
      uint8_t *seq = (uint8_t *)entries[qq]->seq;
      sfs_t **out = output[qq];
      int l = entries[qq]->l;
      rb3_char2nt6(l, seq);

      int n = ping_pong_search(&fmd, seq, l, qq,
                               out); // query idx is local to the batch

      /* if (n == 0) */
      /*   continue; */
      n = assemble(out, n, d);
      output_n[qq] = n;

      anchor(sketch, graph, out, n, entries[qq]->seq, l, klen, NA);

      // strand stuff
      int strands[2] = {0, 0};
      for (int ss = 0; ss < n; ++ss) {
        if (out[ss]->qidx != -1)
          ++strands[out[ss]->strand];
      }
      // + strand if tie
      int strand = strands[0] > strands[1] ? 0 : 1;
      for (int ss = 0; ss < n; ++ss) {
        if (out[ss]->qidx == -1)
          continue;
        if (out[ss]->strand == 0)
          // reverse
          out[ss]->s = l - (out[ss]->s + out[ss]->l);
        if (out[ss]->strand != strand)
          out[ss]->strand = 2;
      }
    }
    fprintf(stderr, "[M::%s] searched in %.3f sec\n", __func__,
            realtime() - rt1);
    rt1 = realtime();

    // output
    // TODO: use a thread to do this
    for (int qq = 0; qq < x; ++qq) {
      assembled_n += output_n[qq];
      sfs_t *s;
      for (int j = 0; j < output_n[qq]; ++j) {
        s = output[qq][j];
        if (s->qidx == -1) {
          ++unanchored_n;
          printf("X %d %s %d %d . . .:.:. .:.:.\n", qq, entries[qq]->idx, s->s,
                 s->l);
        } else {
          ++anchored_n;
          char t = s->strand == 2 ? 'S' : 'O';
          printf("%c %d %s %d %d %d %s %ld:%d:%ld %ld:%d:%ld\n", t, s->qidx,
                 entries[qq]->idx, s->s, s->l, s->strand, ".", s->a->v,
                 s->a->offset, s->a->seq, s->b->v, s->b->offset, s->b->seq);
        }
      }
    }

    fprintf(stderr, "[M::%s] output in %.3f sec\n", __func__, realtime() - rt1);
  }

  /* At this point, specific strings are anchored. Anchors follow + strand on
   * graph. If read was on -, we have reversed the specific strings so that
   * everything is on + strand */

  fprintf(stderr,
          "[M::%s] Computed %d specific strings (%d anchored, %d unanchored) "
          "in %.3f sec\n",
          __func__, assembled_n, anchored_n, unanchored_n, realtime() - rt);

  for (int i = 0; i < bsize; ++i)
    destroy_read(entries[i]);
  free(entries);

  free(output_n);

  for (int i = 0; i < bsize; ++i) {
    for (int j = 0; j < 8192; ++j) {
      destroy_sfs(output[i][j]);
    }
    free(output[i]);
  }
  free(output);

  kseq_destroy(seq);
  gzclose(fp);
  sk_destroy(sketch);
  rb3_fmi_free(&fmd);
  destroy_graph(graph);

  fprintf(stderr, "[M::%s] done in %.3f sec\n", __func__, realtime() - rt0);

  return 0;
}
