#include <algorithm>
#include <assert.h>
#include <bit>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <zlib.h>

#include "abpoa.h"
#include "fm-index.h"
#include "ketopt.h"
#include "kseq.h"

#include "sfs.h"

#include "cluster.hpp"

extern "C" {
#include "graph.h"
#include "ksw2.h" // XXX: redefinitions. This has to be included after graph.h
#include "usage.h"
}
#include "sketch.hpp"
#include "utils.h"

// KSEQ_INIT(gzFile, gzread) // we already init kstream in graph.h
// XXX: there should be a better way to do this
__KSEQ_TYPE(gzFile)
__KSEQ_BASIC(static, gzFile)
__KSEQ_READ(static)

using namespace std;

string decode(const uint8_t *s, int l, int shift) {
  if (s == NULL)
    return "";
  char ds[l + 1];
  for (int i = 0; i < l; ++i)
    ds[i] = s[i] + shift <= 5 ? "NACGTN"[s[i] + shift] : 'N';
  ds[l] = '\0';
  return ds;
}

/* Add specific string s to cluster c and update its anchors */
void cl_insert(cluster_t &c, sfs_t *s) {
  c.specifics.push_back(s);
  if (c.va == -1 || s->a.v < c.va) {
    c.va = s->a.v;
    c.offa = s->a.offset;
    c.ka = s->a.seq;
  }
  if (c.vb == -1 || s->b.v >= c.vb) {
    if (s->b.v == c.vb) {
      if (s->b.offset > c.offb) {
        c.offb = s->b.offset;
        c.kb = s->b.seq;
      }
    } else {
      c.vb = s->b.v;
      c.offb = s->b.offset;
      c.kb = s->b.seq;
    }
  }
}

/* Sweep line clustering of specific strings based on their anchors. Assuming
 * topological ordering of vertices */
vector<cluster_t> cluster(const vector<sfs_t *> SS, int klen) {
  vector<cluster_t> clusters(1);
  cl_insert(clusters.back(), SS[0]);
  for (int i = 1; i < SS.size(); ++i) {
    if (SS[i]->a.v > clusters.back().vb ||
        (SS[i]->a.v == clusters.back().vb &&
         SS[i]->a.offset >= clusters.back().offb + klen)) {
      // no overlap, so new cluster
      clusters.push_back({});
    }
    cl_insert(clusters.back(), SS[i]);
  }
  return clusters;
}

int merge(cluster_t &C) {
  vector<sfs_t *> newC;
  map<int, vector<sfs_t *>> byread;
  for (sfs_t *s : C.specifics)
    byread[s->qidx].push_back(s);

  int skipped = 0;
  sfs_t *s;
  for (const auto &c : byread) {
    int mins = 100000, maxe = 0; // FIXME: assuming HiFi
    int first = -1, last = -1;
    for (int i = 0; i < c.second.size(); ++i) {
      if (c.second[i]->s < mins) {
        mins = c.second[i]->s;
        first = i;
      }
      if (c.second[i]->s + c.second[i]->l > maxe) {
        maxe = c.second[i]->s + c.second[i]->l;
        last = i;
      }
    }
    if (first == -1 || last == -1) {
      exit(1);
    }
    if (c.second[first]->a.v < c.second[last]->b.v ||
        (c.second[first]->a.v == c.second[last]->b.v &&
         c.second[first]->a.offset < c.second[last]->b.offset)) {

      s = init_sfs();
      s->qidx = c.first;
      s->s = c.second[first]->s;
      s->l = c.second[last]->s + c.second[last]->l - c.second[first]->s;
      s->a = c.second[first]->a;
      s->b = c.second[last]->b;
      s->strand = c.second[first]->strand;
      newC.push_back(s);

    } else {
      // TODO
      ++skipped;
    }
  }
  C.specifics = newC;
  return skipped;
}

vector<sfs_t *> load_sfs(char *fn) {
  vector<sfs_t *> SS;
  FILE *fp;
  char *line = NULL;
  size_t len = 0;
  ssize_t read;

  fp = fopen(fn, "r");
  if (fp == NULL)
    exit(EXIT_FAILURE);

  sfs_t *s;
  while ((read = getline(&line, &len, fp)) != -1) {
    if (line[0] != 'O')
      continue;
    s = init_sfs();
    parse_sfs_line(line + 2, s);
    if (s->b.v - s->a.v < 50) // FIXME
      SS.push_back(s);
    else
      destroy_sfs(s);
  }
  fclose(fp);
  free(line);
  return SS;
}

vector<vector<sfs_t *>> split_by_len(cluster_t &cluster, float lr) {
  vector<int> subclusters_avgl;
  vector<vector<sfs_t *>> subclusters;

  // move to first good specific string and create the first cluster
  int i = 0;
  while (i < cluster.specifics.size() && !cluster.specifics[i]->good)
    ++i;
  if (i == cluster.specifics.size())
    return subclusters;
  subclusters_avgl.push_back(cluster.specifics[i]->l);
  subclusters.push_back({{cluster.specifics[i]}});

  // iterate over remaining specific strings
  for (i = i + 1; i < cluster.specifics.size(); ++i) {
    if (!cluster.specifics[i]->good)
      continue;
    int j = 0;
    for (j = 0; j < subclusters.size(); ++j) {
      if (min(cluster.specifics[i]->l, subclusters_avgl[j]) /
              (float)max(cluster.specifics[i]->l, subclusters_avgl[j]) >=
          lr)
        break;
    }
    if (j == subclusters.size()) {
      subclusters.push_back({cluster.specifics[i]});
      subclusters_avgl.push_back(cluster.specifics[i]->l);
    } else {
      subclusters_avgl[j] =
          (subclusters_avgl[j] * subclusters.size() + cluster.specifics[i]->l) /
          (subclusters.size() + 1);
      subclusters[j].push_back(cluster.specifics[i]);
    }
  }

  // sort subclusters by size
  sort(subclusters.begin(), subclusters.end(),
       [](const vector<sfs_t *> &a, const vector<sfs_t *> &b) {
         return a.size() > b.size();
       });
  return subclusters;
}

int main_call(int argc, char *argv[]) {
  double rt0 = realtime();
  double rt = rt0, rt1;

  int klen = 27;   // kmer size
  int min_w = 2;   // minimum support for cluster
  float lr = 0.97; // Ratio to cluster specific strings inside same cluster
  int hd = 0;      // hamming distance for fixing anchors
  int verbose = 0;
  // char *bed_fn = NULL; // bed file

  static ko_longopt_t longopts[] = {/*{"bed", ko_required_argument, 301},*/
                                    {NULL, 0, 0}};
  ketopt_t opt = KETOPT_INIT;
  int _c;
  while ((_c = ketopt(&opt, argc, argv, 1, "k:w:l:r:vh", longopts)) >= 0) {
    if (_c == 'k')
      klen = atoi(opt.arg);
    else if (_c == 'w')
      min_w = atoi(opt.arg);
    else if (_c == 'r')
      lr = atof(opt.arg);
    else if (_c == 'v')
      verbose = 1;
    else if (_c == 'h') {
      fprintf(stderr, "%s", CALL_USAGE_MESSAGE);
      return 0;
    }
  }
  if (argc - opt.ind != 4) {
    fprintf(stderr, "%s", CALL_USAGE_MESSAGE);
    return 1;
  }
  char *gfa_fn = argv[opt.ind++];
  char *skt_fn = argv[opt.ind++];
  char *sfs_fn = argv[opt.ind++];
  char *fq_fn = argv[opt.ind++];

  // FILE *bed_f = NULL;
  // if (bed_fn != NULL)
  //   bed_f = fopen(bed_fn, "w");

  // Sketch loading
  sketch_t sketch;
  sk_load(sketch, skt_fn);
  fprintf(stderr, "[M::%s] loaded %ld sketches in %.3f sec\n", __func__,
          sketch.size(), realtime() - rt);
  rt = realtime();

  // Graph loading
  graph_t *graph = init_graph(gfa_fn);
  load_vertices(graph);
  load_paths(graph);
  fprintf(stderr, "[M::%s] loaded %d paths in %.3f sec\n", __func__, graph->np,
          realtime() - rt);
  rt = realtime();

  // Extracting vertex to reference position
  // (XXX: only to first path, assuming reference sequence)
  // map<int, int> positions;
  // path_t *path = graph->paths[0];
  // char *ref_name = path->idx;
  // int refpath_len = 0;
  // seg_t *seg;
  // for (int i = 0; i < path->l; ++i) {
  //   seg = graph->vertices[path->vertices[i]];
  //   positions[seg->idx] = refpath_len;
  //   refpath_len += seg->l;
  // }

  vector<sfs_t *> SS = load_sfs(sfs_fn);
  fprintf(stderr, "[M::%s] loaded %d specific strings in %.3f sec\n", __func__,
          SS.size(), realtime() - rt);
  rt = realtime();

  pair<int64_t, int16_t> p1, p2;
  uint64_t v1, v2;
  int pos1, pos2;

  int nreads = 0;
  for (sfs_t *ss : SS) {
    if (ss->qidx > nreads)
      nreads = ss->qidx;
    /* // uncomment if we want specific strings in the bed file
    if (bed_f != NULL) {
      p1 = sk_get(sketch, ss.a.seq);
      v1 = p1.first;
      p2 = sk_get(sketch, ss.b.seq);
      v2 = p2.first;

      if (v1 != -1 && v2 != -1) {
        // anchored specific string
        if (positions.find(v1) == positions.end() ||
            positions.find(v2) == positions.end()) {
          // TODO: find best path by intersecting paths passing trough the
          two
          // vertices and report over it
          // printf("%s\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n", ss.rname);
        } else {
          pos1 = positions.at(v1) + p1.second;
          pos2 = positions.at(v2) + p2.second + klen;
          // d = pos2 - (pos1 + klen);

          fprintf(bed_f, "%s\t%d\t%d\t%s:%d-%d\n", ref_name, pos1, pos2,
                  ss.rname, ss.s, ss.s + ss.l);
        }
      }
    }
    */
  }
  nreads += 1; // since qidx starts from 0

  std::sort(SS.begin(), SS.end(), [](const sfs_t *a, const sfs_t *b) {
    return a->a.v < b->a.v || (a->a.v == b->a.v && a->a.offset < b->a.offset);
  });
  fprintf(stderr, "[M::%s] sorted specific strings in %.3f sec\n", __func__,
          realtime() - rt);
  rt = realtime();

  // ---

  // Clustering
  vector<cluster_t> Cs = cluster(SS, klen);
  fprintf(stderr, "[M::%s] created %ld clusters in %.3f sec\n", __func__,
          Cs.size(), realtime() - rt);
  rt = realtime();
  rt1 = rt;
  // ---

  // Merging clusters
  int lowsc_n = 0;
  int lowsc_am_n = 0;
  int sc_n = 0;
  int total = 0;
  int skipped = 0;
  for (cluster_t &c : Cs) {
    if (c.specifics.size() < min_w) {
      ++lowsc_n;
      c.specifics.clear();
    } else {
      total += c.specifics.size();
      skipped += merge(c);
      if (c.specifics.size() < min_w) {
        ++lowsc_am_n;
      } else {
        ++sc_n;
      }
    }
  }
  for (sfs_t *ss : SS)
    destroy_sfs(ss);

  fprintf(stderr,
          "[M::%s] merged clusters in %.3f sec (%d+%d=%d clusters filtered - "
          "%d/%d strings skipped)\n",
          __func__, realtime() - rt, lowsc_n, lowsc_am_n, lowsc_n + lowsc_am_n,
          skipped, total);
  rt = realtime();

  /* We now need to clean specific strings to make them start/end with "correct"
   * kmers wrt cluster */

  // Split specific strings by read
  vector<vector<sfs_t *>> ss_byread(nreads);
  for (cluster_t &c : Cs) {
    if (c.specifics.size() < min_w)
      continue;
    for (sfs_t *s : c.specifics) {
      s->esk = c.ka;
      s->eek = c.kb;
      s->good = (s->a.seq == c.ka && s->b.seq == c.kb);
      ss_byread[s->qidx].push_back(s);
    }
  }

  gzFile fp = gzopen(fq_fn, "r");
  kseq_t *seq = kseq_init(fp);
  uint8_t *s;
  int l;
  int qidx = 0;

  char *kmer = (char *)malloc((klen + 1) * sizeof(char));
  kmer[klen] = '\0';
  uint64_t kmer_d;       // kmer
  uint64_t rckmer_d = 0; // reverse and complemented kmer
  uint64_t ckmer_d = 0;  // canonical kmer
  int c;                 // current char
  int p;
  int pc;
  vector<string> qnames;
  while ((l = kseq_read(seq)) >= 0) {
    qnames.push_back(seq->name.s);
    if (qidx >= ss_byread.size() || ss_byread[qidx].empty()) {
      ++qidx;
      continue;
    }

    s = (uint8_t *)seq->seq.s;
    rb3_char2nt6(seq->seq.l, s);
    if (ss_byread[qidx][0]->strand == 0) {
      // Reverse and complement the sequence
      int i;
      for (i = 0; i < (l >> 1); ++i) {
        int tmp = s[l - 1 - i];
        tmp = (tmp >= 1 && tmp <= 4) ? 5 - tmp : tmp;
        s[l - 1 - i] = (s[i] >= 1 && s[i] <= 4) ? 5 - s[i] : s[i];
        s[i] = tmp;
      }
      if (l & 1)
        s[i] = (s[i] >= 1 && s[i] <= 4) ? 5 - s[i] : s[i];
    }

    // Extend each specific string to match the expected anchors (from cluster)
    for (sfs_t *s : ss_byread[qidx]) {
      // left anchor
      if (s->a.seq != s->esk && s->s > 0) {
        p = s->s - 1;
        memcpy(kmer, seq->seq.s + p, klen);
        kmer_d = k2d(kmer, klen);
        rckmer_d = rc(kmer_d, klen);
        ckmer_d = std::min(kmer_d, rckmer_d);
        pc = hd + 1;
        while (p > 0 && (pc = popcount(ckmer_d ^ s->esk)) > hd) {
          --p;
          c = seq->seq.s[p] < 5 ? seq->seq.s[p] - 1 : rand() % 4;
          kmer_d = rsprepend(kmer_d, c, klen);
          rckmer_d = lsappend(rckmer_d, reverse_char(c), klen);
          ckmer_d = std::min(kmer_d, rckmer_d);
        }
        if (pc > hd) {
          s->good = 0;
        } else {
          s->a.seq = ckmer_d; // TODO: update also the vertex
          s->l += s->s - p;
          s->s = p;
        }
      }

      // right anchor
      if (s->b.seq != s->eek && s->s + s->l <= l - klen) {
        p = s->s + s->l;
        memcpy(kmer, seq->seq.s + p, klen);
        kmer_d = k2d(kmer, klen);
        rckmer_d = rc(kmer_d, klen);
        ckmer_d = std::min(kmer_d, rckmer_d);
        pc = hd + 1;
        while (p < l - klen + 1 && (pc = popcount(ckmer_d ^ s->eek)) > hd) {
          ++p;
          c = seq->seq.s[p + klen - 1] < 5 ? seq->seq.s[p + klen - 1] - 1
                                           : rand() % 4;
          kmer_d = lsappend(kmer_d, c, klen);
          rckmer_d = rsprepend(rckmer_d, reverse_char(c), klen);
          ckmer_d = std::min(kmer_d, rckmer_d);
        }
        if (pc > hd) {
          s->good = 0;
        } else {
          s->b.seq = ckmer_d; // XXX: update also the vertex
          s->l += p + klen - (s->s + s->l);
        }
      }
      s->good = s->a.seq == s->esk && s->b.seq == s->eek;

      if (s->good) {
        s->seq = (uint8_t *)malloc((s->l + 1) * sizeof(uint8_t));
        memcpy(s->seq, seq->seq.s + s->s, s->l);
        s->seq[s->l] = '\0';
        // XXX: convert to 0123 since it seems abpoa works better this way in
        // diploid mode
        for (int i = 0; i < s->l; ++i)
          --s->seq[i];
      }

      // if (bed_f != NULL) {
      //   p1 = sk_get(sketch, s->esk);
      //   v1 = p1.first;
      //   p2 = sk_get(sketch, s->eek);
      //   v2 = p2.first;
      //   pos1 = positions.at(v1) + p1.second;
      //   pos2 = positions.at(v2) + p2.second + klen;
      //   fprintf(bed_f, "%s\t%d\t%d\t%s:%d-%d:%d/2\n", ref_name, pos1, pos2,
      //           seq->name.s, s->s, s->s + s->l, s->good);
      // }
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

  // --- Analyzing clusters

  path_t *path = NULL;
  char *pathseq = (char *)malloc(8192 * sizeof(char));
  int pathseq_c = 8192;
  int pathseq_len = 0;

  int bestpath_idx;
  path_t *bestpath = NULL;
  uint8_t *bestpathseq = (uint8_t *)malloc(8192 * sizeof(uint8_t));
  int bestpathseq_c = 8192;
  int bestpathseq_len = 0;
  int bestpath_lastprefix;

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

  uint32_t *cigar = (uint32_t *)malloc(1 * sizeof(uint32_t));
  int n_cigar = 0;
  int m_cigar = 1;

  // === INIT ABPOA
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

  // XXX: Assuming clusters of size <= 128
  int *ssseqs_lens = (int *)malloc(sizeof(int) * 128);
  uint8_t **ssseqs = (uint8_t **)malloc(sizeof(uint8_t *) * 128);
  int goods = 0;

  // XXX
  char *cigar_s = (char *)malloc(16384 * 2 * sizeof(char)); // FIXME: hardcoded

  int cc_idx = 0;
  for (auto &cc : Cs) {
    if (cc.specifics.size() < min_w)
      continue;

    ++cc_idx;
    // if (cc_idx != 48552)
    //   continue;

    if (cc_idx % 5000 == 0) {
      fprintf(stderr, "[M::%s] analyzed %d clusters (%d left) in %.3f sec\n",
              __func__, cc_idx, sc_n - cc_idx, realtime() - rt1);
      rt1 = realtime();
    }

    int va = cc.va, vb = cc.vb;
    int offa = cc.offa, offb = cc.offb;

    // if (vb - va > 100)
    //   continue;
    if (verbose)
      fprintf(stderr, "\n\n%d : %ld %d>%d\n", cc_idx, cc.specifics.size(),
              cc.va, cc.vb);

    // split cluster based on specific strings length
    vector<vector<sfs_t *>> subclusters = split_by_len(cc, lr);
    if (subclusters.size() == 0)
      continue;
    int potential_diploid = subclusters.size() == 1;

    // ---

    vector<pair<path_t *, string>> collapsed_subpaths;
    for (int px = 0; px < graph->np; ++px) {
      path = extract_subpath(graph, graph->paths[px], va, vb);
      if (path == NULL)
        continue;
      int i = 0;
      for (const auto &csp : collapsed_subpaths) {
        if (path->l != csp.first->l)
          continue;
        int j = 0;
        for (j = 0; j < path->l; ++j) {
          if (path->vertices[j] != csp.first->vertices[j])
            break;
        }
        if (j == path->l)
          break;
        ++i;
      }
      if (i == collapsed_subpaths.size()) {
        collapsed_subpaths.push_back(make_pair(path, path->idx));
      } else {
        collapsed_subpaths[i].second += "," + string(path->idx);
      }
    }
    if (path == NULL) {
      if (verbose)
        fprintf(stderr, "No path for cluster %d (%d>%d)\n", cc_idx, va, vb);
      continue;
    }

    // ---

    for (int sci = 0; sci < (potential_diploid ? 1 : 2); ++sci) {
      vector<sfs_t *> specifics = subclusters[sci];
      if (specifics.size() < min_w)
        continue;

      // abpoa
      goods = 0;
      for (const sfs_t *s : specifics) {
        if (s->good == 0)
          continue;
        assert(s->l > 0);
        ssseqs_lens[goods] = s->l;
        ssseqs[goods] = s->seq;
        // printf(">%d.%d\n%s\n", goods, s->l, decode(s->seq, s->l, 1).c_str());
        ++goods;
      }
      assert(goods < 128);
      // printf(">>> %d\n", goods);
      abpoa_msa(ab, abpt, goods, NULL, ssseqs_lens, ssseqs, NULL, NULL);

      abpoa_cons_t *abc = ab->abc;
      for (int ci = 0; ci < (potential_diploid ? abc->n_cons : 1); ++ci) {
        int cons_len = abc->cons_len[ci];
        uint8_t *cons_seq = abc->cons_base[ci];
        int score = -INT32_MAX;
        int abpoa_supp = abc->clu_n_seq[ci];
        if (abpoa_supp < min_w)
          continue;
        // printf(">>> %d.%d.%d - %d - %d : %d %d\n", cc_idx, sci, ci,
        //        specifics.size(), potential_diploid,
        //        potential_diploid ? abc->n_cons : 1, abpoa_supp);

        // uint8_t *PPSEQ = NULL;
        // int ppseq_c = 0;
        // string pnames;
        // int plen = 0;
        //       int pprefix = 0;
        //       int PPSEQ_IDX = -1;

        int path_idx = -1;

        for (pair<path_t *, string> collp : collapsed_subpaths) {
          ++path_idx;
          path = collp.first;
          pathseq_len = get_sequence(graph, path, &pathseq, &pathseq_c);
          for (int i = 0; i < pathseq_len; ++i)
            pathseq[i] = to_int[pathseq[i]] - 1;
          int prefix_len_onlastvertex =
              graph->vertices[path->vertices[path->l - 1]]->l - offb - klen;

          // if (verbose) {
          //   fprintf(stderr, "--- %s ---\n", path->idx);
          //   fprintf(stderr, "C: %d %s\n", cons_l,
          //           decode(cons, cons_l, 1).c_str());
          //   fprintf(stderr, "P: %d %s\n",
          //           pathseq_len - offa - prefix_len_onlastvertex,
          //           decode((uint8_t *)pseq + offa,
          //                  pathseq_len - offa - prefix_len_onlastvertex, 1)
          //               .c_str());
          // }

          ksw_extd2_sse(0, cons_len, cons_seq,
                        pathseq_len - offa - prefix_len_onlastvertex,
                        (uint8_t *)pathseq + offa, 5, mat, gapo, gape, gapo2,
                        gape2, -1, -1, -1, 0, &ez);

          if (ez.score <= score)
            continue;

          score = ez.score;
          bestpath_idx = path_idx;

          if (ez.m_cigar > m_cigar) {
            cigar = (uint32_t *)realloc(cigar, ez.m_cigar * sizeof(uint32_t));
            m_cigar = ez.m_cigar;
            // TODO: check reallocation
          }
          memcpy(cigar, ez.cigar, ez.m_cigar * sizeof(uint32_t));
          n_cigar = ez.n_cigar;

          if (pathseq_len + 1 > bestpathseq_c) {
            bestpathseq = (uint8_t *)realloc(bestpathseq, pathseq_len + 1);
            bestpathseq_c = pathseq_len + 1;
          }
          bestpathseq_len = pathseq_len;
          bestpath_lastprefix = prefix_len_onlastvertex;
          memcpy(bestpathseq, (uint8_t *)pathseq, pathseq_len);
          bestpathseq[pathseq_len] = '\0';

          // pnames = collp.second;
        }

        int clipped = 0;
        if ((cigar[0] & 0xf) != 0 || (cigar[n_cigar - 1] & 0xf) != 0) {
          clipped = 1;
          continue; // XXX: do we want these?
        }

        // OUTPUT
        int opl;
        int tot_cigar_len = 0;
        int tot_res_matches = 0;

        // XXX: avoid reallocation at each iteration
        char *cs = (char *)malloc((cons_len + bestpathseq_len) * 4);
        int cs_p = 0;
        int cons_p = 0;
        int pseq_p = offa;
        int cp = 0;

        // We could do this while building difference strings, but it's better
        // for debugging purpose to keep them separated
        // for (int i = 0; i < n_cigar; ++i) {
        //   opl = cigar[i] >> 4;
        //   tot_cigar_len += opl;
        //   cp += sprintf(cigar_s + cp, "%d%c", opl, "MID"[cigar[i] & 0xf]);
        // }
        // if (verbose)
        //   fprintf(stderr, "%s\n", cigar_s);

        // XXX: avoid string copies
        for (int i = 0; i < n_cigar; ++i) {
          opl = cigar[i] >> 4;
          tot_cigar_len += opl;

          // cigar
          cp += sprintf(cigar_s + cp, "%d%c", opl, "MID"[cigar[i] & 0xf]);

          // difference string, residues, cigar length
          if ((cigar[i] & 0xf) == 0) {
            // M
            int l_tmp = 0;
            for (int j = 0; j < opl; ++j) {
              if (cons_seq[cons_p + j] != bestpathseq[pseq_p + j]) {
                if (l_tmp > 0) {
                  cs_p += sprintf(cs + cs_p, ":%d", l_tmp);
                  tot_res_matches += l_tmp;
                  l_tmp = 0;
                }
                cs_p += sprintf(cs + cs_p, "*%c%c",
                                bestpathseq[pseq_p + j] <= 3
                                    ? "ACGT"[bestpathseq[pseq_p + j]]
                                    : 'N',
                                cons_seq[cons_p + j] <= 3
                                    ? "ACGT"[cons_seq[cons_p + j]]
                                    : 'N');
              } else
                ++l_tmp;
            }
            if (l_tmp > 0) {
              tot_res_matches += l_tmp;
              cs_p += sprintf(cs + cs_p, ":%d", l_tmp);
            }
            cons_p += opl;
            pseq_p += opl;
          } else if ((cigar[i] & 0xf) == 1) {
            // I
            string tmp = decode(cons_seq + cons_p, opl, 1);
            // char *tmp = (char *)malloc(opl+1);
            // tmp[opl] = '\0';
            // strncpy(tmp, cons + cons_p, opl);
            cons_p += opl;
            // printf("+%s\n", tmp.c_str());
            cs_p += sprintf(cs + cs_p, "+%s", tmp.c_str());
            // free(tmp);
          } else if ((cigar[i] & 0xf) == 2) {
            // D
            string tmp = decode((uint8_t *)bestpathseq + pseq_p, opl, 1);
            // char *tmp = (char *)malloc(opl+1);
            // tmp[opl] = '\0';
            // strncpy(tmp, cons + cons_p, opl);
            pseq_p += opl;
            // printf("-%s\n", tmp.c_str());
            cs_p += sprintf(cs + cs_p, "-%s", tmp.c_str());
            // free(tmp);
          } else {
            fprintf(stderr,
                    "Cluster %d --- We shouldn't be here while parsing ksw "
                    "cigar (%s). Halting.\n",
                    cc_idx, cigar_s);
            exit(EXIT_FAILURE);
          }
        }
        cs[cs_p] = '\0';

        // print GAF line
        printf("%d.%d.%d\t%d\t%d\t%d\t+\t", cc_idx, sci, ci, cons_len, 0,
               cons_len);
        printf(
            ">%d",
            graph->vertices[collapsed_subpaths[bestpath_idx].first->vertices[0]]
                ->idx);
        for (int i = 1; i < path->l; ++i)
          printf(">%d", graph
                            ->vertices[collapsed_subpaths[bestpath_idx]
                                           .first->vertices[i]]
                            ->idx);
        printf("\t%d\t%d\t%d\t%d\t%d\t%d\tAS:i:%d\tcg:Z:%s\tcs:Z:%s\tcl:i:%d"
               // "\tpn:Z:%s"
               "\tcw:i:%d\tqs:Z:%s\tps:Z:%s\n",
               pathseq_len, offa, bestpathseq_len - bestpath_lastprefix,
               tot_res_matches, tot_cigar_len, 60, score, cigar_s, cs, clipped,
               // pnames.c_str(),
               abpoa_supp, decode(cons_seq, cons_len, 1).c_str(),
               decode(bestpathseq + offa,
                      bestpathseq_len - bestpath_lastprefix - offa, 1)
                   .c_str());

        free(cs);
      }
    }
    for (pair<path_t *, string> collp : collapsed_subpaths) {
      destroy_path(collp.first);
    }
  }
  fprintf(stderr, "[M::%s] analyzed %d clusters (%d left) in %.3f sec\n",
          __func__, cc_idx, sc_n - cc_idx, realtime() - rt1);

  // ---

  for (cluster_t &c : Cs) {
    for (sfs_t *s : c.specifics) {
      destroy_sfs(s);
    }
  }

  // Cleaning up
  // if (bed_f != NULL)
  //   fclose(bed_f);

  free(ez.cigar);
  free(cigar);

  abpoa_free_para(abpt);
  abpoa_free(ab);
  free(ssseqs);
  free(ssseqs_lens);
  free(cigar_s);
  free(pathseq);
  free(bestpathseq);

  destroy_graph(graph);

  fprintf(stderr, "[M::%s] done in %.3f sec\n", __func__, realtime() - rt0);

  return 0;
}
