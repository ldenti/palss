#include <assert.h>
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

extern "C" {
#include "graph.h"
}
#include "sketch.hpp"
#include "utils.h"

#include "ksw2.h" // XXX: redefinitions

// KSEQ_INIT(gzFile, gzread) // we already init kstream in graph.h
// XXX: there should be a better way to do this
__KSEQ_TYPE(gzFile)
__KSEQ_BASIC(static, gzFile)
__KSEQ_READ(static)

using namespace std;

/* Add specific string s to cluster c and update its anchors */
void cl_insert(cluster_t &c, const sfs_t &s) {
  c.specifics.push_back(s);
  if (c.va == -1 || s.a.v < c.va) {
    c.va = s.a.v;
    c.offa = s.a.offset;
    c.ka = s.a.seq;
  }
  if (c.vb == -1 || s.b.v >= c.vb) {
    if (s.b.v == c.vb) {
      if (s.b.offset > c.offb) {
        c.offb = s.b.offset;
        c.kb = s.b.seq;
      }
    } else {
      c.vb = s.b.v;
      c.offb = s.b.offset;
      c.kb = s.b.seq;
    }
  }
}

/* Sweep line clustering of specific strings based on their anchors. Assuming
 * topological ordering of vertices */
vector<cluster_t> cluster(const vector<sfs_t> SS, int klen) {
  vector<cluster_t> clusters(1);
  cl_insert(clusters.back(), SS[0]);
  for (int i = 1; i < SS.size(); ++i) {
    if (SS[i].a.v > clusters.back().vb ||
        (SS[i].a.v == clusters.back().vb &&
         SS[i].a.offset >= clusters.back().offb + klen)) {

      // no overlap, so new cluster
      clusters.push_back({});
    }
    cl_insert(clusters.back(), SS[i]);
  }
  return clusters;
}

set<uint64_t> get_kmers(uint8_t *seq, int l, int klen) {
  // XXX: seq is on + strand. We may use sequences as they are (without r&c) and
  // use canonical kmers here
  set<uint64_t> kmers;
  char *kmer = (char *)malloc(sizeof(char) *
                              (klen + 1)); // first kmer on sequence (plain)
  uint64_t kmer_d = 0;                     // kmer
  // uint64_t rckmer_d = 0; // reverse and complemented kmer
  // uint64_t ckmer_d = 0;  // canonical kmer
  uint8_t c; // new character to append
  int p = 0; // current position on seq
  memcpy(kmer, seq, klen);
  kmer[klen] = '\0';
  kmer_d = k2d(kmer, klen);
  kmers.insert(kmer_d);
  for (p = klen; p < l; ++p) {
    c = seq[p] - 1; // A is 1 but it should be 0
    kmer_d = lsappend(kmer_d, c, klen);
    kmers.insert(kmer_d);
  }
  free(kmer);
  return kmers;
}

float jaccard(set<uint64_t> A, set<uint64_t> B) {
  set<uint64_t> intersect;
  set_intersection(A.begin(), A.end(), B.begin(), B.end(),
                   inserter(intersect, intersect.begin()));
  int den = A.size() + B.size() - intersect.size();
  if (den == 0)
    return 1.0;
  return intersect.size() / (float)den;
}

vector<set<int>> get_clusters(float **M, int size) {
  vector<set<int>> clusters;
  map<int, int> membership;
  int mb;
  for (int ri = 0; ri < size; ++ri) {
    // XXX: if ri is already in membership, skip the row
    for (int ci = ri + 1; ci < size; ++ci) {
      if (M[ri][ci] == 1.0) {
        if (clusters.size() == 0 || membership.find(ri) == membership.end()) {
          clusters.push_back(set<int>());
          mb = clusters.size() - 1;
          clusters[mb].insert(ri);
          clusters[mb].insert(ci);
          membership[ri] = mb;
          membership[ci] = mb;
        } else {
          mb = membership[ri];
          clusters[mb].insert(ci);
          membership[ci] = mb;
        }
      }
    }
  }
  return clusters;
}

int build_consensus(abpoa_t *ab, abpoa_para_t *abpt,
                    const vector<sfs_t *> &specifics, vector<int> tokeep,
                    char **cons, int *cons_c) {

  return 0; // cons_l;
}

int main_call(int argc, char *argv[]) {
  double rt0 = realtime();
  double rt = rt0, rt1;

  int klen = 27;                    // kmer size
  int min_w = 2;                    // minimum support for cluster
  int min_l = 50;                   // minimum SV length
  static ko_longopt_t longopts[] = {/*{ "search-only", ko_no_argument, 301 },*/
                                    /*{"sfs", ko_required_argument, 302},*/
                                    {NULL, 0, 0}};
  ketopt_t opt = KETOPT_INIT;
  int _c;
  while ((_c = ketopt(&opt, argc, argv, 1, "k:w:l:", longopts)) >= 0) {
    if (_c == 'k')
      klen = atoi(opt.arg);
    else if (_c == 'w')
      min_w = atoi(opt.arg);
    else if (_c == 'l')
      min_l = atoi(opt.arg);
  }
  if (argc - opt.ind != 3) {
    fprintf(stderr, "Argh");
    return 1;
  }
  char *gfa_fn = argv[opt.ind++];
  char *sfs_fn = argv[opt.ind++];
  char *fq_fn = argv[opt.ind++];

  // Graph sketching and path extraction
  sketch_t sketch;
  // sk_load(sketch, skt_fn);
  // fprintf(stderr, "[M::%s] loaded %ld sketches in %.3f sec\n", __func__,
  //         sketch.size(), realtime() - rt);
  // rt = realtime();

  graph_t *graph = init_graph(gfa_fn);
  load_vertices(graph);
  load_paths(graph);

  fprintf(stderr, "[M::%s] loaded %d paths in %.3f sec\n", __func__, graph->np,
          realtime() - rt);
  rt = realtime();
  rt1 = rt;
  // ---

  // Specific strings computation and anchoring
  /* At this point, specific strings are anchored. Anchors follow + strand on
   * graph. If read was on -, we have reversed the specific strings so that
   * everything is on + strand
   */
  // fprintf(stderr,
  //         "%d specific strings, %d anchored, %d unanchored. After merging: %d
  //         " "anchored\n", specifics_n, anchored_n, unanchored_n, SS.size());

  // rt = realtime();

  // std::sort(SS.begin(), SS.end(), [](const sfs_t &a, const sfs_t &b) {
  //   return a.a.v < b.a.v || (a.a.v == b.a.v && a.a.offset < b.a.offset);
  // });
  // fprintf(stderr, "[M::%s] sorted specific strings in %.3f sec\n", __func__,
  //         realtime() - rt);
  // rt = realtime();
  // // ---

  // // Clustering
  // vector<cluster_t> Cs = cluster(SS, klen);
  // fprintf(stderr, "[M::%s] created %ld clusters in %.3f sec\n", __func__,
  //         Cs.size(), realtime() - rt);
  // rt = realtime();
  // // ---

  // // Calling SVs
  // // char *pseq = (char *)malloc(8192 * sizeof(char));
  // // int pseq_c = 8192;
  // // int pseq_l = 0;
  // char *cons = (char *)malloc(16384 * sizeof(char));
  // int cons_c = 16384;
  // int cons_l = 0;

  // //
  // https://github.com/lh3/minimap2/blob/69e36299168d739dded1c5549f662793af10da83/options.c#L144
  // int sc_mch = 1, sc_mis = -9, gapo = 41, gape = 1;
  // int8_t a = (int8_t)sc_mch,
  //   b = sc_mis < 0 ? (int8_t)sc_mis : -(int8_t)sc_mis; // a>0 and b<0
  // int8_t mat[25] = {a, b, b, b, 0, b, a, b, b, 0, b, b, a,
  // 		    b, 0, b, b, b, a, 0, 0, 0, 0, 0, 0};

  // int vuidx = 0;
  // char *cigar = (char *)malloc(16384 * sizeof(char)); // FIXME: hardcoded
  // char *ins = (char *)malloc(100000 * sizeof(char));  // FIXME: hardcoded
  // char *del = (char *)malloc(100000 * sizeof(char));  // FIXME: hardcoded
  // int cp, cc = 0;

  // FILE *cls_f = NULL;
  // if (cls_fn != NULL)
  //   cls_f = fopen(cls_fn, "w");

  // cluster_t C;
  // int w;
  // for (int cc = 0; cc < Cs.size(); ++cc) {
  //   C = Cs[cc];
  //   w = C.specifics.size();
  //   if (cls_f != NULL)
  //     fprintf(cls_f, "%d,%d", cc, w);
  //   if (w < min_w)
  //     continue;

  //   vector<set<uint64_t>> kmers(w);
  //   for (int ss = 0; ss < w; ++ss)
  //     kmers[ss] = get_kmers(C.specifics[ss].seq, C.specifics[ss].l, klen);

  //   float **M = (float **)malloc(w * sizeof(float *));
  //   for (int i = 0; i < w; ++i)
  //     M[i] = (float *)malloc(w * sizeof(float));

  //   for (int ss1 = 0; ss1 < w; ++ss1) {
  //     int ss2;
  //     for (ss2 = 0; ss2 < ss1; ++ss2)
  //       M[ss1][ss2] = 0.0;
  //     M[ss1][ss1] = 1.0;
  //     for (ss2 = ss1 + 1; ss2 < w; ++ss2) {
  //       M[ss1][ss2] = jaccard(kmers[ss1], kmers[ss2]);
  //     }
  //   }

  //   if (cls_f != NULL) {
  //     for (int ss1 = 0; ss1 < w; ++ss1) {
  // 	fprintf(cls_f, "%f", M[ss1][0]);
  // 	for (int ss2 = 1; ss2 < w; ++ss2) {
  // 	  fprintf(cls_f, "-%f", M[ss1][ss2]);
  // 	}
  // 	if (ss1 < w-1)
  // 	  fprintf(cls_f, "|");
  //     }
  //   }

  //   vector<set<int>> subclusters = get_clusters(M, w);
  //   if (cls_f != NULL)
  //     fprintf(cls_f, ",%d\n", subclusters.size());

  //   int sci = 0;
  //   for (const auto &subcluster : subclusters) {
  //     // XXX: selecting first as representative since jaccard is 1 for each
  //     subcluster.
  //     // TODO: abpoa if we change clustering
  //     sfs_t repr = C.specifics[*(subcluster.begin())];
  //     for(int i=0; i<repr.l; ++i)
  // 	repr.seq[i] -= 1;
  //     int va = repr.a.v, vb = repr.b.v;
  //     int offa = repr.a.offset, offb = repr.b.offset;

  //     if (cls_f != NULL)
  // 	fprintf(cls_f, "%d.%d,%d:%d>%d:%d", cc, sci, va, offa, vb, offb);

  //     // XXX: calling variations wrt reference path only
  //     path_t *path = extract_subpath(graph, graph->paths[0], va, vb);
  //     if (path == NULL) {
  // 	fprintf(stderr, "Skipping %d.%d\n", cc, sci);
  // 	continue;
  //     }

  //     char *pseq = (char *)malloc(8192 * sizeof(char));
  //     int pseq_c = 8192;
  //     int pseq_l = 0;
  //     pseq_l = get_sequence(graph, path, &pseq, &pseq_c);
  //     for (int i = 0; i < pseq_l; ++i)
  // 	pseq[i] = to_int[pseq[i]] - 1;
  //     int suffix = graph->vertices[path->vertices[path->l-1]]->l - offb -
  //     klen;

  //     if (cls_f != NULL) {
  // 	fprintf(cls_f, ",%d,%s", repr.l, decode(repr.seq, repr.l, 1).c_str());
  // 	fprintf(cls_f, ",%d,%s", pseq_l - offa - suffix, decode((uint8_t*)pseq +
  // offa, pseq_l - offa - suffix, 1).c_str());
  //     }
  //     ksw_extz_t ez;
  //     memset(&ez, 0, sizeof(ksw_extz_t));
  //     ksw_extz2_sse(0, repr.l, repr.seq,
  //      		    pseq_l - offa - suffix,
  //      		    (uint8_t*)pseq + offa,
  // 		    5, mat, gapo, gape, -1, -1,
  //      		    200, 0 /*0x40*/, &ez);

  //     // OUTPUT
  //     int opl;
  //     int op;
  //     int qp = 0;
  //     int tp = offa;
  //     cp = 0;
  //     for (int i = 0; i < ez.n_cigar; ++i) {
  // 	opl = ez.cigar[i] >> 4;
  // 	cp += sprintf(cigar + cp, "%d%c", opl, "MID"[ez.cigar[i] & 0xf]);
  //     }
  //     if(cls_f != NULL)
  // 	fprintf(cls_f, ",%s\n", cigar);

  //     for (int i = 0; i < ez.n_cigar; ++i) {
  // 	l = ez.cigar[i] >> 4;
  // 	op = ez.cigar[i] & 0xf; // 0:M, 1:I, 2:D
  // 	if (op == 0) {
  // 	  // M
  // 	  qp += l;
  // 	  tp += l;
  // 	} else if (op == 1) {
  // 	  // I
  // 	  if (l >= min_l) {
  // 	    printf("%d\t%d.%d\tINS\t%d\t%s\t", vuidx, cc, sci, l, path->idx);
  // 	    printf("%d", graph->vertices[path->vertices[0]]->idx);
  // 	    for (int i = 1; i < path->l; ++i)
  // 	      printf(">%d", graph->vertices[path->vertices[i]]->idx);
  // 	    printf("\t100\t%c\t", "ACGT"[pseq[tp - 1]]);
  // 	    memcpy(ins, repr.seq + qp, l);
  // 	    ins[l] = '\0';
  // 	    printf("%s\t%d\t%s\n", decode((uint8_t*)ins, l, 1).c_str(), offa +
  // tp, cigar);
  // 	    ++vuidx;
  // 	  }
  // 	  qp += l;
  // 	} else if (op == 2) {
  // 	  // D
  // 	  if (l >= min_l) {
  // 	    printf("%d\t%d.%d\tDEL\t%d\t%s\t", vuidx, cc, sci, l, path->idx);
  // 	    printf("%d", graph->vertices[path->vertices[0]]->idx);
  // 	    for (int i = 1; i < path->l; ++i)
  // 	      printf(">%d", graph->vertices[path->vertices[i]]->idx);
  // 	    memcpy(del, pseq + tp, l);
  // 	    del[l] = '\0';
  // 	    printf("\t100\t%s\t", decode((uint8_t*)del, l, 1).c_str());
  // 	    printf("%c\t%d\t%s\n", "ACGT"[repr.seq[qp - 1]], offa + tp, cigar);
  // 	    ++vuidx;
  // 	  }
  // 	  tp += l;
  // 	}
  //     }
  //     free(ez.cigar);
  //     ++sci;
  //     destroy_path(path);

  //     //graph_t *subgraph = extract_subgraph(graph, va, vb);
  //     //destroy_graph(subgraph);

  //     // --------------- ABPOA STUFF ---------------
  //     // ---
  //     // int *seq_lens = (int *)malloc(sizeof(int) * subcluster.size());
  //     // uint8_t **bseqs =
  //     //     (uint8_t **)malloc(sizeof(uint8_t *) * subcluster.size());
  //     // int i = 0, ii = 0;
  //     // for (const int &s_idx : subcluster) {
  //     //   seq_lens[ii] = C.specifics[s_idx].l;
  //     //   bseqs[ii] = C.specifics[s_idx].seq;
  //     //   for (i = 0; i < C.specifics[s_idx].l; ++i)
  //     //     bseqs[ii][i] = C.specifics[s_idx].seq[i] - 1;
  //     //   bseqs[ii][C.specifics[s_idx].l] = '\0';
  //     //   ++ii;
  //     // }
  //     // // INIT ABPOA
  //     // abpoa_t *ab = abpoa_init();
  //     // abpoa_para_t *abpt = abpoa_init_para();
  //     // abpt->disable_seeding = 1;
  //     // abpt->align_mode = 0; // global
  //     // abpt->out_msa = 0;
  //     // abpt->out_cons = 1;
  //     // abpt->out_gfa = 0;
  //     // // abpt->is_diploid = 1; // TODO: maybe this works now
  //     // abpt->progressive_poa = 1;
  //     // abpt->amb_strand = 0;
  //     // abpt->wb = -1;
  //     // abpt->max_n_cons = 1; // to generate 1 consensus sequences
  //     //                       // abpt->match = 2;      // match score
  //     //                       // abpt->mismatch = 4;   // mismatch penalty
  //     // // abpt->gap_mode = ABPOA_CONVEX_GAP; // gap penalty mode
  //     // // abpt->gap_open1 = 4;  // gap open penalty #1
  //     // // abpt->gap_ext1 = 2;   // gap extension penalty #1
  //     // // abpt->gap_open2 = 24; // gap open penalty #2
  //     // // abpt->gap_ext2 = 1;   // gap extension penalty #2
  //     // // gap_penalty = min{gap_open1 + gap_len*gap_ext1,
  //     // // gap_open2+gap_len*gap_ext2}

  //     // abpoa_post_set_para(abpt);

  //     // for (i = 0; i < ii; ++i)
  //     //   printf("%s\n", decode(bseqs[i], seq_lens[i], 1).c_str());
  //     // abpoa_msa(ab, abpt, ii, NULL, seq_lens, bseqs, NULL, NULL);
  //     // for (i = 0; i < ii; ++i)
  //     //   printf("%s\n", decode(bseqs[i], seq_lens[i], 1).c_str());
  //     // abpoa_cons_t *abc = ab->abc;
  //     // if (abc->n_cons > 0) {
  //     //   cons_l = abc->cons_len[0];
  //     //   if (cons_l + 1 > cons_c) {
  //     //     char *temp = (char *)realloc(cons, (cons_l + 1) * 2 *
  //     //     sizeof(char)); if (temp == NULL) {
  //     //       free(cons);
  //     //       fprintf(stderr,
  //     //               "Error while reallocating memory for consensus
  //     //               string\n");
  //     //       exit(2);
  //     //     } else {
  //     //       cons = temp;
  //     //     }
  //     //     cons_c = (cons_l + 1) * 2;
  //     //   }
  //     //   for (i = 0; i < cons_l; ++i) {
  //     //     cons[i] = "ACGTN"[abc->cons_base[0][i]];
  //     //   }
  //     //   cons[i] = '\0';
  //     // }
  //     // printf("%s\n", cons);
  //     // free(bseqs);
  //     // free(seq_lens);
  //     // ---
  //     // --------------- ABPOA STUFF ---------------
  //   }

  //   free(M);
  // }
  // if (cls_f != NULL)
  //   fclose(cls_f);
  // free(cons);
  // free(cigar);
  // free(ins);
  // free(del);

  // fprintf(stderr, "[M::%s] called %d variations in %.3f sec\n", __func__,
  // vuidx, realtime() - rt0);

  // destroy_graph(graph);

  // for (sfs_t s : SS)
  //   free(s.seq);

  fprintf(stderr, "[M::%s] done in %.3f sec\n", __func__, realtime() - rt0);

  return 0;
}
