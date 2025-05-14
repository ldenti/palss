#include <string>
#include <vector>
#include <zlib.h>

extern "C" {
#include "abpoa.h"
#include "io.h"
#include "kseq.h"
#include "ksw2.h"
#include "sfs.h"
}
#include "graph.hpp"
// #include "usage.hpp"

KSEQ_INIT(gzFile, gzread)

int main_call(int argc, char *argv[]) {
  double rt = realtime(), rt1;

  int klen = 27; // kmer size
  int min_w = 2; // minimum support for cluster
  // float lr = 0.97; // Ratio to cluster specific strings inside same cluster

  int _c;
  while ((_c = getopt(argc, argv, "k:h")) != -1) {
    switch (_c) {
    case 'k':
      klen = std::stoi(optarg);
      break;
    case 'h':
      fprintf(stderr, "%s", "CALL_USAGE_MESSAGE");
      return 0;
    default:
      fprintf(stderr, "Error\n");
      return 1;
    }
  }
  if (argc - optind != 3) {
    fprintf(stderr, "Error!\n");
    return 1;
  }

  std::string gbz_fn = argv[optind++];
  std::string sfs_fn = argv[optind++];
  std::string fq_fn = argv[optind++];

  // Graph
  Graph graph(gbz_fn);
  graph.load();
  // graph.print_stats();
  positions_t positions = graph.get_positions();

  // cluster by starting/ending solid anchors
  int nreads = 0;
  FILE *fp;
  char *line = NULL;
  size_t len = 0;
  ssize_t read;
  std::vector<sfs_t> specifics;
  fp = fopen(sfs_fn.c_str(), "r");
  if (fp == NULL)
    exit(EXIT_FAILURE);
  while ((read = getline(&line, &len, fp)) != -1) {
    if (line[0] != 'O')
      continue;
    sfs_t ss = parse_sfs_line(line + 2);
    if (ss.qidx > nreads)
      nreads = ss.qidx;
    specifics.push_back(ss);
  }
  ++nreads; // since qidx starts from 0
  fclose(fp);

  std::sort(specifics.begin(), specifics.end(),
            [](const sfs_t &a, const sfs_t &b) { return a.a.seq < b.a.seq; });

  std::vector<std::vector<sfs_t>> clusters;
  clusters.push_back({specifics[0]});
  int last_c = 0;
  for (int ss = 1; ss < specifics.size(); ++ss) {
    sfs_t sfs = specifics[ss];
    if (sfs.a.seq == clusters[last_c][0].a.seq) {
      int c = last_c;
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

  // Split specific strings by read
  std::vector<std::vector<sfs_t *>> ss_byread(nreads);
  for (auto &cluster : clusters) {
    if (cluster.size() < min_w)
      continue;
    for (sfs_t &sfs : cluster) {
      ss_byread[sfs.qidx].push_back(&sfs);
    }
  }

  // Fill specific strings sequence (from reads)
  gzFile fq_fp = gzopen(fq_fn.c_str(), "r");
  kseq_t *seq = kseq_init(fq_fp);
  uint8_t *s;
  int l;
  int qidx = 0;
  std::vector<std::string> qnames;
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
    for (sfs_t *ss : ss_byread[qidx]) {
      ss->seq = (uint8_t *)malloc((ss->l + 1) * sizeof(uint8_t));
      memcpy(ss->seq, seq->seq.s + ss->s, ss->l);
      ss->seq[ss->l] = '\0';
      // XXX: convert to 0123 since it seems abpoa works better this way in
      // diploid mode
      for (int i = 0; i < ss->l; ++i)
        --ss->seq[i];
    }
    ++qidx;
  }
  kseq_destroy(seq);
  gzclose(fq_fp);

  // --- Analyzing clusters

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

  // XXX: Assuming clusters of size <= 64
  uint8_t **cseqs = (uint8_t **)malloc(sizeof(uint8_t *) * 64);
  int *cseqs_lens = (int *)malloc(sizeof(int) * 64);
  int goods = 0;

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

  int cidx = 0;
  for (auto &cluster : clusters) {
    if (cluster.size() < min_w)
      continue;
    if (cluster.size() > 63)
      continue;
    printf("=== %d (%ld) ===\n", cidx, cluster.size());

    // print reference locus
    sfs_t ss = cluster[0];
    int v1 = ss.a.v;
    int v2 = ss.b.v;
    std::string ref1 = positions.v2ref[v1];
    std::string ref2 = positions.v2ref[v2];
    if (ref1.compare(ref2) != 0)
      continue;
    if (positions.offsets.find(v1) == positions.offsets.end() ||
        positions.offsets.find(v2) == positions.offsets.end()) {
      // TODO: find best path by intersecting paths passing trough the two
      // vertices and report over it
      continue;
    }
    int pos1 = positions.offsets.at(v1) + ss.a.offset;
    int pos2 = positions.offsets.at(v2) + ss.b.offset;
    printf("%s:%d-%d\n", ref1.c_str(), pos1 + 1, pos2 + klen);

    // === Get paths in between v1 and v2
    graph.get_subpaths(v1, v2);

    // === Consensus via abpoa
    goods = 0;
    for (sfs_t &s : cluster) {
      cseqs_lens[goods] = s.l;
      cseqs[goods] = s.seq;
      // printf(">%s.%d.%d.%d %d:%d %d\n", s.rname, goods, s.l, s.strand, s.s,
      //        s.s + s.l, s.l);
      // for (int ii = 0; ii < s.l; ++ii)
      //   printf("%c", "ACGT"[s.seq[ii]]);
      // printf("\n");
      ++goods;
    }
    assert(goods < 64);
    abpoa_msa(ab, abpt, goods, NULL, cseqs_lens, cseqs, NULL, NULL);

    abpoa_cons_t *abc = ab->abc;
    printf("#consensus: %d\n", abc->n_cons);
    for (int ci = 0; ci < abc->n_cons; ++ci) {
      int cons_len = abc->cons_len[ci];
      uint8_t *cons_seq = abc->cons_base[ci];
    }
    // === === ===

    for (sfs_t &s : cluster)
      free(s.seq);
  }

  abpoa_free_para(abpt);
  abpoa_free(ab);
  free(cseqs);
  free(cseqs_lens);

  return 0;
}
