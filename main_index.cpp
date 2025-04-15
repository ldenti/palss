#include <cstdint>
#include <getopt.h>
#include <omp.h>
#include <string>

#include "gbwt/fast_locate.h"
#include "gbwt/gbwt.h"
#include "gbwtgraph/gbz.h"
#include "sdsl/simple_sds.hpp"

extern "C" {
#include "fm-index.h"
#include "kmer.h"
#include "ksort.h"
#include "sketch.h"
}
#include "misc.hpp"
#include "usage.hpp"

KSORT_INIT_GENERIC(uint64_t)

// Backward search until qinterval size is > than min_size
int search(const rb3_fmi_t *fmd, rb3_sai_t ik, uint8_t *kmer, int k,
           int min_size) {
  int begin = k;
  while (ik.size != 0 && begin > 0) {
    --begin;
    rb3_sai_t ok[RB3_ASIZE]; // output SA intervals (one for each symbol
                             // between 0 and 5)
    rb3_fmd_extend(fmd, &ik, ok, 1);
    ik = ok[kmer[begin]];
    if (ik.size <= min_size)
      return ik.size;
  }
  return ik.size;
}

// Backward search a given kmer and set qinterval ik
void set_qint(const rb3_fmi_t *fmd, rb3_sai_t *ik, uint8_t *kmer, int k) {
  int begin;
  begin = k - 1;
  rb3_fmd_set_intv(fmd, kmer[begin], ik);

  while (ik->size != 0 && begin > 0) {
    --begin;
    rb3_sai_t ok[RB3_ASIZE]; // output SA intervals (one for each symbol
                             // between 0 and 5)
    rb3_fmd_extend(fmd, ik, ok, 1);
    memcpy(ik, &ok[kmer[begin]], sizeof(rb3_sai_t));
  }
}

// Q-intervals memoization (all mklen-mers)
rb3_sai_t *memoize(const rb3_fmi_t *fmd, const int mklen) {
  int n = 1 << (mklen * 2);
  rb3_sai_t *qints = (rb3_sai_t *)malloc(n * sizeof(rb3_sai_t));

  char *kmer = (char *)malloc(mklen + 1);
  uint8_t *qkmer;
  kmer[mklen] = '\0';
  for (int i = 0; i < n; ++i) {
    d23(i, mklen, kmer);
    qkmer = (uint8_t *)kmer;
    rb3_char2nt6(mklen, qkmer);
    set_qint(fmd, &qints[i], qkmer, mklen);
  }
  free(kmer);
  return qints;
}

int main_index(int argc, char *argv[]) {
  double rt = realtime(), rt1;

  uint klen = 27;     // kmer size
  int mklen = 9;      // memoization kmer size
  int nh = INT32_MAX; // expected number of haplotypes
  int txt_f = 0;      // dump sketch in txt
  int big_f = 0;      // do we need big sketch?
  int nth = 4;        // number of threads

  int _c;
  while ((_c = getopt(argc, argv, "k:m:g:@:tbh")) != -1) {
    switch (_c) {
    case 'k':
      klen = std::stoi(optarg);
      break;
    case 'm':
      mklen = std::stoi(optarg);
      break;
    case 'g':
      nh = std::stoi(optarg);
      break;
    case 't':
      txt_f = 1;
      break;
    case 'b':
      big_f = 1;
      break;
    case '@':
      nth = std::stoi(optarg);
      break;
    case 'h':
      fprintf(stderr, "%s", SKETCH_USAGE_MESSAGE);
      return 0;
    default:
      fprintf(stderr, "%s", SKETCH_USAGE_MESSAGE);
      return 1;
    }
  }

  if (argc - optind != 2) {
    fprintf(stderr, "%s", SKETCH_USAGE_MESSAGE);
    return 1;
  }
  std::string gbz_fn = argv[optind++];
  std::string fmd_fn = argv[optind++];

  // R-index
  rt = realtime();
  gbwtgraph::GBZ gbz;
  sdsl::simple_sds::load_from(gbz, gbz_fn);
  fprintf(stderr, "[M::%s] restored GBZ in %.3f sec\n", __func__,
          realtime() - rt);

  rt = realtime();
  gbwt::FastLocate fl;
  fl = gbwt::FastLocate(gbz.index);
  fprintf(stderr, "[M::%s] built R-index in %.3f sec\n", __func__,
          realtime() - rt);

  rt = realtime();
  std::ofstream out;
  out.open(gbz_fn + ".ri", std::ofstream::out | std::ofstream::app);
  fl.serialize(out);
  out.close();
  fprintf(stderr, "[M::%s] stored R-index in %.3f sec\n", __func__,
          realtime() - rt);

  // Path lengths
  // XXX: is this info really not somewhere in gbwt/gbwtgraph?
  rt = realtime();
  uint npaths = gbz.index.metadata.paths();
  sdsl::int_vector<0> plens(npaths, 0, 32); // FIXME: is 32 enough?
  for (uint p = 0; p < npaths; ++p) {
    int pidx = gbwt::Path::encode(p, 0);
    plens[p] = gbz.index.extract(pidx).size();
  }
  out.open(gbz_fn + ".pl", std::ofstream::out | std::ofstream::app);
  plens.simple_sds_serialize(out);
  out.close();
  fprintf(stderr, "[M::%s] computed and stored path lengths in %.3f sec\n",
          __func__, realtime() - rt);

  // FMD-index loading
  rb3_fmi_t fmd;
  rb3_fmi_restore(&fmd, fmd_fn.c_str(), 0);
  if (fmd.e == 0 && fmd.r == 0) {
    fprintf(stderr, "[E::%s] FMD-index cannot be restored", __func__);
    exit(EXIT_FAILURE);
  }
  fprintf(stderr, "[M::%s] restored FMD index in %.3f sec\n", __func__,
          realtime() - rt);

  // Q-intervals memoization
  rt = realtime();
  int memo_bm = (1 << (mklen * 2)) - 1;
  rb3_sai_t *qints = memoize(&fmd, mklen);
  fprintf(stderr, "[M::%s] memoization (all %d-mers) in %.3f sec\n", __func__,
          mklen, realtime() - rt);
  rt = realtime();

  // First pass over graph. Store all kmers in an array
  sketch_t *sketch =
      sk_init((uint64_t)1 << (big_f ? 32 : 30), klen, mklen); // XXX: hardcoded
  uint64_t *kmers = sketch->vls; // we are going to "reuse" values array
  uint64_t totkmers = 0;         // index for insertion in kmers

  char *kmer = (char *)malloc(sizeof(char) *
                              (klen + 1)); // first kmer on sequence (plain)
  uint64_t kmer_d = 0;                     // kmer
  uint64_t rckmer_d = 0;                   // reverse and complemented kmer
  uint64_t ckmer_d = 0;                    // canonical kmer
  uint8_t c;                               // new character to append
  uint p = 0;                              // current position on segment
  int nvertices = 0;                       // number of vertices processed

  rt1 = realtime();
  // vmin is usually 1
  gbwtgraph::nid_t vmin = gbz.graph.min_node_id(),
                   vmax = gbz.graph.max_node_id();
  fprintf(stderr, "[M::%s] we have %ld vertices (%lld..%lld)\n", __func__,
          gbz.graph.get_node_count(), vmin, vmax);

  for (gbwtgraph::nid_t v = vmin; v <= vmax; ++v) {
    if (!gbz.graph.has_node(v))
      continue;
    ++nvertices;
    if (nvertices % 1000000 == 0) {
      fprintf(stderr, "[M::%s] parsed 1M vertices in %.3f sec\n", __func__,
              realtime() - rt1);
      rt1 = realtime();
    }

    gbwtgraph::handle_t vh = gbz.graph.get_handle(v);
    if (gbz.graph.get_length(vh) < klen)
      continue;
    std::string seg = gbz.graph.get_sequence(vh);
    // std::pair<std::string, std::pair<gbwtgraph::nid_t, gbwtgraph::nid_t>> x =
    //     gbz.graph.get_segment(vh);
    // std::cout << x.first << " " << seg << std::endl;

    strncpy(kmer, seg.c_str(), klen);
    kmer_d = k2d(kmer, klen);
    rckmer_d = rc(kmer_d, klen);
    ckmer_d = MIN(kmer_d, rckmer_d);
    kmers[totkmers++] = ckmer_d;

    for (p = klen; p < seg.size(); ++p) {
      c = to_int[(int)seg[p]] - 1; // A is 1 but it should be 0
      kmer_d = lsappend(kmer_d, c, klen);
      rckmer_d = rsprepend(rckmer_d, reverse_char(c), klen);
      ckmer_d = MIN(kmer_d, rckmer_d);
      kmers[totkmers++] = ckmer_d;
    }
  }
  fprintf(stderr, "[M::%s] loaded %ld kmers (from %d vertices) in %.3f sec\n",
          __func__, totkmers, nvertices, realtime() - rt);
  rt = realtime();

  // Sort kmers
  rt = realtime();
  ks_introsort(uint64_t, totkmers, kmers);
  fprintf(stderr, "[M::%s] sorted %ld kmers in %.3f sec\n", __func__, totkmers,
          realtime() - rt);

  // Flag all repeated anchors
  rt = realtime();
  for (uint64_t i = 0; i < totkmers;) {
    kmer_d = kmers[i];
    ++i;
    while (i < totkmers && kmers[i] == kmer_d) {
      kmers[i - 1] = 0;
      kmers[i] = 0;
      ++i;
    }
  }
  fprintf(stderr, "[M::%s] flagged kmers in %.3f sec\n", __func__,
          realtime() - rt);

  fprintf(stderr, "[M::%s] searching kmers in the FMD-index using %d threads\n",
          __func__, nth);
  rt = realtime();

  // Flag all non-solid anchors by querying the FMD-index
  char **kmers_t = (char **)malloc(nth * sizeof(char *));
  for (int i = 0; i < nth; ++i)
    kmers_t[i] = (char *)malloc(sizeof(char) * (klen + 1));
  uint8_t **qkmers_t = (uint8_t **)malloc(nth * sizeof(uint8_t *));
#pragma omp parallel for num_threads(nth)
  for (uint64_t i = 0; i < totkmers; ++i) {
    if (kmers[i] == 0)
      continue;
    int tidx = omp_get_thread_num();

    d23(kmers[i], klen, kmers_t[tidx]);
    qkmers_t[tidx] = (uint8_t *)kmers_t[tidx];
    rb3_char2nt6(klen, qkmers_t[tidx]);

    rb3_sai_t qint = qints[kmers[i] & memo_bm]; // memoization interval
    int hits = 1;                               // hits in the FMD index

    // XXX: assuming nh>0
    if (qint.size > nh) {
      hits = search(&fmd, qint, qkmers_t[tidx], klen - mklen, nh);
      assert(hits > 0);
    }
    if (hits > nh)
      kmers[i] = 0;
  }
  free(qkmers_t);
  for (int i = 0; i < nth; ++i)
    free(kmers_t[i]);
  free(kmers_t);
  fprintf(stderr, "[M::%s] searched kmers in %.3f sec\n", __func__,
          realtime() - rt);
  rt = realtime();

  // Add solid anchors to sketch
  for (uint64_t i = 0; i < totkmers; ++i) {
    kmer_d = kmers[i];
    if (kmer_d == 0)
      continue;
    sk_add(sketch, kmer_d);
  }
  fprintf(stderr, "[M::%s] added %ld kmers to sketch in %.3f sec\n", __func__,
          sketch->n, realtime() - rt);

  // Reiterate over graph to assign values to solid anchors
  rt = realtime();
  rt1 = rt;
  nvertices = 0;
  for (gbwtgraph::nid_t v = vmin; v <= vmax; ++v) {
    if (!gbz.graph.has_node(v))
      continue;
    ++nvertices;
    if (nvertices % 1000000 == 0) {
      fprintf(stderr, "[M::%s] parsed 1M vertices in %.3f sec\n", __func__,
              realtime() - rt1);
      rt1 = realtime();
    }

    gbwtgraph::handle_t vh = gbz.graph.get_handle(v);
    if (gbz.graph.get_length(vh) < klen)
      continue;
    std::string seg = gbz.graph.get_sequence(vh);

    strncpy(kmer, seg.c_str(), klen);
    kmer_d = k2d(kmer, klen);
    rckmer_d = rc(kmer_d, klen);
    ckmer_d = MIN(kmer_d, rckmer_d);
    sk_add_v(sketch, ckmer_d, v, 0);

    for (p = klen; p < gbz.graph.get_length(vh); ++p) {
      c = to_int[(int)seg[p]] - 1; // A is 1 but it should be 0
      kmer_d = lsappend(kmer_d, c, klen);
      rckmer_d = rsprepend(rckmer_d, reverse_char(c), klen);
      ckmer_d = MIN(kmer_d, rckmer_d);
      sk_add_v(sketch, ckmer_d, v, p - klen + 1);
    }
  }
  fprintf(stderr, "[M::%s] sketched %ld kmers in %.3f sec\n", __func__,
          sketch->n, realtime() - rt);

  // Dump sketch
  rt = realtime();
  if (txt_f)
    sk_dump(sketch, "-");
  else
    sk_store(sketch, (gbz_fn + ".skt").c_str());
  fprintf(stderr, "[M::%s] dumped sketch in %.3f sec\n", __func__,
          realtime() - rt);

  // Clean everything
  free(kmer);
  rb3_fmi_free(&fmd);
  free(qints);
  sk_destroy(sketch);

  return 0;
}
