#include <fstream>
#include <getopt.h>
#include <string>

#include "gbwt/fast_locate.h"
#include "gbwt/gbwt.h"
#include "gbwtgraph/gbz.h"
#include "sdsl/simple_sds.hpp"

#include "kmer.hpp"
#include "misc.hpp"
#include "sketch.hpp"
#include "usage.hpp"

typedef struct {
  uint64_t kmer;
  uint32_t vertex; // gbwt identifier, without strand information
  bool good;       // TODO: make this last bit of vertex
} solid_anchor_t;

// std::string print_kmer(uint64_t kmer_d, int klen) {
//   char *kmer = (char *)malloc(sizeof(char) * (klen + 1));
//   d23(kmer_d, klen, kmer);
//   for (int i = 0; i < klen; ++i)
//     kmer[i] = "0ACGT0"[kmer[i]];
//   return kmer;
// }

int main_sketch(int argc, char *argv[]) {
  double rt;
  uint klen = 31; // kmer size
  // int nth = 4;       // number of threads

  int _c;
  while ((_c = getopt(argc, argv, "k:@:h")) != -1) {
    switch (_c) {
    case 'k':
      klen = std::stoi(optarg);
      break;
    // case '@':
    //   nth = std::stoi(optarg);
    //   break;
    case 'h':
      fprintf(stderr, "%s", SKETCH_USAGE_MESSAGE);
      return 0;
    }
  }

  if (argc - optind != 1) {
    fprintf(stderr, "%s", SKETCH_USAGE_MESSAGE);
    return 1;
  }
  std::string gbz_fn = argv[optind++];

  rt = realtime();
  gbwtgraph::GBZ gbz;
  sdsl::simple_sds::load_from(gbz, gbz_fn);
  // gbwt::printStatistics(gbz.index, gbz_fn, std::cerr);
  fprintf(stderr, "[M::%s] Loaded GBZ in %.3fs\n", __func__, realtime() - rt);
  // R-index
  rt = realtime();
  gbwt::FastLocate fl;
  fl = gbwt::FastLocate(gbz.index);
  std::ofstream out;
  out.open(gbz_fn + ".ri", std::ofstream::out);
  fl.serialize(out);
  out.close();
  fprintf(stderr, "[M::%s] Built R-index in %.3fs\n", __func__,
          realtime() - rt);

  std::vector<solid_anchor_t> ALL_KMERS;

  const gbwt::Metadata &metadata = gbz.index.metadata;
  for (gbwt::Metadata::size_type sample_id = 0; sample_id < metadata.samples();
       ++sample_id) {
    rt = realtime();
    std::string sample_name = metadata.sample(sample_id);

    fprintf(stderr, "[M::%s] === %s ===\n", __func__, sample_name.c_str());

    std::vector<gbwt::size_type> paths = metadata.pathsForSample(sample_id);
    std::vector<int> contigs(2);
    std::vector<std::vector<solid_anchor_t>> anchors(2);
    for (size_t pp = 0; pp < paths.size(); ++pp) {
      gbwt::size_type path_id = paths[pp];
      gbwt::size_type seq_id = gbwt::Path::encode(path_id, false);
      // std::string sample_name =
      //     gbz.index.metadata.fullPath(seq_id).sample_name;
      std::string contig_name =
          gbz.index.metadata.fullPath(path_id).contig_name;
      int haplotype = gbz.index.metadata.fullPath(path_id).haplotype;
      // int offset = gbz.index.metadata.fullPath(path_id).offset;
      // std::cerr << contig_name << " " << haplotype << std::endl;

      std::string sequence = ""; // full path sequence
      std::vector<gbwt::short_type> values;
      // XXX: we do not need offsets. If we have same (curr.first>>1), then we
      // are good, each path can have the anchors at different position on same
      // vertex (due to strand)
      gbwt::edge_type curr = gbz.index.start(seq_id);
      while (curr.first != gbwt::ENDMARKER) {
        // if vertex is too long, it gets split. So we have different curr.first
        // but still same gfa_idx

        gbwtgraph::handle_t handle =
            gbwtgraph::GBWTGraph::node_to_handle(curr.first);
        std::string gfa_idx = gbz.graph.get_segment_name(handle);
        // view_type: in-place view of the sequence: (start, length)
        gbwtgraph::view_type view = gbz.graph.get_sequence_view(handle);
        sequence.append(view.first, view.second);
        // seq is already reverse&complemented if vertex was on - strand
        // int strand = (curr.first & 1); // 0 means +, 1 means -
        for (uint64_t pos = 0; pos < view.second; ++pos)
          values.push_back(curr.first >> 1);
        curr = gbz.index.LF(curr);
      }

      if (sequence.size() < klen)
        continue;

      // reference has haplotype 0, other samples have haplotypes 1 and 2, P
      // lines have haplotype -1
      haplotype = haplotype - 1 < 0 ? 0 : haplotype - 1;
      ++contigs[haplotype];

      char *kmer = (char *)malloc(sizeof(char) *
                                  (klen + 1)); // first kmer on sequence (plain)
      uint64_t kmer_d = 0;                     // kmer
      uint64_t rckmer_d = 0;                   // reverse and complemented kmer
      uint64_t ckmer_d = 0;                    // canonical kmer
      uint8_t c;                               // new character to append
      size_t pos;                              // current position

      // first kmer
      pos = 0;
      strncpy(kmer, sequence.c_str(), klen);
      kmer_d = k2d(kmer, klen);
      rckmer_d = rc(kmer_d, klen);
      ckmer_d = std::min(kmer_d, rckmer_d);
      anchors[haplotype].push_back({ckmer_d, values[0], true});

      // all other kmers
      for (pos = 1; pos < sequence.size() - klen + 1; ++pos) {
        c = to_int[(uint8_t)sequence[pos + klen - 1]] -
            1; // A is 1 but it should be 0
        kmer_d = lsappend(kmer_d, c, klen);
        rckmer_d = rsprepend(rckmer_d, reverse_char(c), klen);
        ckmer_d = std::min(kmer_d, rckmer_d);
        anchors[haplotype].push_back({ckmer_d, values[pos], true});
      }

      free(kmer);
    }

    fprintf(
        stderr,
        "[M::%s] Extracted %ld/%ld kmers from %s (%d|%d=%d contigs) in %.3fs\n",
        __func__, anchors[0].size(), anchors[1].size(), sample_name.c_str(),
        contigs[0], contigs[1], contigs[0] + contigs[1], realtime() - rt);

    rt = realtime();
    for (int hh = 0; hh < 2; ++hh) {
      // Flag all repeated anchors by sorting
      std::sort(anchors[hh].begin(), anchors[hh].end(),
                [](const solid_anchor_t &a, const solid_anchor_t &b) {
                  return a.kmer < b.kmer;
                });

      for (size_t a = 0; a < anchors[hh].size();) {
        solid_anchor_t &anchor = anchors[hh][a];
        size_t aa = a;
        ++a;
        while (a < anchors[hh].size() && anchors[hh][a].kmer == anchor.kmer) {
          anchors[hh][aa].good = false;
          if (a - 1 != aa)
            anchors[hh][a - 1].kmer = -1;
          anchors[hh][a].kmer = -1;
          ++a;
        }
        // add anchor even if it's not good (we need it for filtering out the
        // same anchor on other paths)
        ALL_KMERS.push_back(anchor);
      }
    }
    fprintf(stderr, "[M::%s] Flagged and pushed kmers in %.3fs\n", __func__,
            realtime() - rt);

    rt = realtime();
    std::sort(ALL_KMERS.begin(), ALL_KMERS.end(),
              [](const solid_anchor_t &a, const solid_anchor_t &b) {
                return a.kmer < b.kmer;
              });
    fprintf(stderr, "[M::%s] Sorted %ld kmers in %.3fs\n", __func__,
            ALL_KMERS.size(), realtime() - rt);

    rt = realtime();
    for (uint64_t a = 0; a < ALL_KMERS.size();) {
      solid_anchor_t &anchor = ALL_KMERS[a];
      uint64_t first_a = a;
      ++a;
      std::string gfa_idx = gbz.graph.get_segment_name(
          gbwtgraph::GBWTGraph::node_to_handle(anchor.vertex << 1));
      bool keep = anchor.good;
      while (a < ALL_KMERS.size() && ALL_KMERS[a].kmer == anchor.kmer) {
        keep &= gbz.graph
                        .get_segment_name(gbwtgraph::GBWTGraph::node_to_handle(
                            ALL_KMERS[a].vertex << 1))
                        .compare(gfa_idx) == 0 &&
                ALL_KMERS[a].good;
        ++a;
      }
      ALL_KMERS[first_a].good = keep;
      uint64_t j = first_a + 1;
      while (j < a) {
        ALL_KMERS[j].kmer = -1;
        ++j;
      }
    }
    // Remove gaps by shifting left
    int new_n = 0;
    uint64_t a = 0;
    while (a < ALL_KMERS.size()) {
      if (ALL_KMERS[a].kmer != (uint64_t)-1) {
        ALL_KMERS[new_n].kmer = ALL_KMERS[a].kmer;
        ALL_KMERS[new_n].vertex = ALL_KMERS[a].vertex;
        ALL_KMERS[new_n].good = ALL_KMERS[a].good;
        ++new_n;
      }
      ++a;
    }
    ALL_KMERS.resize(new_n);
    fprintf(stderr, "[M::%s] Cleaned and shifted %ld kmers in %.3fs\n",
            __func__, ALL_KMERS.size(), realtime() - rt);
  }

  fprintf(stderr, "[M::%s] === === === === === === ===\n", __func__);

  // Build and dump sketch
  rt = realtime();

  int shift = std::ceil(log2(ALL_KMERS.size()));
  std::cerr << "Shift: " << shift << std::endl;
  assert(shift < 33);
  sketch_t *sketch = sk_init((uint64_t)1 << shift, klen, 9); // XXX: hardcoded
  for (uint64_t i = 0; i < ALL_KMERS.size(); ++i) {
    if (ALL_KMERS[i].good && ALL_KMERS[i].kmer != (uint64_t)-1)
      sk_insert(sketch, ALL_KMERS[i].kmer, ALL_KMERS[i].vertex);
  }
  fprintf(stderr, "[M::%s] Added %ld kmers (%.4f) to sketch in %.3fs\n",
          __func__, sketch->n, (float)sketch->n / ALL_KMERS.size(),
          realtime() - rt);
  rt = realtime();
  sk_store(sketch, "-");
  fprintf(stderr, "[M::%s] Dumped sketch in %.3fs\n", __func__,
          realtime() - rt);
  sk_destroy(sketch);

  return 0;
}
