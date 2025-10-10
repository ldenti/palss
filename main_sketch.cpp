#include <filesystem>
#include <fstream>
#include <getopt.h>
#include <string>

#include "gbwt/fast_locate.h"
#include "gbwt/gbwt.h"
#include "gbwtgraph/gbz.h"
#include "kmc_api/kmc_file.h"
#include "sdsl/simple_sds.hpp"

#include "kmer.hpp"
#include "misc.hpp"
#include "sketch.hpp"
#include "usage.hpp"

// std::string print_kmer(uint64_t kmer_d, int klen) {
//   char *kmer = (char *)malloc(sizeof(char) * (klen + 1));
//   d23(kmer_d, klen, kmer);
//   for (int i = 0; i < klen; ++i)
//     kmer[i] = "0ACGT0"[kmer[i]];
//   return kmer;
// }

int main_sketch(int argc, char *argv[]) {
  double rt;
  uint klen = 31;   // kmer size
  int nThreads = 4; // number of threads
  std::string wd = ".";

  std::string kmc_bin = "kmc";

  int _c;
  while ((_c = getopt(argc, argv, "k:@:w:h")) != -1) {
    switch (_c) {
    case 'k':
      klen = std::stoi(optarg);
      break;
    case '@':
      nThreads = std::stoi(optarg);
      break;
    case 'w':
      wd = optarg;
      break;
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
  // rt = realtime();
  // gbwt::FastLocate fl;
  // fl = gbwt::FastLocate(gbz.index);
  // std::ofstream out;
  // out.open(gbz_fn + ".ri", std::ofstream::out);
  // fl.serialize(out);
  // out.close();
  // fprintf(stderr, "[M::%s] Built R-index in %.3fs\n", __func__,
  //         realtime() - rt);

  std::filesystem::create_directory(wd);

  const gbwt::Metadata &metadata = gbz.index.metadata;
  int i = 0, j = 1;
  bool first = true;
  std::vector<std::string> fas_fn = {wd + "/paths-h1.fa", wd + "/paths-h2.fa"};
  std::vector<std::string> kmc_dbs = {wd + "/kmc-db-1", wd + "/kmc-db-2",
                                      wd + "/kmc-db-12", wd + "/kmc-db-last",
                                      wd + "/kmc-db-curr"};
  std::vector<bool> has_haplotype = {false, false};
  std::vector<std::ofstream> outfas(2);
  std::string kmc_final_db;
  std::string cmd;
  int ret;
  for (gbwt::Metadata::size_type sample_id = 0; sample_id < metadata.samples();
       ++sample_id) {

    std::string sample_name = metadata.sample(sample_id);
    fprintf(stderr, "[M::%s] Extracting %s (%ld)\n", __func__,
            sample_name.c_str(), sample_id);

    outfas[0].open(fas_fn[0], std::ios::out);
    outfas[1].open(fas_fn[1], std::ios::out);
    has_haplotype = {false, false};
    std::vector<gbwt::size_type> paths = metadata.pathsForSample(sample_id);
    for (size_t pp = 0; pp < paths.size(); ++pp) {
      gbwt::size_type path_id = paths[pp];
      gbwt::size_type seq_id = gbwt::Path::encode(path_id, false);
      int haplotype = gbz.index.metadata.fullPath(path_id).haplotype;

      haplotype = haplotype <= 0 ? 0 : haplotype - 1;
      has_haplotype[haplotype] = true;

      std::string sequence = ""; // full path sequence
      gbwt::edge_type curr = gbz.index.start(seq_id);
      while (curr.first != gbwt::ENDMARKER) {
        gbwtgraph::handle_t handle =
            gbwtgraph::GBWTGraph::node_to_handle(curr.first);
        // view_type: in-place view of the sequence: (start, length)
        gbwtgraph::view_type view = gbz.graph.get_sequence_view(handle);
        sequence.append(view.first, view.second);
        curr = gbz.index.LF(curr);
      }
      outfas[haplotype] << ">" << pp << std::endl;
      outfas[haplotype] << sequence << std::endl;
    }
    outfas[0].close();
    outfas[1].close();

    for (int h = 0; h < 2; ++h) {
      fprintf(stderr, "[M::%s] Counting kmers from %s.h%d\n", __func__,
              sample_name.c_str(), h + 1);
      cmd = kmc_bin + " -hp -t" + std::to_string(nThreads) + " -k" +
            std::to_string(klen) + " -ci1 -cx1 -fm " + fas_fn[h] + " " +
            kmc_dbs[h] + " " + wd + " &> /dev/null";
      // std::cerr << cmd << std::endl;
      ret = std::system(cmd.c_str());
      if (ret != 0) {
        std::cerr << "Error in running kmc. Return code: " << ret << std::endl;
        exit(EXIT_FAILURE);
      }
    }

    kmc_final_db = kmc_dbs[2];
    if (first) {
      kmc_final_db = kmc_dbs[i + 3];
    }

    fprintf(stderr, "[M::%s] Intersecting haplotypes from %s\n", __func__,
            sample_name.c_str());
    cmd = "kmc_tools -hp simple " + kmc_dbs[0] + " " + kmc_dbs[1] + " union " +
          kmc_final_db + " -ocmin";
    // std::cerr << cmd << std::endl;
    ret = std::system(cmd.c_str());
    if (ret != 0) {
      std::cerr << "Error in running kmc_tools. Return code: " << ret
                << std::endl;
      exit(EXIT_FAILURE);
    }

    if (first) {
      first = false;
      continue;
    }

    j = (i + 1) % 2;

    fprintf(stderr, "[M::%s] Intersecting kmc databases\n", __func__);
    cmd = "kmc_tools -hp simple " + kmc_dbs[i + 3] + " " + kmc_dbs[2] +
          " union " + kmc_dbs[j + 3] + " -ocmin";
    // std::cerr << cmd << std::endl;
    ret = std::system(cmd.c_str());
    if (ret != 0) {
      std::cerr << "Error in running kmc_tools. Return code: " << ret
                << std::endl;
      exit(EXIT_FAILURE);
    }

    i = j;
  }

  CKMCFile kmer_db;
  if (!kmer_db.OpenForListing(kmc_dbs[j + 3]))
    return -1;
  uint64_t totkmers = kmer_db.KmerCount();
  fprintf(stderr, "[M::%s] Total kmers from %s: %ld\n", __func__,
          kmc_dbs[j + 3].c_str(), totkmers);

  int shift = std::ceil(log2(totkmers));
  sketch_t *sketch = sk_init((uint64_t)1 << shift, klen, 9); // XXX: hardcoded

  uint32_t counter;
  CKmerAPI kmer_obj(klen);
  char kmer_str[klen + 1];
  while (kmer_db.ReadNextKmer(kmer_obj, counter)) {
    // XXX: can we get uint64_t from kmer_obj?
    kmer_obj.to_string(kmer_str);
    sk_add(sketch, k2d(kmer_str, klen));
  }
  kmer_db.Close();

  for (gbwt::Metadata::size_type sample_id = 0; sample_id < metadata.samples();
       ++sample_id) {
    //   rt = realtime();
    std::string sample_name = metadata.sample(sample_id);

    fprintf(stderr, "[M::%s] === %s ===\n", __func__, sample_name.c_str());

    std::vector<gbwt::size_type> paths = metadata.pathsForSample(sample_id);
    for (size_t pp = 0; pp < paths.size(); ++pp) {
      gbwt::size_type path_id = paths[pp];
      gbwt::size_type seq_id = gbwt::Path::encode(path_id, false);

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
      sk_add_v(sketch, ckmer_d, values[0]);

      // all other kmers
      for (pos = 1; pos < sequence.size() - klen + 1; ++pos) {
        c = to_int[(uint8_t)sequence[pos + klen - 1]] -
            1; // A is 1 but it should be 0
        kmer_d = lsappend(kmer_d, c, klen);
        rckmer_d = rsprepend(rckmer_d, reverse_char(c), klen);
        ckmer_d = std::min(kmer_d, rckmer_d);
        sk_add_v(sketch, ckmer_d, values[pos]);
      }
      free(kmer);
    }
  }

  sk_store(sketch, "-");
  sk_destroy(sketch);

  return 0;
}
