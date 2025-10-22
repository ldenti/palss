#include <filesystem>
#include <fstream>
#include <getopt.h>
#include <string>

#include "kmc_api/kmc_file.h"

#include "graph.hpp"
#include "kmer.hpp"
#include "misc.hpp"
#include "sketch.hpp"
#include "usage.hpp"

int main_sketch(int argc, char *argv[]) {
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

  Graph graph(gbz_fn);
  graph.load();
  graph.build_fl();

  std::filesystem::create_directory(wd);

  const gbwt::Metadata &metadata = graph.get_metadata();

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
      int haplotype = metadata.fullPath(path_id).haplotype;

      haplotype = haplotype <= 0 ? 0 : haplotype - 1;
      has_haplotype[haplotype] = true;
      path_t path = graph.get_path(path_id);
      outfas[haplotype] << ">" << pp << std::endl;
      outfas[haplotype] << path.sequence << std::endl;
    }
    outfas[0].close();
    outfas[1].close();

    for (int h = 0; h < 2; ++h) {
      fprintf(stderr, "[M::%s] Counting kmers from %s.h%d\n", __func__,
              sample_name.c_str(), h + 1);
      cmd = kmc_bin + " -hp -t" + std::to_string(nThreads) + " -k" +
            std::to_string(klen) + " -ci1 -cs7 -fm " + fas_fn[h] + " " +
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
          kmc_final_db + " -ocmax";
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
          " union " + kmc_dbs[j + 3] + " -ocmax";
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
  kmer_db.SetMinCount(1);
  kmer_db.SetMaxCount(1);

  uint32_t counter;
  CKmerAPI kmer_obj(klen);
  uint64_t totkmers = 0;
  while (kmer_db.ReadNextKmer(kmer_obj, counter))
    ++totkmers;
  kmer_db.Close();

  fprintf(stderr, "[M::%s] Total kmers from %s: %ld\n", __func__,
          kmc_dbs[j + 3].c_str(), totkmers);

  int shift = std::ceil(log2(totkmers));
  sketch_t *sketch = sk_init((uint64_t)1 << shift, klen, 9); // XXX: hardcoded

  kmer_db.OpenForListing(kmc_dbs[j + 3]);
  kmer_db.SetMinCount(1);
  kmer_db.SetMaxCount(1);

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
      path_t path = graph.get_path(path_id);

      if (path.sequence.size() < klen)
        continue;

      // Get for each character, the vertex identifier
      std::vector<gbwt::size_type> info;
      for (size_t i = 0; i < path.vertices.size(); ++i) {
        size_t l = graph.get_vertex_len(path.vertices[i] >> 1);
        for (size_t x = 0; x < l; ++x)
          info.push_back(path.vertices[i] >> 1);
      }

      // we need start/end vertex for each anchor
      std::vector<uint64_t> values(path.sequence.size());
      for (size_t i = 0; i < path.sequence.size() - klen + 1; ++i) {
        uint64_t start = info[i];
        uint64_t end = info[i + klen - 1];
        values[i] = (start << 32) | end;
      }

      char *kmer = (char *)malloc(sizeof(char) *
                                  (klen + 1)); // first kmer on sequence (plain)
      uint64_t kmer_d = 0;                     // kmer
      uint64_t rckmer_d = 0;                   // reverse and complemented kmer
      uint64_t ckmer_d = 0;                    // canonical kmer
      uint8_t c;                               // new character to append
      size_t pos;                              // current position

      // first kmer
      pos = 0;
      strncpy(kmer, path.sequence.c_str(), klen);
      kmer_d = k2d(kmer, klen);
      rckmer_d = rc(kmer_d, klen);
      ckmer_d = std::min(kmer_d, rckmer_d);
      sk_add_v(sketch, ckmer_d, values[0]);

      // all other kmers
      for (pos = 1; pos < path.sequence.size() - klen + 1; ++pos) {
        c = to_int[(uint8_t)path.sequence[pos + klen - 1]] -
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
