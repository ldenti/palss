#ifndef PS_SFS_HPP
#define PS_SFS_HPP

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <kmer.hpp>

typedef struct {
  std::string rname; // plain read name
  //
  int s; // start on query
  // int rs; // start on query (reverse)
  int l; // length
  //
  uint8_t flag;
  //
  bool strand;
  //
  uint32_t sv, ev;
  uint32_t soff, eoff;
  uint64_t skmer, ekmer;
  std::vector<uint64_t> paths;
  std::vector<std::pair<uint64_t, uint64_t>> qualities;
  //
  bool swapped;
  bool reverse;
  //
  std::string plain_seq;
  uint8_t *seq; // 0-3 encoded sequence
} sfs_t;

std::string sfs_to_string(const sfs_t &s, const std::string &v1,
                          const std::string &v2);
sfs_t read_sfs_line(const std::string &line);
std::vector<sfs_t> load_sfs(const std::string &fn);

#endif
