#ifndef PS_SFS_HPP
#define PS_SFS_HPP

#include <cstdint>
#include <fstream>
#include <map>
#include <set>
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
  //
  std::string plain_seq;
  uint8_t *seq; // 1-4encoded sequence
} sfs_t;

// std::map<std::string, std::vector<sfs_t>> load_sfs(const std::string &fn);
sfs_t parse_sfs_line(const std::string &line);
std::vector<sfs_t> load_sfs(const std::string &fn);

#endif
