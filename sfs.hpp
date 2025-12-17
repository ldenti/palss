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
  int l; // length
  //
  uint8_t flag;
  //
  // bool strand;
  uint32_t sv, ev;
  uint32_t soff, eoff;
  uint64_t skmer, ekmer;
  std::vector<uint64_t> paths;
  //
  std::string plain_seq;
  uint8_t *seq; // 1-4encoded sequence
} sfs_t;

// std::map<std::string, std::vector<sfs_t>> load_sfs(const std::string &fn);
std::vector<sfs_t> load_sfs(const std::string &fn);

#endif