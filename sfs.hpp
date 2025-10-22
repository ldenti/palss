#ifndef PS_SFS_HPP
#define PS_SFS_HPP

#include <cstdint>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

typedef struct {
  std::string rname; // plain read name
  //
  int s; // start on query
  int l; // length
  //
  uint8_t flag;
  //
  uint64_t sv1, sv2;
  uint64_t ev1, ev2;
  uint64_t skmer, ekmer;
  //
  uint8_t *seq; // 1-4encoded sequence
  //
  std::string plain_seq;
} sfs_t;

std::map<std::string, std::vector<sfs_t>> load_sfs(const std::string &fn);

#endif