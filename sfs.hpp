#ifndef PS_SFS_HPP
#define PS_SFS_HPP

#include <cstdint>
#include <string>

typedef struct {
  int s; // start on query
  int l; // length
  uint8_t flag;
  //
  uint64_t sv1, sv2;
  uint64_t ev1, ev2;
  uint64_t skmer, ekmer;
  //
  std::string rname; // plain read name
  uint8_t *seq;      // 1-4encoded sequence
  //
  std::string plain_seq;
  bool reversed;
} sfs_t;

#endif