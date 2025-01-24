#ifndef PSCL_HPP
#define PSCL_HPP

#include <stdint.h>
#include <vector>

struct cluster_t {
  std::vector<sfs_t *> specifics; // TODO: kvec
  int va = -1, vb = -1;           // starting and ending vertices
  int offa = -1, offb = -1;       // offsets on the two vertices
  uint64_t ka = -1, kb = -1;      // starting and ending kmers
};

#endif
