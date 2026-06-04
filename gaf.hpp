#ifndef PS_GAF_HPP
#define PS_GAF_HPP

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

class GAFREC {
public:
  std::string qname;
  size_t qlen, qs, qe;
  bool strand; // true: +; false: -
  std::vector<std::string> path;
  size_t plen, ps, pe;
  size_t tot_res_matches, tot_cigar_len, mapq;
  int as;
  std::string cigar, cs;
  std::vector<std::string> reads;
  std::string qseq, pseq;
  bool clipped;

  GAFREC();
  void write();
};

#endif
