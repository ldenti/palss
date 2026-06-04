#include "gaf.hpp"

GAFREC::GAFREC() {}

void GAFREC::write() {
  std::ostringstream oss;
  oss << qname << "\t" << qlen << "\t" << qs << "\t" << qe << "\t"
      << (strand ? "+" : "-") << "\t";
  for (const std::string &v : path)
    oss << v;

  oss << "\t" << plen << "\t" << ps << "\t" << pe << "\t" << tot_res_matches
      << "\t" << tot_cigar_len << "\t" << mapq << "\t"
      << "AS:i:" << as << "\t"
      << "cg:Z:" << cigar << "\t"
      << "cs:Z:" << cs
      << "\t"
      // << "cl:Z:" << (clipped ? 1 : 0) << "\t"
      << "rs:Z:";
  oss << reads[0];
  for (size_t x = 1; x < reads.size(); ++x)
    oss << "|" << reads[x];
  oss << "\t"
      << "qs:Z:" << qseq << "\t"
      << "ps:Z:" << pseq << "\n";

  std::cout << oss.str();
}
