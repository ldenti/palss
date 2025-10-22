#include <algorithm>
#include <getopt.h>
#include <omp.h>
#include <vector>

extern "C" {
#include "fm-index.h"
}

#include "misc.hpp"
#include "reads.hpp"
#include "sfs.hpp"
// #include "sketch.hpp"
#include "usage.hpp"

// Compute SFS strings from P and store them into solutions
std::vector<sfs_t> ping_pong_search(const rb3_fmi_t *index, uint8_t *P, int l,
                                    char *name) {
  std::vector<sfs_t> out;
  rb3_sai_t ik;
  int begin = l - 1;
  while (begin >= 0) {
    // Backward search. Stop at first mismatch.
    // int bmatches = 0;
    rb3_fmd_set_intv(index, P[begin], &ik);
    while (ik.size != 0 && begin > 0) {
      --begin;
      // ++bmatches;
      rb3_sai_t ok[RB3_ASIZE]; // output SA intervals (one for each symbol
                               // between 0 and 5)
      rb3_fmd_extend(index, &ik, ok, 1);
      ik = ok[P[begin]];
    }
    // last checked char (i.e., first of the query) was a match
    if (begin == 0 && ik.size != 0) {
      break;
    }

    // Forward search
    int end = begin;
    // int fmatches = 0;
    rb3_fmd_set_intv(index, P[end], &ik);
    while (ik.size != 0) {
      ++end;
      // ++fmatches;
      rb3_sai_t ok[RB3_ASIZE];
      rb3_fmd_extend(index, &ik, ok, 0);
      ik = ok[P[end] >= 1 && P[end] <= 4 ? 5 - P[end] : P[end]];
    }

    out.push_back({name, begin, end - begin + 1});

    if (begin == 0)
      break;
    // if (config->overlap == 0) // Relaxed
    //   begin -= 1;
    // else
    begin = end - 1;
  }
  return out;
}

// Merge specifics strings that are too close (d-bp apart) on read
void assemble(std::vector<sfs_t> &S, int d = 0) {
  // Reverse the vector
  std::reverse(S.begin(), S.end());

  size_t i = 0;
  while (i < S.size()) {
    size_t j;
    for (j = i + 1; j <= S.size(); ++j) {
      if (j == S.size() || S[j - 1].s + S[j - 1].l <= S[j].s - d) {
        // non-overlapping: update first, clean others
        S[i].l = S[j - 1].s + S[j - 1].l - S[i].s;
        for (size_t j2 = i + 1; j2 < j; ++j2)
          S[j2].l = 0;
        break;
      }
    }
    i = j;
  }

  // Remove gaps by shifting left and resize
  int new_n = 0;
  i = 0;
  while (i < S.size()) {
    if (S[i].l > 0) {
      S[new_n].s = S[i].s;
      S[new_n].l = S[i].l;
      ++new_n;
    }
    ++i;
  }
  S.resize(new_n);
}

int main_search(int argc, char *argv[]) {
  double rt;
  int bsize = 10000; // batch size
  int nth = 4;       // number of threads

  int _c;
  while ((_c = getopt(argc, argv, "b:@:h")) != -1) {
    switch (_c) {
    case 'b':
      bsize = std::stoi(optarg);
      break;
    case '@':
      nth = std::stoi(optarg);
      break;
    case 'h':
      fprintf(stderr, "%s", SEARCH_USAGE_MESSAGE);
      return 0;
    }
  }

  if (argc - optind != 2) {
    fprintf(stderr, "%s", SEARCH_USAGE_MESSAGE);
    return 1;
  }
  std::string fmd_fn = argv[optind++];
  std::string fx_fn = argv[optind++];

  // FMD-index loading
  rt = realtime();
  rb3_fmi_t fmd;
  rb3_fmi_restore(&fmd, fmd_fn.c_str(), 0);
  if (fmd.e == 0 && fmd.r == 0) {
    fprintf(stderr, "[E::%s] FMD-index cannot be restored", __func__);
    exit(EXIT_FAILURE);
  }
  fprintf(stderr, "[M::%s] restored FMD index in %.3f sec\n", __func__,
          realtime() - rt);

  rbatch_t *rb = rbx_init(fx_fn.c_str(), bsize);
  std::vector<std::vector<sfs_t>> output(bsize);

  int total = 0;
  int nreads = 0;
  int x;
  rt = realtime();
  while ((x = rbx_load(rb)) > 0) {
    nreads += x;
#pragma omp parallel for num_threads(nth) schedule(static, 1)
    for (int qq = 0; qq < x; ++qq) {
      uint8_t *seq = (uint8_t *)rb->reads[qq]->seq;
      int l = rb->reads[qq]->seq_l;
      rb3_char2nt6(l, seq);
      output[qq] = ping_pong_search(&fmd, seq, l, rb->reads[qq]->name);
      assemble(output[qq]);
    }

    for (int qq = 0; qq < x; ++qq) {
      for (uint j = 0; j < output[qq].size(); ++j) {
        const sfs_t &s = output[qq][j];
        ++total;
        printf("%d\t%s\t%d\t%d\t%d\n", 0, s.rname.c_str(), s.s, s.l, s.s + s.l);
      }
    }
    fprintf(stderr,
            "[M::%s] computed %d specific strings from %d reads in %.3f sec\n",
            __func__, total, nreads, realtime() - rt);
  }

  rbx_destroy(rb);
  rb3_fmi_free(&fmd);

  return 0;
}