#ifndef PSUTILS_HPP
#define PSUTILS_HPP

#include <sys/time.h>

static const uint8_t to_int[128] = {0, 0, 1, 2, 3, 0, 0, 0, 0, 0, // 0
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 10
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 20
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 30
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 40
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 50
                                    0, 0, 0, 0, 0, 1, 0, 2, 0, 0, // 60
                                    0, 3, 0, 0, 0, 0, 0, 0, 0, 0, // 70
                                    0, 0, 0, 0, 4, 0, 0, 0, 0, 0, // 80
                                    0, 0, 0, 0, 0, 0, 0, 1, 0, 2, // 90
                                    0, 0, 0, 3, 0, 0, 0, 0, 0, 0, // 100
                                    0, 0, 0, 0, 0, 0, 4, 0, 0, 0, // 110
                                    0, 0, 0, 0, 0, 0, 0, 0};      // 120

inline uint8_t reverse_char(const uint8_t c) { return ((~c) & 3); }

inline uint64_t k2d(char *kmer, uint8_t k) {
  uint64_t kmer_d = 0;
  uint8_t x;
  for (uint8_t i = 0; i < k; ++i) {
    x = (kmer[i] < 6 ? kmer[i] : to_int[kmer[i]]) -
        1; // we assume sequence to be encoded as A:1 but in kmer A is 0
    kmer_d = (kmer_d << 2) | (x < 4 ? x : rand() % 4);
  }
  return kmer_d;
}

inline uint64_t rc(uint64_t kmer, const uint8_t k) {
  uint64_t rckmer = 0;
  kmer = ~kmer;
  for (uint8_t i = 0; i < k; ++i) {
    rckmer = (rckmer << 2) | (kmer & 3);
    kmer >>= 2;
  }
  return rckmer;
}

inline uint64_t lsappend(const uint64_t kmer, const uint64_t c,
                         const uint64_t k) { // left shift and append
  return ((kmer << 2) | c) & ((1UL << 2 * k) - 1);
}

inline uint64_t rsprepend(const uint64_t kmer, const uint64_t c,
                          const uint64_t k) { // right shift and prepend
  return (kmer >> 2) | (c << (2 * k - 2));
}

static inline double realtime() {
  struct timeval tp;
  struct timezone tzp;
  gettimeofday(&tp, &tzp);
  return (double)tp.tv_sec + (double)tp.tv_usec * 1e-6;
}

#endif
