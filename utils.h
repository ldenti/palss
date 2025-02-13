#ifndef PS_UTILS_HPP
#define PS_UTILS_HPP

#include <string>
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

void d23(uint64_t kmer, int k, char *kk) {
  for (int i = 1; i <= k; ++i)
    kk[i - 1] = ((kmer >> (k - i) * 2) & 3) + 1;
  kk[k] = '\0';
}

std::string d2s(uint64_t kmer, int k) {
  char kk[k + 1];
  for (int i = 1; i <= k; ++i)
    kk[i - 1] = "ACGT"[(kmer >> (k - i) * 2) & 3];
  kk[k] = '\0';
  return kk;
}

uint64_t k2d(char *kmer, uint8_t k) {
  uint64_t kmer_d = 0;
  uint8_t x;
  for (uint8_t i = 0; i < k; ++i) {
    x = (kmer[i] < 6 ? kmer[i] : to_int[kmer[i]]) -
        1; // we assume sequence to be encoded as A:1 but in kmer A is 0
    kmer_d = (kmer_d << 2) | (x < 4 ? x : rand() % 4);
  }
  return kmer_d;
}

uint64_t rc(uint64_t kmer, const uint8_t k) {
  uint64_t rckmer = 0;
  kmer = ~kmer;
  for (uint8_t i = 0; i < k; ++i) {
    rckmer = (rckmer << 2) | (kmer & 3);
    kmer >>= 2;
  }
  return rckmer;
}

uint64_t lsappend(const uint64_t kmer, const uint64_t c,
                  const uint64_t k) { // left shift and append
  return ((kmer << 2) | c) & ((1UL << 2 * k) - 1);
}

uint64_t rsprepend(const uint64_t kmer, const uint64_t c,
                   const uint64_t k) { // right shift and prepend
  return (kmer >> 2) | (c << (2 * k - 2));
}

#endif
