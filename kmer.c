#include "kmer.h"

inline uint8_t reverse_char(const uint8_t c) { return ((~c) & 3); }

void d23(const uint64_t kmer, int k, char *kk) {
  for (int i = 1; i <= k; ++i)
    kk[i - 1] = ((kmer >> (k - i) * 2) & 3) + 1;
  kk[k] = '\0';
}

void d2s(const uint64_t kmer, int k, char *kk) {
  for (int i = 1; i <= k; ++i)
    kk[i - 1] = "ACGT"[(kmer >> (k - i) * 2) & 3];
  kk[k] = '\0';
}

uint64_t k2d(const char *kmer, uint8_t k) {
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
