#ifndef PS_KMER_H
#define PS_KMER_H

#include <stdint.h>
#include <stdlib.h>

static const uint8_t to_int[128] = {1, 2, 3, 4, 0, 0, 0, 0, 0, 0, // 0
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

uint8_t reverse_char(uint8_t c);

void d23(const uint64_t kmer, int k, char *kk);

void d2s(const uint64_t kmer, int k, char *kk);
void d3s(const uint64_t kmer, int k, char *kk);

uint64_t k2d(const char *kmer, uint8_t k);

uint64_t rc(uint64_t kmer, const uint8_t k);

// left shift and append
uint64_t lsappend(const uint64_t kmer, const uint64_t c, const uint64_t k);

// right shift and prepend
uint64_t rsprepend(const uint64_t kmer, const uint64_t c, const uint64_t k);

#endif
