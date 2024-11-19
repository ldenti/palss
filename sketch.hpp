#ifndef SKETCH_HPP
#define SKETCH_HPP

#include <cstdint>
#include <map>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"

using namespace std;

// sketch : 2bit encoded kmer -> (47[vertex, 16[offset, 1[unique)
typedef map<uint64_t, uint64_t> sketch_t;

uint64_t encode(int v, int off, int unique);
inline int8_t decode_unique(uint64_t e) { return e & 1; }
inline int64_t decode_v(uint64_t e) { return (e >> 17); }
inline int16_t decode_off(uint64_t e) { return (e >> 1) & 0xFFFF; }

void add_kmer(sketch_t &sketch, uint64_t kmer_d, uint64_t v, uint16_t offset,
              int good);
int store_sketch(char *fn, sketch_t &sketch);
int load_sketch(char *fn, sketch_t &sketch);

#endif
