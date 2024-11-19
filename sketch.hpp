#ifndef SKETCH_HPP
#define SKETCH_HPP

#include <cstdint>
#include <map>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <utility>

#include "utils.h"

using namespace std;

// sketch : 2bit encoded kmer -> (47[vertex, 16[offset, 1[unique)
typedef map<uint64_t, uint64_t> sketch_t;
typedef pair<int64_t, int16_t> hit_t;

uint64_t sk_encode(int v, int off, int unique);
inline int8_t sk_decode_unique(uint64_t e) { return e & 1; }
inline int64_t sk_decode_v(uint64_t e) { return (e >> 17); }
inline int16_t sk_decode_off(uint64_t e) { return (e >> 1) & 0xFFFF; }

void sk_add(sketch_t &sketch, uint64_t kmer_d, uint64_t v, uint16_t offset,
            int good);
hit_t sk_get(sketch_t &sketch, uint64_t &kmer_d);
int sk_store(sketch_t &sketch, char *fn);
int sk_load(sketch_t &sketch, char *fn);

#endif
