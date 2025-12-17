#ifndef PS_SKETCH_H
#define PS_SKETCH_H

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>

// TODO: add flags for: same vertex/different vertices. In this second case, add also which version we saw ++, +-, -+, --

inline std::pair<uint32_t, uint32_t> decode_vertices(uint64_t v) {
  return std::make_pair(v >> 32, (uint32_t)v);
}

inline std::pair<uint32_t, uint32_t> decode_positions(uint32_t i) {
  return std::make_pair((i >> 17) & 0x7FFF, (i >> 2) & 0x7FFF);
}
inline bool decode_hasboth(uint32_t i) { return (i >> 1) & 1; }
inline bool decode_type(uint32_t i) { return i & 1; }

typedef struct {
  int k, m;      // kmer size, prefix size
  int64_t size;  // allocated size
  int64_t n;     // how many kmers we have
  int64_t np;    // number of prefixes
  uint64_t *pxs; // prefixes // XXX: why pxs is 64 bits?
  // XXX: we may store only suffix and then merge with prefix
  uint64_t *sxs;  // keys (2bit-encoded kmers)
  uint64_t *vls;  // value array (gbwt-encoded vertex, w/ strand bit)
  uint32_t *info; // info array (positions 15+15, both bit, reference bit)
} sketch_t;

typedef struct {
  uint64_t value;
  uint32_t info;
} hit_t;

// Init the sketch
sketch_t *sk_init(int64_t n, int k, int m);
// Init sketch by loading from file
sketch_t *sk_load(const std::string &fn);
// Dump sketch to file in .txt format
int sk_dump(sketch_t *sk, const char *fn);
// Dump sketch to file in binary format
int sk_store(sketch_t *sk, const char *fn);
// Destroy the sketch
void sk_destroy(sketch_t *sk);

// Get position of kmer
int64_t sk_get_p(sketch_t *sk, uint64_t kmer_d);
// Get value corresponding to kmer (reference only or all)
hit_t sk_get(sketch_t *sk, uint64_t kmer_d, uint8_t ref_only);

// Add kmer to sketch
void sk_insert(sketch_t *sk, uint64_t kmer_d, uint32_t v1, uint32_t v2,
               uint32_t p1, uint32_t p2, uint8_t has_both,
               uint8_t is_reference);

#endif
