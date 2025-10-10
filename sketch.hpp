#ifndef PS_SKETCH_H
#define PS_SKETCH_H

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>

// typedef __int128 int128_t;
typedef unsigned __int128 uint128_t;

/* sketch: 2bit encoded kmer -> (32[offset, 32[vertex in GBWT space)
 * vertex in GBWT space: [31 for ID, [1 for strand
 */

typedef struct {
  int k, m;            // kmer size, prefix size
  int64_t size;        // allocated size
  int64_t n;           // how many kmers we have
  int64_t np;          // number of prefixes
  uint64_t *pxs, *sxs; // arrays (prefixes, keys, values)
  uint32_t *vls;       // value array
} sketch_t;

/* Init the sketch */
sketch_t *sk_init(int64_t n, int k, int m);
/* Init sketch by loading from file */
sketch_t *sk_load(const std::string &fn);
/* Dump sketch to file in .txt format */
int sk_dump(sketch_t *sk, const char *fn);
/* Dump sketch to file in binary format */
int sk_store(sketch_t *sk, const char *fn);
/* Destroy the sketch */
void sk_destroy(sketch_t *sk);

/* Get position of kmer */
int64_t sk_get_p(sketch_t *sk, uint64_t kmer_d);
/* Get vertex, offset corresponding to kmer */
uint64_t sk_get(sketch_t *sk, uint64_t kmer_d);

/* Add kmer to sketch, no value */
void sk_insert(sketch_t *sk, uint64_t kmer_d, uint64_t value);

/* Add kmer to sketch, no value */
void sk_add(sketch_t *sk, uint64_t kmer_d);
/* Add value to corresponding kmer */
void sk_add_v(sketch_t *sk, uint64_t kmer_d, uint32_t value);

#endif
