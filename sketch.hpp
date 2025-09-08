#ifndef PS_SKETCH_H
#define PS_SKETCH_H

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>

/* sketch: 2bit encoded kmer -> (48[vertex in GFA space, 15[offset, 1[reference)
 */

typedef struct {
  int k, m;                  // kmer size, prefix size
  int64_t size;              // allocated size
  int64_t n;                 // how many kmers we have
  int64_t np;                // number of prefixes
  uint64_t *pxs, *sxs, *vls; // arrays (prefixes, keys, values)
} sketch_t;

typedef struct {
  int64_t first;
  int64_t second;
  int64_t third;
} hit_t;

/* Encode/Decode value */
uint64_t sk_encode(int64_t v, int16_t offset, uint8_t is_ref);
int64_t sk_decode_v(uint64_t e);
int16_t sk_decode_off(uint64_t e);
uint8_t sk_decode_isref(uint64_t e);

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
hit_t sk_get(sketch_t *sk, uint64_t kmer_d);

/* Add kmer to sketch, no value */
void sk_add(sketch_t *sk, uint64_t kmer_d);
/* Add value to corresponding kmer */
void sk_add_v(sketch_t *sk, uint64_t kmer_d, int64_t v, int16_t offset,
              uint8_t isref);

#endif
