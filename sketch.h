#ifndef SKETCH_HPP
#define SKETCH_HPP

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "khash.h"

// sketch : 2bit encoded kmer -> (47[vertex, 16[offset, 1[unique)
KHASH_MAP_INIT_INT64(m64, uint64_t)
typedef khash_t(m64) sketch_t;

typedef struct {
  int64_t first;
  int64_t second;
} hit_t;

uint64_t sk_encode(uint64_t v, uint16_t offset, int unique);
int8_t sk_decode_unique(uint64_t e);
int64_t sk_decode_v(uint64_t e);
int16_t sk_decode_off(uint64_t e);

sketch_t *sk_init();
void sk_destroy(sketch_t *sk);

void sk_add(sketch_t *sk, uint64_t kmer_d, uint64_t v, uint16_t offset,
            int good);
void sk_clean(sketch_t *sk);
hit_t sk_get(sketch_t *sk, uint64_t kmer_d);
int sk_dump(sketch_t *sk, char *fn);
int sk_store(sketch_t *sk, char *fn);
int sk_load(sketch_t *sk, char *fn);

#endif
