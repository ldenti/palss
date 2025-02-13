#include "sketch.h"

sketch_t *sk_init() { return kh_init(m64); }

void sk_destroy(sketch_t *s) { kh_destroy(m64, s); }

uint64_t sk_encode(uint64_t v, uint16_t off, int unique) {
  uint64_t e = ((uint64_t)v << 17) | (off & 0xFFFF) << 1 | (unique & 1);
  return e;
}

int8_t sk_decode_unique(uint64_t e) { return e & 1; }
int64_t sk_decode_v(uint64_t e) { return (e >> 17); }
int16_t sk_decode_off(uint64_t e) { return (e >> 1) & 0xFFFF; }

void sk_add(sketch_t *sk, uint64_t kmer_d, uint64_t v, uint16_t offset,
            int good) {
  int ret;
  khiter_t k;
  k = kh_get(m64, sk, kmer_d);
  if (k == kh_end(sk)) {
    k = kh_put(m64, sk, kmer_d, &ret);
    kh_value(sk, k) = sk_encode(v, offset, good);
  } else {
    // fprintf(stderr, "XXX\n");
    kh_value(sk, k) = 0;
  }
}

hit_t sk_get(sketch_t *sk, uint64_t kmer_d) {
  khiter_t k;
  k = kh_get(m64, sk, kmer_d);
  hit_t hit = {-1, -1};
  if (k != kh_end(sk)) {
    uint64_t v = kh_value(sk, k);
    if (sk_decode_unique(v)) {
      hit.first = sk_decode_v(v);
      hit.second = sk_decode_off(v);
    }
  }
  return hit;
}

void sk_clean(sketch_t *sk) {
  khiter_t k;
  for (k = kh_begin(sk); k != kh_end(sk); ++k) {
    if (kh_exist(sk, k)) {
      /* fprintf(stderr, "%ld > %d exists\n", kh_key(sk, k), kh_value(sk, k));
       */
      if (kh_value(sk, k) == 0) {
        // fprintf(stderr, "Deleting %d\n", k);
        kh_del(m64, sk, k);
      }
    }
  }
}

int sk_store(sketch_t *sk, char *fn) {
  FILE *fp = strcmp(fn, "-") ? fopen(fn, "wb") : fdopen(fileno(stdout), "wb");
  if (fp == 0)
    return -1;

  if (fwrite(&sk->n_buckets, 4, 1, fp) != 1) {
    fprintf(stderr, "[M::%s] failed to write sketch (n_buckets)\n", __func__);
    exit(1);
  }

  if (fwrite(&sk->size, 4, 1, fp) != 1) {
    fprintf(stderr, "[M::%s] failed to write sketch (size)\n", __func__);
    exit(1);
  }

  if (fwrite(&sk->n_occupied, 4, 1, fp) != 1) {
    fprintf(stderr, "[M::%s] failed to write sketch (n_occupied)\n", __func__);
    exit(1);
  }

  if (fwrite(&sk->upper_bound, 4, 1, fp) != 1) {
    fprintf(stderr, "[M::%s] failed to write sketch (upper_bound)\n", __func__);
    exit(1);
  }

  printf("%d %d\n", __ac_fsize(sk->n_buckets), sk->n_buckets);

  int ret;
  int nf = __ac_fsize(sk->n_buckets);
  if ((ret = fwrite(sk->flags, 4, nf, fp) != nf)) {
    assert(0);
    fprintf(stderr, "[M::%s] failed to write sketch (flags)\n", __func__);
    fprintf(stderr, "Ret: %d\n", ret);
    exit(1);
  }

  if (fwrite(sk->keys, 8, sk->n_buckets, fp) != sk->n_buckets) {
    fprintf(stderr, "[M::%s] failed to write sketch (keys)\n", __func__);
    exit(1);
  }

  if (fwrite(sk->vals, 2, sk->n_buckets, fp) != sk->n_buckets) {
    fprintf(stderr, "[M::%s] failed to write sketch (vals)\n", __func__);
    exit(1);
  }

  fclose(fp);

  return 0;
}

int sk_dump(sketch_t *sk, char *fn) {
  FILE *fp = strcmp(fn, "-") ? fopen(fn, "w") : fdopen(fileno(stdout), "w");
  if (fp == 0)
    return -1;
  khiter_t k;
  for (k = kh_begin(sk); k != kh_end(sk); ++k) {
    if (kh_exist(sk, k)) {
      fprintf(fp, "%ld,%ld\n", kh_key(sk, k), kh_value(sk, k));
    }
  }
  fclose(fp);
  return 0;
}

int sk_load(sketch_t *sk, char *fn) {
  FILE *fp = fopen(fn, "rb");
  if (fread(&sk->n_buckets, sizeof(int32_t), 1, fp) != 1) {
    fprintf(stderr, "[M::%s] failed to read sketch (n_buckets)\n", __func__);
    exit(1);
  }
  if (fread(&sk->size, sizeof(int32_t), 1, fp) != 1) {
    fprintf(stderr, "[M::%s] failed to read sketch (size)\n", __func__);
    exit(1);
  }
  if (fread(&sk->n_occupied, sizeof(int32_t), 1, fp) != 1) {
    fprintf(stderr, "[M::%s] failed to read sketch (n_occupied)\n", __func__);
    exit(1);
  }
  if (fread(&sk->upper_bound, sizeof(int32_t), 1, fp) != 1) {
    fprintf(stderr, "[M::%s] failed to read sketch (upper_bound)\n", __func__);
    exit(1);
  }

  int nf = __ac_fsize(sk->n_buckets);
  sk->flags = malloc(sk->n_buckets * 4);
  if (fread(sk->flags, sizeof(int32_t), nf, fp) != nf) {
    fprintf(stderr, "[M::%s] failed to read sketch (flags)\n", __func__);
    exit(1);
  }
  sk->keys = malloc(sk->n_buckets * 8);
  if (fread(sk->keys, sizeof(int64_t), sk->n_buckets, fp) != sk->n_buckets) {
    fprintf(stderr, "[M::%s] failed to read sketch (keys)\n", __func__);
    exit(1);
  }
  sk->vals = malloc(sk->n_buckets * 8);
  if (fread(sk->vals, sizeof(int64_t), sk->n_buckets, fp) != sk->n_buckets) {
    fprintf(stderr, "[M::%s] failed to read sketch (vals)\n", __func__);
    exit(1);
  }

  fclose(fp);
  return 0;
}
