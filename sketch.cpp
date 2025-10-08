#include "sketch.hpp"

sketch_t *sk_init(int64_t n, int k, int m) {
  sketch_t *sk = (sketch_t *)malloc(sizeof(sketch_t));
  sk->k = k;
  sk->m = m;
  sk->n = 0;
  sk->size = n;
  sk->np = (int64_t)1 << (2 * m);
  sk->pxs = (uint64_t *)malloc(sk->np * sizeof(uint64_t));
  for (int i = 0; i < sk->np; ++i)
    sk->pxs[i] = -1;
  sk->sxs = (uint64_t *)malloc(n * sizeof(uint64_t));
  sk->vls = (uint64_t *)malloc(n * sizeof(uint64_t));
  return sk;
}

void sk_destroy(sketch_t *sk) {
  free(sk->vls);
  free(sk->sxs);
  free(sk->pxs);
  free(sk);
}

uint64_t sk_encode(uint32_t v, uint32_t offset, uint8_t strand) {
  uint64_t e = ((uint64_t)v << 32) | (offset << 1) | (strand & 1);
  return e;
}
uint32_t sk_decode_v(uint64_t e) { return e >> 32; }
uint32_t sk_decode_off(uint64_t e) { return ((uint32_t)e >> 1); }
uint8_t sk_decode_strand(uint64_t e) { return e & 1; }
uint8_t sk_decode_isref(uint64_t e) { return 0; }

/* Binary search sx in arr[l:r+1]*/
int64_t binary_search(uint64_t *arr, int64_t l, int64_t r, uint64_t sx) {
  int64_t mid;
  while (l <= r) {
    mid = l + (r - l) / 2;
    if (arr[mid] == sx)
      return mid;
    if (arr[mid] < sx)
      l = mid + 1;
    else
      r = mid - 1;
  }
  return -1;
}

int64_t sk_get_p(sketch_t *sk, uint64_t kmer_d) {
  uint64_t px = (kmer_d >> (2 * (sk->k - sk->m))) & (sk->np - 1);
  if (sk->pxs[px] == -1UL)
    return -1;

  int64_t npx = px + 1;
  // find next prefix
  while (npx < sk->np && sk->pxs[npx] == -1UL)
    ++npx;
  int64_t end;
  if (npx == sk->np)
    end = sk->n - 1;
  else
    end = sk->pxs[npx];
  return binary_search(sk->sxs, sk->pxs[px], end, kmer_d);
}

void sk_insert(sketch_t *sk, uint64_t kmer_d, uint64_t value) {
  assert(sk->n < sk->size);
  uint64_t px = (kmer_d >> (2 * (sk->k - sk->m))) & (sk->np - 1);
  if (sk->pxs[px] == -1UL) {
    // first kmer with this prefix will be at this position
    sk->pxs[px] = sk->n;
  }
  sk->sxs[sk->n] = kmer_d;
  sk->vls[sk->n] = value;
  ++sk->n;
}

// void sk_add(sketch_t *sk, uint64_t kmer_d) {
//   uint64_t px = (kmer_d >> (2 * (sk->k - sk->m))) & (sk->np - 1);
//   if (sk->pxs[px] == -1UL) {
//     // first kmer with this prefix will be at this position
//     sk->pxs[px] = sk->n;
//   }
//   sk->sxs[sk->n] = kmer_d;
//   ++sk->n;
// }

// void sk_add_v(sketch_t *sk, uint64_t kmer_d, int64_t v, int16_t offset,
//               uint8_t isref) {
//   int64_t p = sk_get_p(sk, kmer_d);
//   if (p == -1)
//     return;
//   sk->vls[p] = sk_encode(v, offset, isref);
// }

hit_t sk_get(sketch_t *sk, uint64_t kmer_d) {
  int64_t p = sk_get_p(sk, kmer_d);
  hit_t hit = {-1, -1, -1};
  if (p != -1) {
    hit.first = sk_decode_v(sk->vls[p]);
    hit.second = sk_decode_off(sk->vls[p]);
    hit.third = sk_decode_strand(sk->vls[p]);
  }
  return hit;
}

// int sk_dump(sketch_t *sk, const char *fn) {
//   FILE *fp = strcmp(fn, "-") ? fopen(fn, "w") : fdopen(fileno(stdout), "w");
//   if (fp == 0)
//     return -1;
//   fprintf(fp, "=== k: %d == m: %d == size: %ld ===\n", sk->k, sk->m, sk->n);
//   fprintf(fp, "=== PXS: ");
//   for (int i = 0; i < sk->np; ++i)
//     fprintf(fp, "%ld ", sk->pxs[i]);
//   fprintf(fp, "\n");
//   for (int i = 0; i < sk->n; ++i) {
//     fprintf(fp, "%ld,%ld,%lld\n", sk->sxs[i], sk->vls[i], sk->paths[i]);
//   }
//   fclose(fp);
//   return 0;
// }

int sk_store(sketch_t *sk, const char *fn) {
  FILE *fp = strcmp(fn, "-") ? fopen(fn, "wb") : fdopen(fileno(stdout), "wb");
  if (fp == 0)
    return -1;

  if (fwrite(&sk->k, sizeof sk->k, 1, fp) != 1) {
    fprintf(stderr, "[M::%s] failed to store sketch (k)\n", __func__);
    exit(1);
  }
  if (fwrite(&sk->m, sizeof sk->m, 1, fp) != 1) {
    fprintf(stderr, "[M::%s] failed to store sketch (m)\n", __func__);
    exit(1);
  }
  if (fwrite(&sk->n, sizeof sk->n, 1, fp) != 1) {
    fprintf(stderr, "[M::%s] failed to store sketch (n)\n", __func__);
    exit(1);
  }
  if ((int64_t)fwrite(sk->pxs, sizeof sk->pxs, sk->np, fp) != sk->np) {
    fprintf(stderr, "[M::%s] failed to store sketch (pxs)\n", __func__);
    exit(1);
  }
  if ((int64_t)fwrite(sk->sxs, sizeof sk->sxs, sk->n, fp) != sk->n) {
    fprintf(stderr, "[M::%s] failed to store sketch (sxs)\n", __func__);
    exit(1);
  }
  if ((int64_t)fwrite(sk->vls, sizeof sk->vls, sk->n, fp) != sk->n) {
    fprintf(stderr, "[M::%s] failed to store sketch (vls)\n", __func__);
    exit(1);
  }

  fclose(fp);

  return 0;
}

sketch_t *sk_load(const std::string &fn) {
  sketch_t *sk = (sketch_t *)malloc(sizeof(sketch_t));

  FILE *fp = fopen(fn.c_str(), "rb");

  if (fread(&sk->k, sizeof sk->k, 1, fp) != 1) {
    fprintf(stderr, "[M::%s] failed to read sketch (k)\n", __func__);
    exit(1);
  }

  if (fread(&sk->m, sizeof sk->m, 1, fp) != 1) {
    fprintf(stderr, "[M::%s] failed to read sketch (m)\n", __func__);
    exit(1);
  }
  if (fread(&sk->n, sizeof sk->n, 1, fp) != 1) {
    fprintf(stderr, "[M::%s] failed to read sketch (n)\n", __func__);
    exit(1);
  }
  sk->np = 1 << (2 * sk->m);
  sk->size = sk->n;

  sk->pxs = (uint64_t *)malloc(sk->np * sizeof sk->pxs);
  if ((int64_t)fread(sk->pxs, sizeof sk->pxs, sk->np, fp) != sk->np) {
    fprintf(stderr, "[M::%s] failed to read sketch (pxs)\n", __func__);
    exit(1);
  }

  sk->sxs = (uint64_t *)malloc(sk->n * sizeof sk->sxs);
  if ((int64_t)fread(sk->sxs, sizeof sk->sxs, sk->n, fp) != sk->n) {
    fprintf(stderr, "[M::%s] failed to read sketch (sxs)\n", __func__);
    exit(1);
  }

  sk->vls = (uint64_t *)malloc(sk->n * sizeof sk->vls);
  if ((int64_t)fread(sk->vls, sizeof sk->vls, sk->n, fp) != sk->n) {
    fprintf(stderr, "[M::%s] failed to read sketch (vls)\n", __func__);
    exit(1);
  }

  fclose(fp);

  return sk;
}
