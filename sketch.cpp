#include "sketch.hpp"

sketch_t *sk_init(int64_t n, int k, int m) {
  sketch_t *sk = (sketch_t *)malloc(sizeof(sketch_t));
  sk->k = k;
  sk->m = m;
  sk->n = 0;
  sk->size = n;
  fprintf(stderr, "[M::%s] Allocating ~%ldGB (size: %ld)\n", __func__,
          n * 128 / 8 / 1024 / 1024 / 1024, sk->size);
  sk->np = (int64_t)1 << (2 * m);
  sk->pxs = (uint64_t *)malloc(sk->np * sizeof(uint64_t));
  for (int i = 0; i < sk->np; ++i)
    sk->pxs[i] = -1;
  sk->sxs = (uint64_t *)malloc(sk->size * sizeof(uint64_t));
  sk->vls = (uint64_t *)calloc(sk->size, sizeof(uint64_t));
  sk->info = (uint32_t *)calloc(sk->size, sizeof(uint32_t));
  return sk;
}

void sk_destroy(sketch_t *sk) {
  free(sk->info);
  free(sk->vls);
  free(sk->sxs);
  free(sk->pxs);
  free(sk);
}

// Binary search sx in arr[l:r+1]
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

void sk_insert(sketch_t *sk, uint64_t kmer_d, uint32_t v1, uint32_t v2,
               uint32_t p1, uint32_t p2, uint8_t has_both,
               uint8_t is_reference) {
  assert(sk->n < sk->size);
  uint64_t px = (kmer_d >> (2 * (sk->k - sk->m))) & (sk->np - 1);
  if (sk->pxs[px] == -1UL) {
    // first kmer with this prefix will be at this position
    sk->pxs[px] = sk->n;
  }
  sk->sxs[sk->n] = kmer_d;
  //
  sk->vls[sk->n] = ((uint64_t)v1 << 32) | (v2);
  sk->info[sk->n] = ((p1 & 0x7FFF) << 17) | ((p2 & 0x7FFF) << 2) |
                    ((has_both & 1) << 1) | (is_reference & 1);
  ++sk->n;
}

hit_t sk_get(sketch_t *sk, uint64_t kmer_d, uint8_t ref_only) {
  int64_t p = sk_get_p(sk, kmer_d);
  if (p == -1)
    return {-1UL, -1U};
  // if ((sk->vls[p] & 1) == 0)
  //   return -1UL;
  return {sk->vls[p], sk->info[p]};
}

int sk_store(sketch_t *sk, const char *fn) {
  FILE *fp = strcmp(fn, "-") ? fopen(fn, "wb") : fdopen(fileno(stdout), "wb");
  if (fp == 0)
    return -1;

  if (fwrite(&sk->k, sizeof(int), 1, fp) != 1) {
    fprintf(stderr, "[M::%s] failed to store sketch (k)\n", __func__);
    exit(1);
  }
  if (fwrite(&sk->m, sizeof(int), 1, fp) != 1) {
    fprintf(stderr, "[M::%s] failed to store sketch (m)\n", __func__);
    exit(1);
  }
  if (fwrite(&sk->n, sizeof(int64_t), 1, fp) != 1) {
    fprintf(stderr, "[M::%s] failed to store sketch (n)\n", __func__);
    exit(1);
  }
  if ((int64_t)fwrite(sk->pxs, sizeof(uint64_t), sk->np, fp) != sk->np) {
    fprintf(stderr, "[M::%s] failed to store sketch (pxs)\n", __func__);
    exit(1);
  }
  if ((int64_t)fwrite(sk->sxs, sizeof(uint64_t), sk->n, fp) != sk->n) {
    fprintf(stderr, "[M::%s] failed to store sketch (sxs)\n", __func__);
    exit(1);
  }
  if ((int64_t)fwrite(sk->vls, sizeof(uint64_t), sk->n, fp) != sk->n) {
    fprintf(stderr, "[M::%s] failed to store sketch (vls)\n", __func__);
    exit(1);
  }
  if ((int64_t)fwrite(sk->info, sizeof(uint32_t), sk->n, fp) != sk->n) {
    fprintf(stderr, "[M::%s] failed to store sketch (info)\n", __func__);
    exit(1);
  }

  fclose(fp);

  return 0;
}

sketch_t *sk_load(const std::string &fn) {
  sketch_t *sk = (sketch_t *)malloc(sizeof(sketch_t));

  FILE *fp = fopen(fn.c_str(), "rb");

  if (fread(&sk->k, sizeof(int), 1, fp) != 1) {
    fprintf(stderr, "[M::%s] failed to read sketch (k)\n", __func__);
    exit(1);
  }

  if (fread(&sk->m, sizeof(int), 1, fp) != 1) {
    fprintf(stderr, "[M::%s] failed to read sketch (m)\n", __func__);
    exit(1);
  }
  if (fread(&sk->n, sizeof(int64_t), 1, fp) != 1) {
    fprintf(stderr, "[M::%s] failed to read sketch (n)\n", __func__);
    exit(1);
  }
  sk->np = 1 << (2 * sk->m);
  sk->size = sk->n;

  sk->pxs = (uint64_t *)malloc(sk->np * sizeof(uint64_t));
  if ((int64_t)fread(sk->pxs, sizeof(uint64_t), sk->np, fp) != sk->np) {
    fprintf(stderr, "[M::%s] failed to read sketch (pxs)\n", __func__);
    exit(1);
  }

  sk->sxs = (uint64_t *)malloc(sk->n * sizeof(uint64_t));
  if ((int64_t)fread(sk->sxs, sizeof(uint64_t), sk->n, fp) != sk->n) {
    fprintf(stderr, "[M::%s] failed to read sketch (sxs)\n", __func__);
    exit(1);
  }

  sk->vls = (uint64_t *)malloc(sk->n * sizeof(uint64_t));
  if ((int64_t)fread(sk->vls, sizeof(uint64_t), sk->n, fp) != sk->n) {
    fprintf(stderr, "[M::%s] failed to read sketch (vls)\n", __func__);
    exit(1);
  }

  sk->info = (uint32_t *)malloc(sk->n * sizeof(uint32_t));
  if ((int64_t)fread(sk->info, sizeof(uint32_t), sk->n, fp) != sk->n) {
    fprintf(stderr, "[M::%s] failed to read sketch (info)\n", __func__);
    exit(1);
  }

  fclose(fp);

  return sk;
}
