#include "sketch.hpp"

uint64_t sk_encode(int v, int off, int unique) {
  uint64_t e = ((uint64_t)v << 17) | (off & 0xFFFF) << 1 | (unique & 1);
  return e;
}

void sk_add(sketch_t &sketch, uint64_t kmer_d, uint64_t v, uint16_t offset,
            int good) {
  auto x = sketch.find(kmer_d);
  sketch[kmer_d] = sk_encode(v, offset, good && x == sketch.end());
}

pair<int64_t, int16_t> sk_get(sketch_t &sketch, uint64_t &kmer_d) {
  auto x = sketch.find(kmer_d);
  pair<int64_t, int16_t> hit = make_pair((int64_t)-1, (int16_t)-1);
  if (x != sketch.end() && sk_decode_unique(x->second))
    hit = make_pair(sk_decode_v(x->second), sk_decode_off(x->second));
  return hit;
}

int sk_store(sketch_t &sketch, char *fn) {
  double rt = realtime();
  FILE *fp = strcmp(fn, "-") ? fopen(fn, "wb") : fdopen(fileno(stdout), "wb");
  if (fp == 0)
    return -1;
  uint64_t total = 0;
  uint64_t skipped = 0;
  for (auto &it : sketch) {
    ++total;
    if (!sk_decode_unique(it.second)) {
      ++skipped;
      continue;
    }
    // fprintf(stderr, "%lu -> %lu (%d %d:%d)\n", (uint64_t)it.first,
    // (uint64_t)it.second, (int)decode_unique(it.second),
    // (int)decode_v(it.second), (int)decode_off(it.second));
    if (fwrite(&it.first, 8, 1, fp) != 1) {
      fprintf(stderr, "[M::%s] failed to write sketch (1). Aborting...\n",
              __func__);
      exit(1);
    }
    if (fwrite(&it.second, 8, 1, fp) != 1) {
      fprintf(stderr, "[M::%s] failed to write sketch (2). Aborting...\n",
              __func__);
      exit(1);
    }
  }
  fclose(fp);
  fprintf(
      stderr,
      "[M::%s] dumped sketch (%ld kmers out of %ld, %ld skipped) in %.3f sec\n",
      __func__, total - skipped, total, skipped, realtime() - rt);

  return 0;
}

int sk_load(sketch_t &sketch, char *fn) {
  FILE *fp = fopen(fn, "rb");
  uint64_t x, y;
  while (!feof(fp)) {
    if (fread(&x, sizeof(uint64_t), 1, fp) != 1) {
      break;
      // fprintf(stderr, "[M::%s] failed to read sketch (1). Aborting...\n",
      // __func__); exit(1);
    }
    if (fread(&y, sizeof(uint64_t), 1, fp) != 1) {
      fprintf(stderr, "[M::%s] failed to read sketch. Aborting...\n", __func__);
      exit(1);
    }
    sketch[x] = y;
  }
  fclose(fp);
  return 0;
}
