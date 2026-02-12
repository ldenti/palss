#include "reads.hpp"

__KS_BASIC(static, gzFile, 16384)
__KS_GETUNTIL(static, gzread)
__KS_INLINED(gzread)
__KSEQ_BASIC(static, gzFile)
__KSEQ_READ(static)

read_t *rd_init() {
  read_t *r = (read_t *)malloc(sizeof(read_t));
  r->name = (char *)malloc(128);
  r->name_m = 128;
  r->name_l = 0;
  r->seq = (char *)malloc(32768);
  r->seq_m = 32768;
  r->seq_l = 0;
  return r;
}

void rd_destroy(read_t *r) {
  free(r->name);
  free(r->seq);
  free(r);
}

void rd_load(read_t *r, kseq_t *seq) {
  if (seq->name.l >= r->name_m) {
    r->name = (char *)realloc(r->name, seq->name.m);
    r->name_m = seq->name.m;
  }
  memcpy(r->name, seq->name.s, seq->name.l);
  r->name[seq->name.l] = '\0';
  r->name_l = seq->name.l;

  if (seq->seq.l >= r->seq_m) {
    r->seq = (char *)realloc(r->seq, seq->seq.m);
    r->seq_m = seq->seq.m;
  }
  memcpy(r->seq, seq->seq.s, seq->seq.l);
  r->seq[seq->seq.l] = '\0';
  r->seq_l = seq->seq.l;
}

rbatch_t *rbx_init(const char *fn, int m) {
  rbatch_t *rb = (rbatch_t *)malloc(sizeof(rbatch_t));
  rb->fp = gzopen(fn, "r");
  rb->seq = kseq_init(rb->fp);
  rb->m = m;
  rb->n = 0;
  rb->reads = (read_t **)malloc(rb->m * sizeof(read_t *));
  for (size_t i = 0; i < rb->m; ++i) {
    rb->reads[i] = rd_init();
  }
  return rb;
}

void rbx_destroy(rbatch_t *rb) {
  for (size_t i = 0; i < rb->m; ++i) {
    rd_destroy(rb->reads[i]);
  }
  free(rb->reads);
  kseq_destroy(rb->seq);
  gzclose(rb->fp);
  free(rb);
}

int rbx_load(rbatch_t *rb) {
  int l = 0;
  rb->n = 0;
  while (rb->n < rb->m && (l = kseq_read(rb->seq)) >= 0) {
    rd_load(rb->reads[rb->n], rb->seq);
    ++rb->n;
  }
  return rb->n;
}
