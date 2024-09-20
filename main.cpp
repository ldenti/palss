#include <algorithm>
#include <cstdint>
#include <iostream>
#include <vector>
#include <zlib.h>

#include "kmc_api/kmc_file.h"
#include "kseq.h"
#include "ropebwt3/fm-index.h"

// KSTREAM_INIT(gzFile, gzread, 65536)
KSEQ_INIT(gzFile, gzread)

using namespace std;

struct sfs_t {
  int s;
  int l;
};

// struct seg_t {
//   string idx;
//   string seq;
// };

// seg_t gfa_parse_S(char *s) {
//   seg_t ret;
//   int i, is_ok = 0;
//   char *p, *q, *seg = 0, *seq = 0, *rest = 0;
//   uint32_t sid, len = 0;
//   for (i = 0, p = q = s + 2;; ++p) {
//     if (*p == 0 || *p == '\t') {
//       int c = *p;
//       *p = 0;
//       if (i == 0) {
//         ret.idx = q;
//       } else if (i == 1) {
//         ret.seq = q;
//         is_ok = 1, rest = c ? p + 1 : 0;
//         break;
//       }
//       ++i, q = p + 1;
//       if (c == 0)
//         break;
//     }
//   }
//   // if (!is_ok) { // something is missing
//   return ret;
// }

// int extract_segments_from_gfa(int argc, char *argv[]) {
//   char *gfa_fn = argv[1];
//   int l = atoi(argv[2]);

//   kstring_t s = {0, 0, 0}, fa_seq = {0, 0, 0};
//   int dret;
//   gzFile fp = gzopen(
//       gfa_fn, "r"); // && strcmp(fn, "-") ? gzopen(fn, "r") : gzdopen(0,
//       "r");
//   if (fp == 0)
//     return 0;
//   kstream_t *ks = ks_init(fp);
//   while (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) >= 0) {
//     if (s.s[0] == 'S') {
//       seg_t seg = gfa_parse_S(s.s);
//       if (seg.seq.size() >= l)
//         cout << ">" << seg.idx << "\n" << seg.seq << endl;
//     }
//   }
//   free(fa_seq.s);
//   free(s.s);
//   ks_destroy(ks);
//   gzclose(fp);

//   return 0;
// }

vector<sfs_t> assemble_lr(const vector<sfs_t> &sfs) {
  vector<sfs_t> assembled_sfs;
  int i = 0;
  while (i < sfs.size()) {
    int j;
    for (j = i + 1; j < sfs.size(); ++j) {
      if (sfs[j - 1].s + sfs[j - 1].l <= sfs[j].s) {
        // non-overlapping
        int l = sfs[j - 1].s + sfs[j - 1].l - sfs[i].s;
        assembled_sfs.push_back({sfs[i].s, l});
        i = j;
        break;
      }
    }
    if (j == sfs.size()) {
      int l = sfs[j - 1].s + sfs[j - 1].l - sfs[i].s;
      assembled_sfs.push_back({sfs[i].s, l});
      i = j;
    }
  }
  return assembled_sfs;
}

vector<sfs_t> assemble_rl(const vector<sfs_t> &sfs) {
  vector<sfs_t> assembled_sfs;
  int i = sfs.size() - 1;
  while (i >= 0) {
    int j;
    for (j = i - 1; j >= 0; --j) {
      if (sfs[j + 1].s + sfs[j + 1].l <= sfs[j].s) {
        // non-overlapping
        int l = sfs[j + 1].s + sfs[j + 1].l - sfs[i].s;
        assembled_sfs.push_back({sfs[i].s, l});
        i = j;
        break;
      }
    }
    if (j < 0) {
      int l = sfs[j + 1].s + sfs[j + 1].l - sfs[i].s;
      assembled_sfs.push_back({sfs[i].s, l});
      i = j;
    }
  }
  return assembled_sfs;
}

/* Compute SFS strings from P and store them into solutions*/
vector<sfs_t> ping_pong_search(const rb3_fmi_t *index, uint8_t *P, int l) {
  vector<sfs_t> S;
  rb3_sai_t ik;
  int begin = l - 1;
  while (begin >= 0) {
    // Backward search. Stop at first mismatch.
    int bmatches = 0;
    rb3_fmd_set_intv(index, P[begin], &ik);
    while (ik.size != 0 && begin > 0) {
      --begin;
      ++bmatches;
      rb3_sai_t ok[RB3_ASIZE]; // output SA intervals (one for each symbol
                               // between 0 and 5)
      rb3_fmd_extend(index, &ik, ok, 1);
      ik = ok[P[begin]];
    }
    // last checked char (i.e., first of the query) was a match
    if (begin == 0 && ik.size != 0) {
      break;
    }

    // Forward search
    int end = begin;
    int fmatches = 0;
    rb3_fmd_set_intv(index, P[end], &ik);
    while (ik.size != 0) {
      ++end;
      ++fmatches;
      rb3_sai_t ok[RB3_ASIZE];
      rb3_fmd_extend(index, &ik, ok, 0);
      ik = ok[P[end] >= 1 && P[end] <= 4 ? 5 - P[end] : P[end]];
    }
    S.push_back({begin, end - begin + 1});

    if (begin == 0)
      break;

    //   if (config->overlap == 0) // Relaxed
    //     begin -= 1;
    //   else
    begin = end - 1;
  }
  return S;
}

vector<sfs_t> extend(const vector<sfs_t> &sfs, char *P, int l,
                     CKMCFile &kmc_db) {
  vector<sfs_t> extended_sfs;
  int k = kmc_db.KmerLength();
  int b, e;
  char *kmer = (char *)malloc((k + 1) * sizeof(char));
  kmer[k] = '\0';
  CKmerAPI ckmer(k);
  for (const sfs_t &s : sfs) {
    // cerr << s.s << " " << s.l << endl;
    b = s.s - k;
    b = b < 0 ? 0 : b;
    e = s.s + s.l;
    e = e > l - k + 1 ? l - k + 1 : e;
    // cerr << b << ":" << e << endl;

    strncpy(kmer, P + b, k);
    ckmer.from_string(kmer);
    while (b > 0 && !kmc_db.IsKmer(ckmer)) {
      --b;
      strncpy(kmer, P + b, k);
      ckmer.from_string(kmer);
    }

    strncpy(kmer, P + e, k);
    ckmer.from_string(kmer);
    while (e < l - k + 1 && !kmc_db.IsKmer(ckmer)) {
      ++e;
      strncpy(kmer, P + e, k);
      ckmer.from_string(kmer);
    }
    // cerr << b << ":" << e << endl;
    // cerr << endl;
    extended_sfs.push_back({b, e - b + k + 1});
  }
  free(kmer);
  return extended_sfs;
}

int main(int argc, char *argv[]) {
  char *fmd_fn = argv[1];
  char *fq_fn = argv[2];
  char *kmc_fn = argv[3];

  rb3_fmi_t f;
  rb3_fmi_restore(&f, fmd_fn, 0);
  if (f.e == 0 && f.r == 0) {
    cerr << "Error restoring index" << endl;
    return 1;
  }

  CKMCFile kmc_db;
  kmc_db.OpenForRA(kmc_fn);

  gzFile fp = gzopen(fq_fn, "r");
  kseq_t *seq = kseq_init(fp);
  int l;
  char *ss = (char *)malloc(100000); // FIXME
  uint8_t *s;
  vector<sfs_t> S;
  while ((l = kseq_read(seq)) >= 0) {
    strncpy(ss, seq->seq.s, seq->seq.l + 1);
    s = (uint8_t *)seq->seq.s;
    rb3_char2nt6(seq->seq.l, s);
    S = ping_pong_search(&f, s, seq->seq.l);
    // for (const auto &s : S)
    //   cerr << s.s << "\t" << s.l << endl;
    // cerr << "---" << endl;
    S = assemble_rl(S);
    // for (const auto &s : S)
    //   cerr << s.s << "\t" << s.l << endl;
    // S = extend(S, ss, l, kmc_db);
    // for (const auto &s : S)
    //   cerr << s.s << "\t" << s.l << endl;
    // S = assemble_lr(S);
    for (const auto &s : S)
      cout << seq->name.s << "\t" << s.s << "\t" << s.l << endl;
  }
  kseq_destroy(seq);
  gzclose(fp);
  rb3_fmi_free(&f);
  kmc_db.Close();

  cerr << "Done" << endl;
  return 0;
}
