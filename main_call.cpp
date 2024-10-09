#include <iostream>
#include <string>
#include <vector>
#include <zlib.h>

#include "abpoa.h"
#include "kseq.h"

// KSTREAM_INIT(gzFile, gzread, 65536)
KSEQ_INIT(gzFile, gzread)

using namespace std;

// AaCcGgTtNn ==> 0,1,2,3,4
static unsigned char _char26_table[256] = {
    0, 1,         2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4,         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4 /*'-'*/, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0,
    4, 1,         4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4,
    4, 4,         4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4,
    4, 4,         4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4,         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4,         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4,         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4,         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4,         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4,         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};


int main_call(int argc, char *argv[]) {
  char *gfa_fn = argv[1];
  char *fq_fn = argv[2];
  int mode = atoi(argv[3]);

  // INIT ABPOA
  abpoa_t *ab = abpoa_init();
  abpoa_para_t *abpt = abpoa_init_para();
  abpt->disable_seeding = 0;
  abpt->align_mode = 0; // global
  abpt->out_msa = 0;
  abpt->out_cons = 1;
  abpt->out_gfa = 0;
  // abpt->is_diploid = 1; // TODO: maybe this works now
  abpt->progressive_poa = 1;
  abpt->amb_strand = 1;

  abpoa_post_set_para(abpt);
  // abpt->match = 2;      // match score
  // abpt->mismatch = 4;   // mismatch penalty
  // abpt->gap_mode = ABPOA_CONVEX_GAP; // gap penalty mode
  // abpt->gap_open1 = 4;  // gap open penalty #1
  // abpt->gap_ext1 = 2;   // gap extension penalty #1
  // abpt->gap_open2 = 24; // gap open penalty #2
  // abpt->gap_ext2 = 1;   // gap extension penalty #2
  // gap_penalty = min{gap_open1 + gap_len*gap_ext1, gap_open2+gap_len*gap_ext2}

  gzFile fp = gzopen(fq_fn, "r");
  kseq_t *seq = kseq_init(fp);
  int l;
  int nseqs = 0;
  while ((l = kseq_read(seq)) >= 0)
    ++nseqs;
  gzclose(fp);

  int *seq_lens = (int *)malloc(sizeof(int) * nseqs);
  uint8_t **bseqs = (uint8_t **)malloc(sizeof(uint8_t *) * nseqs);
  fp = gzopen(fq_fn, "r");
  seq = kseq_init(fp);
  int i = 0;
  while ((l = kseq_read(seq)) >= 0) {
    seq_lens[i] = l;
    bseqs[i] = (uint8_t *)malloc(sizeof(uint8_t) * (l + 1));
    for (int j = 0; j < l; ++j)
      bseqs[i][j] = _char26_table[(int)seq->seq.s[j]];
    bseqs[i][l] = '\0';
    ++i;
  }
  gzclose(fp);
  kseq_destroy(seq);

  abpoa_msa(ab, abpt, nseqs, NULL, seq_lens, bseqs, NULL, NULL);
  abpoa_cons_t *abc = ab->abc;
  string cons = "";
  int *cons_l = (int *)malloc(sizeof(int) * 1);
  uint8_t **bcons = (uint8_t **)malloc(sizeof(uint8_t *) * 1);
  if (abc->n_cons > 0) {
    l = abc->cons_len[0];
    cons_l[0] = l;
    bcons[0] = (uint8_t *)malloc(sizeof(uint8_t) * (l + 1));
    memcpy(bcons[0], abc->cons_base[0], l);
    bcons[0][l] = '\0';
    for (int j = 0; j < abc->cons_len[0]; ++j)
      cons += "ACGTN"[abc->cons_base[0][j]];
  }
  cerr << cons << endl;
  for (uint i = 0; i < nseqs; ++i)
    free(bseqs[i]);
  free(bseqs);
  free(seq_lens);
  abpoa_free(ab);
  abpoa_free_para(abpt);

  // REINIT ABPOA
  abpoa_t *ab1 = abpoa_init();
  abpoa_para_t *abpt1 = abpoa_init_para();
  abpt1->disable_seeding = 1;
  abpt1->align_mode = mode; // 0: global, 1: local, 2:extension
  abpt1->out_msa = 1;
  abpt1->out_cons = 0;
  abpt1->out_gfa = 0;
  // abpt1->is_diploid = 1; // TODO: maybe this works now
  abpt1->progressive_poa = 0;
  abpt1->amb_strand = 1;
  abpt1->incr_fn = strdup(gfa_fn);

  abpoa_post_set_para(abpt1);
  abpoa_msa(ab1, abpt1, 1, NULL, cons_l, bcons, NULL, NULL);

  int j;
  abc = ab1->abc;
  vector<char> operations;
  int x, y;
  for (j = 0; j < abc->msa_len; ++j) {
    x = abc->msa_base[0][j];
    y = abc->msa_base[1][j];
    if (x != 5 && y != 5)
      // operations.push_back(x == y ? '=' : 'X');
      operations.push_back('M');
    else
      operations.push_back(x == 5 ? 'I' : 'D');
  }
  vector<pair<int, char>> cigar;
  cigar.push_back(make_pair(1, operations[0]));
  for (j = 1; j < operations.size(); j++) {
    if (operations[j] == cigar.back().second)
      ++cigar.back().first;
    else
      cigar.push_back(make_pair(1, operations[j]));
  }
  if (cigar.front().second == 'D')
    cigar.front().second = 'S';
  else if (cigar.front().second == 'I')
    cerr << "Insertion at start" << endl;
  if (cigar.back().second == 'D')
    cigar.back().second = 'S';
  else if (cigar.back().second == 'I')
    cerr << "Insertion at end" << endl;

  for (const auto &c : cigar)
    cout << c.first << c.second << " ";
  cout << endl;

  free(bcons[0]);
  free(bcons);
  free(cons_l);
  abpoa_free(ab1);
  abpoa_free_para(abpt1);

  return 0;
}
