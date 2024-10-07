#include "gsketch.hpp"

GSK::GSK(char *fn) { gfa_fn = fn; }

void GSK::add_kmer(uint64_t &kmer_d, int &idx) {
  auto x = sketch.find(kmer_d);
  auto y = multi.find(kmer_d);

  // assert((x == sketch.end()) == (y == multi.end()) && y != multi.end());
  // we cannot have false false here
  if (x == sketch.end() && y == multi.end()) {
    sketch[kmer_d] = idx;
  } else {
    if (x != sketch.end()) {
      sketch.erase(x);
      multi.insert(x->first);
    }
  }
}

int GSK::get(uint64_t &kmer_d) {
  auto x = sketch.find(kmer_d);
  return x != sketch.end() ? x->second : -1;
}

int GSK::build_graph() {
  kstring_t s = {0, 0, 0};
  int dret;
  gzFile fp = gzopen(gfa_fn, "r");
  if (fp == 0)
    return 0;
  kstream_t *ks = ks_init(fp);

  link_t *link = (link_t *)malloc(1 * sizeof(link_t));

  // CXXGraph::T_EdgeSet<int> edgeSet;
  cerr << nvertices << endl;
  graph = boost::adjacency_list<>(nvertices);
  int edge_idx = 0;
  while (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) >= 0) {
    if (s.s[0] == 'L') {
      gfa_parse_L(s.s, link);
      add_edge(link->idx1, link->idx2, graph);
      // CXXGraph::Node<int> node1(to_string(link->idx1), link->idx1);
      // CXXGraph::Node<int> node2(to_string(link->idx2), link->idx2);
      // CXXGraph::DirectedWeightedEdge<int> edge(edge_idx, node1, node2, 1);
      // ++edge_idx;
      // edgeSet.insert(make_shared<CXXGraph::DirectedWeightedEdge<int>>(edge));
    }
  }
  free(link);
  free(s.s);
  ks_destroy(ks);
  gzclose(fp);

  // graph = CXXGraph::Graph<int>(edgeSet);
  // graph.writeToFile();
  return 0;
}

int GSK::build_sketch(int klen) {
  k = klen;

  kstring_t s = {0, 0, 0};
  int dret;
  gzFile fp = gzopen(gfa_fn, "r");
  if (fp == 0)
    return 0;
  kstream_t *ks = ks_init(fp);

  char kmer[klen + 1];   // first kmer on sequence (plain)
  uint64_t kmer_d = 0;   // kmer
  uint64_t rckmer_d = 0; // reverse and complemented kmer
  uint64_t ckmer_d = 0;  // canonical kmer
  uint8_t c;             // new character to append
  int p = 0;             // current position on segment
  seg_t *seg = (seg_t *)malloc(1);
  seg->l = 0;
  // seg->idx = (char *)malloc(1024);
  seg->seq = (char *)malloc(4096);

  while (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) >= 0) {
    if (s.s[0] == 'S') {
      ++nvertices;
      gfa_parse_S(s.s, seg);
      if (seg->l < klen)
        continue;

      strncpy(kmer, seg->seq, klen);
      kmer_d = k2d(kmer, klen);
      rckmer_d = rc(kmer_d, klen);
      ckmer_d = std::min(kmer_d, rckmer_d);
      add_kmer(ckmer_d, seg->idx);
      for (p = klen; p < seg->l; ++p) {
        c = to_int[seg->seq[p]] - 1; // A is 1 but it should be 0
        kmer_d = lsappend(kmer_d, c, klen);
        rckmer_d = rsprepend(rckmer_d, reverse_char(c), klen);
        ckmer_d = std::min(kmer_d, rckmer_d);
        add_kmer(ckmer_d, seg->idx);
      }
      // cout << ">" << seg->idx << " " << seg->l << "\n" << seg->seq << endl;
    }
  }
  free(s.s);
  ks_destroy(ks);
  gzclose(fp);

  // for (const auto &x : sketch) {
  //   cout << x.first << " " << x.second << endl;
  // }
  multi.clear();

  return 0;
}

void GSK::gfa_parse_S(char *s, seg_t *ret) {
  int i, is_ok = 0;
  char *p, *q, *seg = 0, *seq = 0, *rest = 0;
  uint32_t sid, len = 0;
  for (i = 0, p = q = s + 2;; ++p) {
    if (*p == 0 || *p == '\t') {
      int c = *p;
      *p = 0;
      if (i == 0) {
        ret->idx = stoi(q);
        // strcpy(ret->idx, q);
      } else if (i == 1) {
        strcpy(ret->seq, q);
        ret->l = p - q;
        is_ok = 1, rest = c ? p + 1 : 0;
        break;
      }
      ++i, q = p + 1;
      if (c == 0)
        break;
    }
  }
  // if (!is_ok) { // something is missing
}

void GSK::gfa_parse_L(char *s, link_t *ret) {
  int i, oriv = -1, oriw = -1, is_ok = 0;
  char *p, *q, *segv = 0, *segw = 0, *rest = 0;
  int32_t ov = INT32_MAX, ow = INT32_MAX;
  for (i = 0, p = q = s + 2;; ++p) {
    if (*p == 0 || *p == '\t') {
      int c = *p;
      *p = 0;
      if (i == 0) {
        ret->idx1 = stoi(q);
      } else if (i == 1) {
        if (*q != '+' && *q != '-')
          ; // return -2;
        oriv = (*q != '+');
      } else if (i == 2) {
        ret->idx2 = stoi(q);
      } else if (i == 3) {
        if (*q != '+' && *q != '-')
          ; // return -2;
        oriw = (*q != '+');
      }
      ++i, q = p + 1;
      if (c == 0)
        break;
    }
  }
  // return 0;
}