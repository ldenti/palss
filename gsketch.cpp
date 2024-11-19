#include "gsketch.hpp"

GSK::GSK(char *fn, uint8_t _k) {
  gfa_fn = fn;
  klen = _k;
}

int GSK::build_sketch() {
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
  seg_t *seg = (seg_t *)malloc(1 * sizeof(seg_t));
  seg->l = 0;
  // seg->idx = (char *)malloc(1024);
  seg->seq = (char *)malloc(4096 * sizeof(char));
  seg->c = 4096;
  nvertices = 0;
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
      add_kmer(ckmer_d, seg->idx, 0);
      for (p = klen; p < seg->l; ++p) {
        c = to_int[seg->seq[p]] - 1; // A is 1 but it should be 0
        kmer_d = lsappend(kmer_d, c, klen);
        rckmer_d = rsprepend(rckmer_d, reverse_char(c), klen);
        ckmer_d = std::min(kmer_d, rckmer_d);
        add_kmer(ckmer_d, seg->idx, p - klen + 1);
      }
    }
  }
  free(seg->seq);
  free(seg);
  free(s.s);
  ks_destroy(ks);
  gzclose(fp);

  return 0;
}

int GSK::load_sketch(char *fn) {
  FILE *f = fopen(fn, "rb");
  uint64_t x, y;
  while (!feof(f)) {
    if (fread(&x, sizeof(uint64_t), 1, f) != 1)
      exit(1);
    if (fread(&y, sizeof(uint64_t), 1, f) != 1)
      exit(1);
    sketch[x] = y;
  }
  fclose(f);
  return 0;
}

int GSK::load_vertices() {
  kstring_t s = {0, 0, 0};
  int dret;
  gzFile fp = gzopen(gfa_fn, "r");
  if (fp == 0)
    return 0;
  kstream_t *ks = ks_init(fp);

  seg_t *seg = (seg_t *)malloc(1 * sizeof(seg_t));
  seg->l = 0;
  seg->seq = (char *)malloc(4096 * sizeof(char));
  seg->c = 4096;
  nvertices = 0;
  while (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) >= 0) {
    if (s.s[0] == 'S') {
      ++nvertices;
      gfa_parse_S(s.s, seg);
      vertices[seg->idx] = seg->seq;
    }
  }
  free(seg->seq);
  free(seg);
  free(s.s);
  ks_destroy(ks);
  gzclose(fp);
  return 0;
}

void GSK::add_kmer(uint64_t kmer_d, uint64_t v, uint16_t offset) {
  auto x = sketch.find(kmer_d);
  sketch[kmer_d] = encode(v, offset, x == sketch.end());
}

int GSK::store_sketch(FILE *f, int fa) {
  for (auto &it : sketch) {
    if (!decode_unique(it.second))
      continue;
    if (!fa) {
      if (fwrite(&it.first, sizeof(uint64_t), 1, f) != 1)
        return 1;
      if (fwrite(&it.second, sizeof(uint64_t), 1, f) != 1)
        return 1;
    } else {
      fprintf(f, ">%ld.%d\n%s\n", decode_v(it.second), decode_off(it.second),
              d2s(it.first, klen).c_str());
    }
  }
  return 0;
}

pair<int64_t, int16_t> GSK::get(uint64_t &kmer_d) {
  auto x = sketch.find(kmer_d);
  pair<int64_t, int16_t> hit = make_pair((int64_t)-1, (int16_t)-1);
  if (x != sketch.end() && decode_unique(x->second))
    hit = make_pair(decode_v(x->second), decode_off(x->second));
  return hit;
}

int GSK::load_paths() {
  kstring_t s = {0, 0, 0};
  int dret;
  gzFile fp = gzopen(gfa_fn, "r");
  if (fp == 0)
    return 0;
  kstream_t *ks = ks_init(fp);

  while (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) >= 0) {
    if (s.s[0] == 'P' || s.s[0] == 'W') {
      path_t *path = init_path(2048);
      if (s.s[0] == 'P')
        gfa_parse_P(s.s, path);
      else
        gfa_parse_W(s.s, path);

      paths.push_back(path);
    }
  }

  free(s.s);
  ks_destroy(ks);
  gzclose(fp);

  return 0;
}

void GSK::destroy_paths() {
  for (path_t *p : paths) {
    destroy_path(p);
  }
}

/* Get positions of vertices in a given path with name pidx */
map<int, int> GSK::get_positions(char *pidx) const {
  map<int, int> positions;
  int cp = 0;
  int get_first = strcmp(pidx, "") == 0;
  for (const path_t *path : paths) {
    if (get_first || strcmp(path->idx, pidx) == 0) {
      for (int i = 0; i < path->l; ++i) {
        positions[path->vertices[i]] = cp;
        cp += vertices.at(path->vertices[i]).size();
      }
      break;
    }
  }
  return positions;
}

// Operations on graph
/*vector<int> GSK::adj(int u) { return graph.at(u); }*/

vector<path_t *> GSK::get_subpaths(int x, int y) {
  vector<path_t *> subpaths(paths.size());
  for (int p = 0; p < paths.size(); ++p) {
    int f = 0, ok = 0;
    subpaths[p] = init_path(y - x + 1);
    strcpy(subpaths[p]->idx, paths[p]->idx);
    for (int i = 0; i < paths[p]->l; ++i) {
      if (paths[p]->vertices[i] == x)
        f = 1;
      if (paths[p]->vertices[i] == y)
        ok = 1;
      if (f)
        add_vertex(subpaths[p], paths[p]->vertices[i]);
      if (paths[p]->vertices[i] > y || ok)
        break;
    }
    if (!ok) {
      destroy_path(subpaths[p]);
      subpaths[p] = NULL;
    }
  }
  return subpaths;
}

int GSK::get_vl(int v) { return vertices.at(v).size(); }

int GSK::get_sequence(const path_t *path, char **pseq, int *pseq_c) {
  int l = 0;
  for (int i = 0; i < path->l; ++i)
    l += vertices[path->vertices[i]].size();
  if (l + 1 > *pseq_c) {
    cerr << "--- Reallocating from " << *pseq_c << " to " << (l + 1) * 2
         << endl;
    char *temp = (char *)realloc(*pseq, (l + 1) * 2 * sizeof(char));
    if (temp == NULL) {
      free(pseq);
      cerr << "Error while reallocating memory for path string" << endl;
      exit(2);
    } else {
      *pseq = temp;
    }
    *pseq_c = (l + 1) * 2;
  }
  int p = 0;
  for (int i = 0; i < path->l; ++i) {
    l = vertices[path->vertices[i]].size();
    strncpy(*pseq + p, vertices[path->vertices[i]].c_str(), l);
    p += l;
  }
  (*pseq)[p] = '\0';

  for (int i = 0; i < p; ++i)
    (*pseq)[i] = to_int[(*pseq)[i]] - 1;
  return p;
}

int GSK::compatible(int x, int y) {
  if (x > y) {
    int tmp = x;
    x = y;
    y = tmp;
  }
  for (int p = 0; p < paths.size(); ++p) {
    int f = 0, ok = 0;
    for (int i = 0; i < paths[p]->l; ++i) {
      if (paths[p]->vertices[i] == x)
        f = 1;
      if (f && paths[p]->vertices[i] == y)
        return f;
    }
  }
  return 0;
}

// Functions to parse GFA lines
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
        // TODO: reallocate if vertex is longer than 4096
        // right now we assume to have a vg chopped graph
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

void GSK::gfa_parse_P(char *s, path_t *path) {
  int x = 0;        // current index for insertion
  int i;            // , oriv = -1, oriw = -1, is_ok = 0;
  char *p, *q, *qq; //, *segv = 0, *segw = 0, *rest = 0;
  // int32_t ov = INT32_MAX, ow = INT32_MAX;
  for (i = 0, p = q = s + 2;; ++p) {
    if (*p == 0 || *p == '\t') {
      *p = 0;
      if (i == 0) {
        // TODO: check for duplicates
        strcpy(path->idx, q);
      } else if (i == 1) {
        char strand = *(p - 1);
        qq = q;
        for (qq = q;; ++qq) {
          if (*qq == 0 || *qq == ',') {
            int c = *qq;
            *qq = 0;
            if (*(qq - 1) != strand) {
              cerr << "Mixed +/- strands in path " << path->idx << endl;
              exit(1);
            }
            *(qq - 1) = 0;
            /*cerr << "Adding " << q << " to " << path->idx << endl;*/
            add_vertex(path, stoi(q));
            q = qq + 1;
            if (c == 0)
              break;
          }
        }
        break;
      }
      ++i, q = p + 1;
    }
  }
}

void GSK::gfa_parse_W(char *s, path_t *path) {
  int x = 0; // current index for insertion
  int i;
  char *p, *q, *qq;
  for (i = 0, p = q = s + 2;; ++p) {
    if (*p == 0 || *p == '\t') {
      *p = 0;
      if (i == 0) {
        strcpy(path->idx, q);
        *(p + 2) = 0;
        path->idx[p - q] = '#';
        path->idx[p - q + 1] = *(p + 1);
        path->idx[p - q + 2] = '\0';
      } else if (i == 5) {
        ++q;
        for (qq = q;; ++qq) {
          if (*qq == 0 || *qq == '>') {
            int c = *qq;
            *qq = 0;
            add_vertex(path, stoi(q));
            q = qq + 1;
            if (c == 0)
              break;
          }
        }
        break;
      }
      ++i, q = p + 1;
    }
  }
}
