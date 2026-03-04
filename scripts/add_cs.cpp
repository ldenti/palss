#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <zlib.h>

#include "gbwtgraph/gbz.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

// XXX: improve
std::string reverseAndComplement(const std::string &input) {
  std::string reversed(input.rbegin(), input.rend());
  std::string complement;
  for (char nucleotide : reversed) {
    switch (nucleotide) {
    case 'A':
      complement += 'T';
      break;
    case 'T':
      complement += 'A';
      break;
    case 'C':
      complement += 'G';
      break;
    case 'G':
      complement += 'C';
      break;
    default:
      complement += nucleotide;
      break;
    }
  }
  return complement;
}

std::vector<std::string> split(const std::string &s, char delimiter = '\t') {
  std::vector<std::string> tokens;
  std::string token;
  std::istringstream tokenStream(s);
  while (std::getline(tokenStream, token, delimiter))
    tokens.push_back(token);
  return tokens;
}

std::vector<std::pair<std::string, bool>> split_path(const std::string &path) {
  std::vector<std::pair<std::string, bool>> p;
  size_t s = 0, e = 1;
  while (e <= path.size()) {
    char c = path[e];
    if (c == '<' || c == '>' || c == '\0') {
      p.push_back(
          std::make_pair(path.substr(s + 1, e - s - 1), path[s] == '>'));
      s = e;
    }
    ++e;
  }
  return p;
}

std::vector<std::pair<int, char>> parse_cigar(const std::string &cigar) {
  std::vector<std::pair<int, char>> cc;
  size_t s = 0, e = 1;
  while (e < cigar.size()) {
    char c = cigar[e];
    if (c == '=' || c == 'X' || c == 'I' || c == 'D') {
      cc.push_back(std::make_pair(std::stoi(cigar.substr(s, e - s)), cigar[e]));
      s = e + 1;
    }
    ++e;
  }
  return cc;
}

std::string build_cs(const std::string &pseq, const std::string &cseq,
                     const std::string &cigar_str) {
  std::vector<std::pair<int, char>> cigar = parse_cigar(cigar_str);

  std::string cs = "";
  int cons_p = 0;
  int pseq_p = 0;
  int opl;
  char op;
  int matches = 0;
  for (size_t i = 0; i < cigar.size(); ++i) {
    opl = cigar[i].first;
    op = cigar[i].second;

    if (op == '=') {
      matches += opl;
      cons_p += opl;
      pseq_p += opl;
    } else {
      if (matches > 0) {
        cs += ":" + std::to_string(matches);
        matches = 0;
      }
      if (op == 'X') {
        for (int j = 0; j < opl; ++j) {
          cs += "*";
          cs += pseq[pseq_p];
          cs += cseq[cons_p];
          pseq_p += 1;
          cons_p += 1;
        }
      } else if (op == 'I') {
        cs += "+" + cseq.substr(cons_p, opl);
        cons_p += opl;
      } else if (op == 'D') {
        cs += "-" + pseq.substr(pseq_p, opl);
        pseq_p += opl;
      } else {
        std::cerr << "Error in CIGAR parsing" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
  }
  if (matches > 0) {
    cs += ":" + std::to_string(matches);
  }

  return cs;
}

std::string join(const std::vector<std::string> &fields,
                 char delimiter = '\t') {
  std::ostringstream oss;
  for (size_t i = 0; i < fields.size(); ++i) {
    oss << fields[i];
    if (i != fields.size() - 1) {
      oss << delimiter;
    }
  }
  return oss.str();
}

int main(int argc, char *argv[]) {
  std::string gfa_fn = argv[1];
  std::string fa_fn = argv[2];
  std::string gaf_fn = argv[3];

  std::map<std::string, std::string> segments;
  std::ifstream gfa(gfa_fn);
  std::string line;
  std::vector<std::string> tokens;
  while (std::getline(gfa, line)) {
    if (line[0] == 'S') {
      tokens = split(line);
      segments[tokens[1]] = tokens[2];
    }
  }
  gfa.close();

  gzFile fp = gzopen(fa_fn.c_str(), "r");
  kseq_t *seq = kseq_init(fp);
  int l;
  std::map<std::string, std::string> sequences;
  while ((l = kseq_read(seq)) >= 0) {
    sequences[seq->name.s] = seq->seq.s;
  }
  kseq_destroy(seq);
  gzclose(fp);

  std::ifstream gaf(gaf_fn);
  // std::string line;
  // std::vector<std::string> tokens;
  std::map<std::string, int> counts;
  while (std::getline(gaf, line)) {
    tokens = split(line);
    std::string qname = tokens[0];

    // size_t ql = std::stoi(tokens[1]);

    size_t qs = std::stoi(tokens[2]), qe = std::stoi(tokens[3]);
    assert(tokens[4] == "+");
    std::string path = tokens[5];
    size_t ps = std::stoi(tokens[7]), pe = std::stoi(tokens[8]);
    std::string cigar(tokens[16], 5);

    std::string cseq = sequences[qname].substr(qs, qe - qs);

    std::vector<std::pair<std::string, bool>> vertices = split_path(path);
    std::string pseq;
    for (const std::pair<std::string, bool> &p : vertices) {
      assert(segments.find(p.first) != segments.end());
      if (p.second)
        pseq += segments[p.first];
      else
        pseq += reverseAndComplement(segments[p.first]);
    }
    pseq = pseq.substr(ps, pe - ps);

    std::string cs = build_cs(pseq, cseq, cigar);
    tokens.push_back("cs:Z:" + cs);
    tokens[0] = "contig." + tokens[0] + "." + std::to_string(counts[qname]);

    std::cout << join(tokens) << std::endl;

    ++counts[qname];
  }
  gaf.close();

  return 0;
}
