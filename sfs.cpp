#include "sfs.hpp"

std::vector<std::string> split(const std::string &s, char delimiter) {
  std::vector<std::string> tokens;
  std::string token;
  std::istringstream tokenStream(s);
  while (std::getline(tokenStream, token, delimiter))
    tokens.push_back(token);
  return tokens;
}

std::map<std::string, std::vector<sfs_t>> load_sfs(const std::string &fn) {
  std::ifstream fp(fn);
  std::string line;
  std::map<std::string, std::vector<sfs_t>> specific_strings;
  if (!fp.is_open()) {
    fprintf(stderr, "Unable to open SFS file!");
    return specific_strings;
  }
  while (getline(fp, line)) {
    if (line[0] != '0')
      continue;

    std::vector<std::string> tokens = split(line, '\t');

    sfs_t s;
    s.flag = 0;
    s.rname = tokens[1];
    s.s = std::stoi(tokens[2]);
    s.l = std::stoi(tokens[3]);
    // s.end = std::stoi(tokens[4]);
    // if (tokens.size() > 5) {
    //   // anchored specific strings
    //   s.sv1 = std::stoi(tokens[5]);
    //   s.sv2 = std::stoi(tokens[6]);
    //   s.ev1 = std::stoi(tokens[7]);
    //   s.ev2 = std::stoi(tokens[8]);
    //   // assuming skmer <= ekmer
    //   s.skmer = std::stoul(tokens[9]);
    //   s.ekmer = std::stoul(tokens[10]);
    //   //
    //   s.plain_seq = tokens[12];
    //   s.seq = (uint8_t *)malloc(s.l + 1);
    //   for (int i = 0; i < s.l; ++i)
    //     s.seq[i] = tokens[12][i] < 128 ? to_int[(int)tokens[12][i]] - 1 : 4;
    //   s.seq[s.l] = '\0';
    // }
    specific_strings[s.rname].push_back(s);
  }
  fp.close();

  return specific_strings;
}

std::vector<sfs_t> load_anchored_sfs(const std::string &fn) {
  std::ifstream fp(fn);
  std::string line;
  std::vector<sfs_t> specific_strings;
  if (!fp.is_open()) {
    fprintf(stderr, "Unable to open SFS file!");
    return specific_strings;
  }
  while (getline(fp, line)) {
    if (line[0] != '0')
      continue;

    std::vector<std::string> tokens = split(line, '\t');

    sfs_t s;
    s.flag = 0;
    s.rname = tokens[1];
    s.s = std::stoi(tokens[2]);
    s.l = std::stoi(tokens[3]);
    // s.end = std::stoi(tokens[4]);
    if (tokens.size() > 5) {
      // anchored specific strings
      s.sv1 = std::stoi(tokens[5]);
      s.sv2 = std::stoi(tokens[6]);
      s.ev1 = std::stoi(tokens[7]);
      s.ev2 = std::stoi(tokens[8]);
      // assuming skmer <= ekmer
      s.skmer = std::stoul(tokens[9]);
      s.ekmer = std::stoul(tokens[10]);
      //
      s.plain_seq = tokens[15];
      s.seq = (uint8_t *)malloc(s.l + 1);
      for (int i = 0; i < s.l; ++i)
        s.seq[i] = tokens[15][i] < 128 ? to_int[(int)tokens[15][i]] - 1 : 4;
      s.seq[s.l] = '\0';
    }
    specific_strings.push_back(s);
  }
  fp.close();

  return specific_strings;
}