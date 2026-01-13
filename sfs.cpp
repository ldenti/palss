#include "sfs.hpp"

std::vector<std::string> split(const std::string &s, char delimiter) {
  std::vector<std::string> tokens;
  std::string token;
  std::istringstream tokenStream(s);
  while (std::getline(tokenStream, token, delimiter))
    tokens.push_back(token);
  return tokens;
}

// std::map<std::string, std::vector<sfs_t>> load_sfs(const std::string &fn) {
//   std::ifstream fp(fn);
//   std::string line;
//   std::map<std::string, std::vector<sfs_t>> specific_strings;
//   if (!fp.is_open()) {
//     fprintf(stderr, "Unable to open SFS file!");
//     return specific_strings;
//   }
//   while (getline(fp, line)) {
//     if (line[0] != '0')
//       continue;

//     std::vector<std::string> tokens = split(line, '\t');

//     sfs_t s;
//     s.flag = 0;
//     s.rname = tokens[1];
//     s.s = std::stoi(tokens[2]);
//     s.l = std::stoi(tokens[3]);
//     // s.end = std::stoi(tokens[4]);
//     // if (tokens.size() > 5) {
//     //   // anchored specific strings
//     //   s.sv1 = std::stoi(tokens[5]);
//     //   s.sv2 = std::stoi(tokens[6]);
//     //   s.ev1 = std::stoi(tokens[7]);
//     //   s.ev2 = std::stoi(tokens[8]);
//     //   // assuming skmer <= ekmer
//     //   s.skmer = std::stoul(tokens[9]);
//     //   s.ekmer = std::stoul(tokens[10]);
//     //   //
//     //   s.plain_seq = tokens[12];
//     //   s.seq = (uint8_t *)malloc(s.l + 1);
//     //   for (int i = 0; i < s.l; ++i)
//     //     s.seq[i] = tokens[12][i] < 128 ? to_int[(int)tokens[12][i]] - 1 :
//     4;
//     //   s.seq[s.l] = '\0';
//     // }
//     specific_strings[s.rname].push_back(s);
//   }
//   fp.close();

//   return specific_strings;
// }

sfs_t parse_sfs_line(const std::string &line) {
  std::vector<std::string> tokens = split(line, '\t');

  sfs_t s;
  s.flag = 0;
  s.rname = tokens[1];
  s.s = std::stoi(tokens[2]);
  s.l = std::stoi(tokens[3]);
  // s.end = std::stoi(tokens[4]);
  //
  s.sv = std::stoi(tokens[5]);
  s.ev = std::stoi(tokens[6]);
  //
  s.soff = std::stoi(tokens[7]);
  s.eoff = std::stoi(tokens[8]);
  //
  // GFA identifiers, 9/10
  //
  s.skmer = std::stoul(tokens[11]);
  s.ekmer = std::stoul(tokens[12]);
  //
  std::vector<std::string> paths = split(tokens[13], ',');
  for (const std::string &p : paths)
    s.paths.push_back(std::stoi(p));
  //
  s.plain_seq = tokens[14];
  s.seq = (uint8_t *)malloc(s.l + 1);
  for (int i = 0; i < s.l; ++i)
    s.seq[i] = s.plain_seq[i] < 128 ? to_int[(int)s.plain_seq[i]] - 1 : 4;
  s.seq[s.l] = '\0';

  return s;
}

std::vector<sfs_t> load_sfs(const std::string &fn) {
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
    sfs_t s = parse_sfs_line(line);
    specific_strings.push_back(s);
  }
  fp.close();

  return specific_strings;
}
