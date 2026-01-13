#include <getopt.h>
#include <iostream>
#include <string>
#include <vector>

#include "graph.hpp"
#include "sfs.hpp"

std::string rc(const std::string &seq) {
  std::string reversed(seq.rbegin(), seq.rend());
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
      complement += nucleotide; // handles non-nucleotide characters
      break;
    }
  }
  return complement;
}

int main_sam(int argc, char *argv[]) {
  char *gbz_fn = argv[1];
  char *sfs_fn = argv[2];

  std::string ref_path = "CHM13";

  size_t klen = 31;

  Graph graph(gbz_fn, ref_path);
  graph.load();
  graph.load_fl();

  const gbwt::Metadata &metadata = graph.get_metadata();
  const gbwtgraph::GBZ &gbz = graph.get_gbz();

  // std::unordered_set<std::string> ref_paths = gbz.get_reference_samples();
  // for (const auto &r : ref_paths)
  //   std::cerr << r << std::endl;

  std::map<std::string, size_t> ref_lengths;
  std::map<gbwt::size_type, std::map<gbwt::short_type, size_t>> offsets;
  gbwt::size_type sample_id = metadata.sample(ref_path);
  for (gbwt::size_type &path_id : metadata.pathsForSample(sample_id)) {
    gbwt::size_type seq_id = gbwt::Path::encode(path_id, false);
    gbwt::edge_type curr = gbz.index.start(seq_id);
    size_t offset = 0;
    while (curr.first != gbwt::ENDMARKER) {
      gbwtgraph::handle_t handle =
          gbwtgraph::GBWTGraph::node_to_handle(curr.first);
      size_t l = gbz.graph.get_length(handle);
      offsets[path_id][curr.first] = offset;
      offset += l;
      curr = gbz.index.LF(curr);
    }
    ref_lengths[metadata.contig(path_id)] = offset;
  }

  std::cout << "@HD\tVN:1.6\tSO:coordinate" << std::endl;
  for (const auto &[k, v] : ref_lengths) {
    std::cout << "@SQ" << "\t"
              << "SN:" << k << "\t"
              << "LN:" << v << std::endl;
  }

  std::ifstream fp(sfs_fn);
  std::string line;
  if (!fp.is_open()) {
    fprintf(stderr, "Unable to open SFS file!");
    return 1;
  }

  while (getline(fp, line)) {
    if (line[0] != '0')
      continue;

    sfs_t s = parse_sfs_line(line);

    size_t path_id;
    bool on_ref = false;
    bool strand = false;
    for (const uint64_t &p : s.paths) {
      if (gbz.index.metadata.fullPath(p >> 2).sample_name.compare(ref_path) ==
          0) {
        path_id = p >> 2;
        strand = p & 1;
        on_ref = true;
        break;
      }
    }

    if (!on_ref)
      continue;

    std::string qidx =
        s.rname + "_" + std::to_string(s.s) + "_" + std::to_string(s.s + s.l);
    std::cerr << qidx << std::endl;
    size_t pos1 = offsets[path_id][s.sv];
    size_t pos2 = offsets[path_id][s.ev];
    size_t reference_start = pos1 + s.soff;
    size_t reference_end = pos2 + s.eoff + 1;
    if (pos1 > pos2) {
      reference_start = pos2 + s.eoff;
      reference_end = pos1 + s.soff + 1;
    }
    int nn = reference_end - reference_start - 2 * klen;
    // NOTE: we can have negative nn if anchors are overlapping on graph but not
    // on read (due to insertion)
    std::string seq;
    std::string cigar;
    if (nn > 0) {
      cigar = std::to_string(klen) + "M" + std::to_string(nn) + "N" +
              std::to_string(klen) + "M";
      seq = s.plain_seq.substr(0, klen) + s.plain_seq.substr(s.l - klen, klen);
    } else {
      nn = -nn;
      cigar = std::to_string(2 * klen - nn) + "M";
      seq = s.plain_seq.substr(0, klen - nn) + s.plain_seq.substr(s.l - klen, klen);
    }

    if (!strand) {
      seq = rc(seq);
    }

    std::cout << qidx << "\t" << (strand ? "0" : "16") << "\t"
              << gbz.index.metadata.fullPath(path_id).contig_name << "\t"
              << reference_start + 1 << "\t" << 60 << "\t" << cigar << "\t"
              << "*" << "\t" << 0 << "\t" << 0 << "\t" << seq << "\t"
              << "*" << std::endl;
  }
  fp.close();

  return 0;
}
