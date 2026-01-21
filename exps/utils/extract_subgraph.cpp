#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "gbwtgraph/gbz.h"

int main(int argc, char *argv[]) {
  std::string gbz_fn = argv[1];
  std::string samples_fn = argv[2];

  std::vector<std::string> samples;
  std::ifstream file(samples_fn);
  std::string line;
  while (std::getline(file, line))
    samples.push_back(line);
  file.close();

  gbwtgraph::GBZ gbz;
  sdsl::simple_sds::load_from(gbz, gbz_fn);

  const gbwt::Metadata &metadata = gbz.index.metadata;

  std::map<gbwt::size_type, size_t> path_lengths;

  std::set<gbwt::node_type> nodes;
  std::set<gbwt::edge_type> edges;

  for (const std::string &sample : samples) {
    gbwt::size_type sample_id = metadata.sample(sample);
    std::vector<gbwt::size_type> paths = metadata.pathsForSample(sample_id);
    for (size_t pp = 0; pp < paths.size(); ++pp) {
      gbwt::size_type path_id = paths[pp];
      gbwt::size_type seq_id = gbwt::Path::encode(path_id, false);

      size_t l = 0;
      gbwt::edge_type curr = gbz.index.start(seq_id);
      l += gbz.graph.get_length(
          gbwtgraph::GBWTGraph::node_to_handle(curr.first));
      nodes.insert(curr.first);
      gbwt::node_type last = curr.first;
      curr = gbz.index.LF(curr);
      while (curr.first != gbwt::ENDMARKER) {
        nodes.insert(curr.first);
        l += gbz.graph.get_length(
            gbwtgraph::GBWTGraph::node_to_handle(curr.first));
        edges.insert(std::make_pair(last, curr.first));
        last = curr.first;
        curr = gbz.index.LF(curr);
      }
      path_lengths[path_id] = l;
    }
  }

  std::cerr << "Extracting " << path_lengths.size() << " paths" << std::endl;

  // XXX: this could be done way better (no need to store strings)

  // Print S lines
  std::set<std::string> segments;
  for (const auto &node : nodes) {
    gbwtgraph::handle_t handle = gbz.graph.node_to_handle(node);
    std::string gfa_name = gbz.graph.get_segment_name(handle);
    if (segments.find(gfa_name) != segments.end())
      continue;
    std::cout << "S" << "\t" << gfa_name << "\t"
              << gbz.graph.get_sequence(handle) << std::endl;
    segments.insert(gfa_name);
  }

  // Print L lines
  std::set<std::string> links;
  for (const auto &edge : edges) {
    gbwt::node_type from = edge.first;
    gbwt::node_type to = edge.second;
    std::string from_gfa_name =
        gbz.graph.get_segment_name(gbz.graph.node_to_handle(from));
    std::string to_gfa_name =
        gbz.graph.get_segment_name(gbz.graph.node_to_handle(to));

    if (from_gfa_name.compare(to_gfa_name) == 0)
      continue;

    std::string line = "L\t" + from_gfa_name + "\t" +
                       (gbwt::Node::is_reverse(from) ? '-' : '+') + "\t" +
                       to_gfa_name + "\t" +
                       (gbwt::Node::is_reverse(to) ? '-' : '+') + "\t" + "0M";

    if (links.find(line) != links.end())
      continue;
    std::cout << line << std::endl;
    links.insert(line);
  }

  // Print W lines
  for (const std::string &sample : samples) {
    gbwt::size_type sample_id = metadata.sample(sample);
    std::vector<gbwt::size_type> paths = metadata.pathsForSample(sample_id);
    for (size_t pp = 0; pp < paths.size(); ++pp) {
      gbwt::size_type path_id = paths[pp];

      gbwt::FullPathName fpn = metadata.fullPath(path_id);

      std::cout << "W" << "\t" << fpn.sample_name << "\t";
      size_t haplotype =
          (fpn.haplotype == gbwtgraph::GBWTGraph::NO_PHASE ? 0 : fpn.haplotype);
      std::cout << haplotype << "\t" << fpn.contig_name << "\t" << fpn.offset
                << "\t" << fpn.offset + path_lengths[path_id] << "\t";

      gbwt::size_type seq_id = gbwt::Path::encode(path_id, false);
      gbwt::edge_type curr = gbz.index.start(seq_id);
      while (curr.first != gbwt::ENDMARKER) {
        std::cout << (gbwt::Node::is_reverse(curr.first) ? '<' : '>')
                  << gbz.graph.get_segment_name(
                         gbz.graph.node_to_handle(curr.first));
        curr = gbz.index.LF(curr);
      }
      std::cout << std::endl;
    }
  }

  return 0;
}
