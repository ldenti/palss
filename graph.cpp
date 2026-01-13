#include "graph.hpp"

Graph::Graph(const std::string &fn, const std::string &reference)
    : fn(fn), reference(reference) {}

int Graph::load() {
  double rt = realtime();
  sdsl::simple_sds::load_from(gbz, fn);
  fprintf(stderr, "[M::%s] restored GBZ in %.3f sec\n", __func__,
          realtime() - rt);
  return 0;
}
int Graph::load_fl() {
  double rt = realtime();
  std::ifstream in;
  in.open(fn + ".ri", std::ifstream::in);
  fl.load(in);
  fl.setGBWT(gbz.index);
  in.close();
  fprintf(stderr, "[M::%s] restored R-index in %.3f sec\n", __func__,
          realtime() - rt);
  return 0;
}

int Graph::build_fl() {
  double rt = realtime();
  fl = gbwt::FastLocate(gbz.index);
  std::ofstream out;
  out.open(fn + ".ri", std::ofstream::out);
  fl.serialize(out);
  out.close();
  fprintf(stderr, "[M::%s] built R-index in %.3f sec\n", __func__,
          realtime() - rt);
  return 0;
}
const gbwt::Metadata &Graph::get_metadata() const { return gbz.index.metadata; }

const gbwtgraph::GBZ &Graph::get_gbz() const { return gbz; }

const gbwt::FastLocate &Graph::get_fl() const { return fl; }

void Graph::print_stats() const {
  gbwt::printStatistics(gbz.index, "gbwt");
  gbwt::printStatistics(fl, "fl");
}

path_t Graph::get_path(gbwt::size_type path_id) const {
  path_t path;

  gbwt::size_type seq_id = gbwt::Path::encode(path_id, false);
  path.id = seq_id;

  gbwt::edge_type curr = gbz.index.start(seq_id);
  while (curr.first != gbwt::ENDMARKER) {
    path.vertices.push_back(curr.first);
    gbwtgraph::handle_t handle =
        gbwtgraph::GBWTGraph::node_to_handle(curr.first);
    // view_type: in-place view of the sequence: (start, length)
    gbwtgraph::view_type view = gbz.graph.get_sequence_view(handle);
    path.sequence.append(view.first, view.second);
    curr = gbz.index.LF(curr);
  }
  return path;
}

// v1 and v2 has strand information
std::vector<path_t> Graph::get_paths(uint32_t v1, uint32_t v2, uint8_t strand,
                                     bool ref_only) const {
  std::vector<path_t> paths;

  int ev1 = v1;
  std::vector<gbwt::size_type> intervals1 = fl.decompressSA(ev1);
  int ev2 = v2;
  std::vector<gbwt::size_type> intervals2 = fl.decompressSA(ev2);

  // std::cerr << get_gfa_name(v1 >> 1) << " " << gbwt::Node::is_reverse(v1)
  //           << std::endl;

  for (size_t i = 0; i < intervals1.size(); ++i) {
    gbwt::FastLocate::size_type int1 = intervals1[i];
    gbwt::size_type seqid1 = fl.seqId(int1);
    // std::cerr << gbz.index.metadata.fullPath(seqid1 >> 1).sample_name << " "
    //           << gbz.index.metadata.fullPath(seqid1 >> 1).haplotype << " "
    //           << gbz.index.metadata.fullPath(seqid1 >> 1).contig_name << " "
    //           << gbwt::Path::is_reverse(seqid1) << std::endl;

    // CHECKME: is this right?
    if (!(strand & (1 << gbwt::Path::is_reverse(seqid1))))
      continue;
    // if (ev1 == ev2 && gbwt::Path::is_reverse(seqid1))
    //   continue;

    std::string sample_name =
        gbz.index.metadata.fullPath(seqid1 >> 1).sample_name;
    if (ref_only && sample_name.compare(reference) != 0)
      continue;

    gbwt::size_type seqoff1 = fl.seqOffset(int1);
    for (size_t j = 0; j < intervals2.size(); ++j) {
      gbwt::FastLocate::size_type int2 = intervals2[j];
      gbwt::size_type seqid2 = fl.seqId(int2);
      if (seqid1 != seqid2)
        continue;

      gbwt::size_type seqoff2 = fl.seqOffset(int2);

      if (seqoff2 > seqoff1) {
        // we do not have v1 --> v2
        continue;
        // std::swap(position, position2);
      }

      gbwt::edge_type position = std::make_pair(ev1, i);
      gbwt::edge_type position2 = std::make_pair(ev2, j);

      assert(fl.index->contains(position));
      assert(fl.index->contains(position2));
      path_t path;
      path.id = seqid1;
      path.is_reference =
          this->get_path_sample(path.id).compare(this->reference) == 0;
      path.reversed = false;
      path.offset1 = seqoff1;
      path.offset2 = seqoff2;
      while (position.first != position2.first) {
        // position.first is encoded with strand
        path.vertices.push_back(position.first);
        gbwtgraph::handle_t handle =
            gbwtgraph::GBWTGraph::node_to_handle(position.first);
        // view_type: in-place view of the sequence: (start, length)
        gbwtgraph::view_type view = gbz.graph.get_sequence_view(handle);
        // std::cerr << gbz.graph.get_segment_name(handle) << " ";
        path.sequence.append(view.first, view.second);
        position = gbz.index.LF(position);
      }
      path.vertices.push_back(position.first);
      gbwtgraph::handle_t handle =
          gbwtgraph::GBWTGraph::node_to_handle(position.first);
      // std::cerr << gbz.graph.get_segment_name(handle) << std::endl;
      // view_type: in-place view of the sequence: (start, length)
      gbwtgraph::view_type view = gbz.graph.get_sequence_view(handle);
      path.sequence.append(view.first, view.second);
      paths.push_back(path);
    }
  }

  return paths;
}

// std::vector<path_t> Graph::get_paths_both(uint32_t v1, uint32_t v2,
//                                           bool ref_only) const {
//   std::vector<path_t> paths;
//   int strand1 = 0;
//   // for (int strand1 = 0; strand1 < 2; ++strand1) {
//   // XXX: assuming to have the same paths from v1+ to v2- and v1- to v2+ (on
//   // reversed path), same for +>+ and ->-
//   int ev1 = gbwt::Node::encode(v1, strand1);
//   std::vector<gbwt::size_type> intervals1 = fl.decompressSA(ev1);
//   for (int strand2 = 0; strand2 < 2; ++strand2) {
//     int ev2 = gbwt::Node::encode(v2, strand2);
//     std::vector<gbwt::size_type> intervals2 = fl.decompressSA(ev2);
//     for (size_t i = 0; i < intervals1.size(); ++i) {
//       gbwt::FastLocate::size_type int1 = intervals1[i];
//       gbwt::size_type seqid1 = fl.seqId(int1);
//       std::string sample_name =
//           gbz.index.metadata.fullPath(seqid1 >> 1).sample_name;
//       if (ref_only && sample_name.compare(reference) != 0)
//         continue;
//       gbwt::size_type seqoff1 = fl.seqOffset(int1);
//       for (size_t j = 0; j < intervals2.size(); ++j) {
//         gbwt::FastLocate::size_type int2 = intervals2[j];
//         gbwt::size_type seqid2 = fl.seqId(int2);
//         gbwt::size_type seqoff2 = fl.seqOffset(int2);
//         if (seqid1 != seqid2)
//           continue;
//         // sample_name = gbz.index.metadata.fullPath(seqid2).sample_name;
//         // if (ref_only && sample_name.compare(reference) != 0)
//         //   continue;
//         // use seqoff to understand where to start/end
//         gbwt::edge_type position = std::make_pair(ev1, i);
//         gbwt::edge_type position2 = std::make_pair(ev2, j);
//         if (seqoff2 > seqoff1)
//           std::swap(position, position2);
//         assert(fl.index->contains(position));
//         assert(fl.index->contains(position2));
//         path_t path;
//         path.id = seqid1;
//         path.offset1 = seqoff1;
//         path.offset2 = seqoff2;
//         while (position.first != position2.first) {
//           // position.first is encoded with strand
//           path.vertices.push_back(position.first);
//           gbwtgraph::handle_t handle =
//               gbwtgraph::GBWTGraph::node_to_handle(position.first);
//           // view_type: in-place view of the sequence: (start, length)
//           gbwtgraph::view_type view = gbz.graph.get_sequence_view(handle);
//           // std::cerr << gbz.graph.get_segment_name(handle) << " ";
//           path.sequence.append(view.first, view.second);
//           position = gbz.index.LF(position);
//         }
//         path.vertices.push_back(position.first);
//         gbwtgraph::handle_t handle =
//             gbwtgraph::GBWTGraph::node_to_handle(position.first);
//         // std::cerr << gbz.graph.get_segment_name(handle) << std::endl;
//         // view_type: in-place view of the sequence: (start, length)
//         gbwtgraph::view_type view = gbz.graph.get_sequence_view(handle);
//         path.sequence.append(view.first, view.second);
//         paths.push_back(path);
//       }
//     }
//   }
//   // }
//   return paths;
// }

std::string Graph::get_gfa_name(uint32_t v) const {
  gbwtgraph::handle_t h = gbwtgraph::GBWTGraph::node_to_handle(v << 1);
  return gbz.graph.get_segment_name(h);
}

size_t Graph::get_vertex_len(gbwt::size_type v) const {
  gbwtgraph::handle_t h = gbwtgraph::GBWTGraph::node_to_handle(v << 1);
  return gbz.graph.get_length(h);
}

// no strand bit here
std::string Graph::get_path_contig(gbwt::size_type pid) const {
  // std::cerr << gbz.graph.get_path_name(gbz.graph.path_to_handle(pid))
  //           << std::endl;
  // std::cerr << gbz.index.metadata.fullPath(pid).contig_name << std::endl;
  return gbz.index.metadata.fullPath(pid).contig_name;
}

std::string Graph::get_path_sample(gbwt::size_type pid) const {
  return gbz.index.metadata.fullPath(pid).sample_name;
}

std::map<std::string, refpath_t> Graph::get_reference_paths() const {
  std::map<std::string, refpath_t> ret;

  gbwt::size_type sample_id = gbz.index.metadata.sample(reference);

  for (const gbwt::size_type path_id :
       gbz.index.metadata.pathsForSample(sample_id)) {
    std::string contig_name = gbz.index.metadata.fullPath(path_id).contig_name;

    size_t seql = 0;
    gbwt::size_type seq_id = gbwt::Path::encode(path_id, false);
    gbwt::edge_type curr = gbz.index.start(seq_id);
    std::vector<std::pair<gbwt::size_type, size_t>> vertices;
    std::map<gbwt::size_type, size_t> offsets;
    while (curr.first != gbwt::ENDMARKER) {
      gbwtgraph::handle_t handle =
          gbwtgraph::GBWTGraph::node_to_handle(curr.first);
      offsets[curr.first] = seql;
      seql += gbz.graph.get_length(handle);
      curr = gbz.index.LF(curr);
    }

    ret[contig_name] = {path_id, contig_name, offsets, seql};
  }

  return ret;
}

// std::string
// Graph::get_path_sequence(const std::vector<gbwt::node_type> &vertices) const
// {
//   std::string sequence;
//   for (const gbwt::node_type &v : vertices) {
//     gbwtgraph::handle_t h = gbwtgraph::GBWTGraph::node_to_handle(v);
//     sequence += gbz.graph.get_sequence(h);
//   }
//   return sequence;
// }

// === paths

bool operator==(const path_t &path1, const path_t &path2) {
  return path1.vertices == path2.vertices;
}

bool operator<(const path_t &path1, const path_t &path2) {
  return path1.vertices < path2.vertices;
}

void p_reverse(path_t &path) {
  path.id = (path.id >> 1) << 1; // XXX: do we want to change this?
  std::reverse(path.vertices.begin(), path.vertices.end());

  // TODO improve this and make function
  std::string reversed(path.sequence.rbegin(), path.sequence.rend());
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
  path.sequence = complement;
  path.offset1 = path.vertices.size() - path.offset2 - 1;
  path.offset2 = path.vertices.size() - path.offset1 - 1;
}
