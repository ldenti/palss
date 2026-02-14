import sys
import re

PREFIX = "new"


def get_reads_support(gaf_fn):
    print("Parsing GAF...", file=sys.stderr)
    np_reads_support = {}
    for line in open(gaf_fn):
        fields = line.strip("\n").split("\t")
        name = fields[0]
        np_reads_support[name] = set(fields[15].split(":")[-1].split("|"))
    return np_reads_support


def parse_gfa(gfa_fn, np_reads_support):
    print("Parsing GFA...", file=sys.stderr)
    knownV, knownE = set(), set()
    new_paths = {}
    for line in open(gfa_fn):
        path = []
        skip = False
        if line.startswith("P"):
            line = line.strip("\n").split("\t")
            name = line[1]
            path = [int(x[:-1]) for x in line[2].split(",")]
            if name in np_reads_support:
                # XXX: assuming new (augmented) paths are always P lines
                new_paths[name] = path
                skip = True
        elif line.startswith("W"):
            line = line.strip("\n").split("\t")
            path = [int(x) for x in re.split("[<>]", line[6][1:])]
        else:
            continue

        if skip:
            # path is a new path, do not store its vertices/edges
            continue

        knownV.add(path[0])
        for v1, v2 in zip(path[:-1], path[1:]):
            knownV.add(v2)
            e = (min(v1, v2), max(v1, v2))
            knownE.add(e)

    return knownV, knownE, new_paths


def dump_gfa(gfa_fn, np_reads_support, v_to_remove, e_to_remove):
    print("Outputting...", file=sys.stderr)
    for line in open(gfa_fn):
        if line.startswith("S"):
            _, v, seq = line.strip("\n").split("\t")
            v = int(v)
            if v in v_to_remove:
                continue
        elif line.startswith("L"):
            _, v1, _, v2, *_ = line.split("\t")
            v1, v2 = int(v1), int(v2)
            if v1 in v_to_remove or v2 in v_to_remove:
                continue
            e = (min(v1, v2), max(v1, v2))
            if e in e_to_remove:
                continue
        elif line.startswith("P"):
            name = line.strip("\n").split("\t")[1]
            if name in np_reads_support or name.startswith(PREFIX):
                # do not print augmented paths referring to consensus or mapping paths
                continue
        print(line, end="")


def parse_gaf(gaf_fn, span_ratio=0.9):  # XXX: hardcoded
    print(f"Getting best alignments from {gaf_fn}...", file=sys.stderr)
    alns = {}
    for line in open(gaf_fn):
        fields = line.strip("\n").split("\t")

        name, ql, qs, qe, _, path = fields[:6]
        ql, qs, qe = int(ql), int(qs), int(qe)
        if (qe - qs) / ql < span_ratio:
            continue
        c = (qe - qs) / ql

        nm = int(fields[12].split(":")[-1])
        # in case of multimappings, keep the one that covers more bases
        # of the read and, in case of tie, the one with lowest NM
        if name in alns:
            old_c, old_nm, old_cigar, old_path = alns[name]
            if old_c > c:
                continue
            if old_c == c and old_nm <= nm:
                continue

        cigar = fields[16].split(":")[-1]
        alns[name] = (c, nm, cigar, path)
    return alns


def support():
    # augmented graph, with paths referring to each consensus added to the graph
    gfa_fn = sys.argv[1]
    # palss gaf for support
    gaf_fn = sys.argv[2]
    # minimum support to keep (>=)
    minw = int(sys.argv[3])

    np_reads_support = get_reads_support(gaf_fn)

    knownV, knownE, new_paths = parse_gfa(gfa_fn, np_reads_support)

    print("Extracting new segments and links...", file=sys.stderr)
    newV, newE = {}, {}
    for line in open(gfa_fn):
        if line.startswith("S"):
            _, v, seq = line.strip("\n").split("\t")
            v = int(v)
            if v not in knownV:
                newV[v] = set()
        elif line.startswith("L"):
            _, v1, _, v2, _, *_ = line.strip("\n").split("\t")
            v1, v2 = int(v1), int(v2)
            e = (min(v1, v2), max(v1, v2))
            if e not in knownE and v1 in knownV and v2 in knownV:
                newE[e] = set()

    print(
        f"Selecting segments/links to remove (due to low support, threshold: {minw})...",
        file=sys.stderr,
    )
    for name, path in new_paths.items():
        if path[0] in newV:
            newV[path[0]].add(name)
        for v1, v2 in zip(path[:-1], path[1:]):
            if v2 in newV:
                newV[v2].add(name)  # just store consensus

            e = (min(v1, v2), max(v1, v2))
            if e in newE:
                newE[e].add(name)

    v_to_remove = set()
    for v, Cs in newV.items():
        reads = set().union(*[np_reads_support[c] for c in Cs])
        if len(reads) < minw:
            v_to_remove.add(v)
    print(
        f"Need to remove {len(v_to_remove)} out of {len(newV)} novel vertices",
        file=sys.stderr,
    )

    e_to_remove = set()
    for e, Cs in newE.items():
        reads = set().union(*[np_reads_support[c] for c in Cs])
        if len(reads) < minw:
            e_to_remove.add(e)
    print(
        f"Need to remove {len(e_to_remove)} out of {len(newE)} novel edges",
        file=sys.stderr,
    )

    dump_gfa(gfa_fn, np_reads_support, v_to_remove, e_to_remove)

    for i, (v, Cs) in enumerate(newV.items()):
        if v not in v_to_remove:
            print(
                "P",
                f"{PREFIX}.v{i}",
                f"{v}+",
                "*",
                "CS:Z:" + ",".join([c for c in Cs]),
                sep="\t",
            )
    for i, (e, Cs) in enumerate(newE.items()):
        if e not in e_to_remove:
            v1, v2 = e
            print(
                "P",
                f"{PREFIX}.e{i}",
                f"{v1}+,{v2}+",
                "*",
                "CS:Z:" + ",".join([c for c in Cs]),
                sep="\t",
            )


def extract_mapping(gfa0_fn, gfa_fn):
    """
    While cleaning the graph (first-pass), we stored as P lines with
    name prefixed by PREFIX some sort of mapping between vertices/edges
    and palss consensuses. This to have a way to map new elements from
    unchopped graph back to consensuses (we needed to unchop the graph
    to facilitate read alignment). Unchopping changes the ID space.
    The path names are the only things that is maintained between the
    two graphs.
    """

    # iterate over chopped/original clean graph
    p2v = {}  # path to old vertex
    p2e = {}  # path to old edge
    p2cs = {}  # path to consensuses
    for line in open(gfa0_fn):
        if line[0] == "P":
            fields = line.strip("\n").split("\t")
            pname = fields[1]
            if pname.startswith(PREFIX):
                p2cs[pname] = fields[4].split(":")[-1].split(",")
                if ".v" in pname:
                    p2v[pname] = int(fields[2][:-1])
                else:
                    p2e[pname] = (int(x[:-1]) for x in fields[2].split(","))

    # iterate over new unchopped clean graph
    v_old2new = {}
    e_old2new = {}
    newv2cs = {}
    newe2cs = {}
    for line in open(gfa_fn):
        if line[0] == "P":
            fields = line.strip("\n").split("\t")
            pname = fields[1]
            if pname.startswith(PREFIX):
                if ".v" in pname:
                    v = int(fields[2][:-1])
                    # v is the new vertex ID
                    # p2v[pname] is the old vertex ID
                    newv2cs[v] = p2cs[pname]
                    v_old2new[p2v[pname]] = v
                else:
                    # v1 and v2 are the new vertices
                    # p2e[pname] are the old vertices
                    v1, v2 = (int(x[:-1]) for x in fields[2].split(","))
                    newe2cs[(v1, v2)] = p2cs[pname]
                    e_old2new[p2e[pname]] = (v1, v2)
    return v_old2new, e_old2new, newv2cs, newe2cs


def post():
    gfa0_fn = sys.argv[1]  # augmented graph, clean/original, with mapping
    gfa_fn = sys.argv[2]  # unchopped augmented graph (alignments refer to this graph)
    palss_gaf_fn = sys.argv[3]  # consensus from palss (to get reads support)
    original_gaf_fn = sys.argv[4]  # consensus aligned to original graph
    augmented_gaf_fn = sys.argv[5]  # consensus aligned tounchopped augmented graph
    minw = int(sys.argv[6])  # minimum support, remove all <
    resulting_c = sys.argv[7]

    np_reads_support = get_reads_support(palss_gaf_fn)

    v_old2new, e_old2new, newv2cons, newe2cons = extract_mapping(gfa0_fn, gfa_fn)
    # NOTE: v_old2new, e_old2new are currently unused

    # invert maps
    cons2newv = {}
    for v, Cs in newv2cons.items():
        for c in Cs:
            if c not in cons2newv:
                cons2newv[c] = set()
            cons2newv[c].add(v)
    cons2newe = {}
    for e, Cs in newe2cons.items():
        for c in Cs:
            if c not in cons2newe:
                cons2newe[c] = set()
            cons2newe[c].add(e)

    # Get alignments
    orig_alns = parse_gaf(original_gaf_fn)
    aug_alns = parse_gaf(augmented_gaf_fn)

    print("Original alignments:", len(orig_alns), file=sys.stderr)
    print("Augmented alignments:", len(aug_alns), file=sys.stderr)
    print("Using", len(set(orig_alns) & set(aug_alns)), "alignments", file=sys.stderr)

    # Save those consensus aligned better to augmented graph (c_tokeep)
    # and get vertices/edges supported by the alignments
    c_tokeep = set()
    supported_vertices = set()
    supported_edges = set()
    info = [0, 0, 0, 0]  # some statistics
    for qname in set(orig_alns) & set(aug_alns):
        oc, onm, ocigar, opath = orig_alns[qname]
        ac, anm, acigar, apath = aug_alns[qname]
        if onm == 0 or onm == anm:
            # no improvement, so maybe we should not consider the corresponding consensus
            i = 0
        elif onm < anm:
            # worse alignment to augmented, why? Some problems with graphaligner?
            # We won't consider this consensus
            print(qname, onm, oc, ocigar, anm, ac, acigar, file=sys.stderr)
            i = 1
        else:
            # Better alignment to augmented, store everything we need
            i = 2
            c_tokeep.add(qname)
            apath = re.split("[<>]", apath[1:])
            supported_vertices |= set([int(x) for x in apath])
            supported_edges |= set(
                [
                    (min(int(x1), int(x2)), max(int(x1), int(x2)))
                    for x1, x2 in zip(apath[:-1], apath[1:])
                ]
            )
        info[i] += 1

    # Save all consensus that were not aligned to original graph but to the augmented one
    for qname in set(aug_alns) - set(orig_alns):
        info[3] += 1
        c_tokeep.add(qname)
        path = re.split("[<>]", aug_alns[qname][3][1:])
        supported_vertices |= set([int(x) for x in path])
        supported_edges |= set(
            [
                (min(int(x1), int(x2)), max(int(x1), int(x2)))
                for x1, x2 in zip(path[:-1], path[1:])
            ]
        )

    print("Supported vertices:", len(supported_vertices), file=sys.stderr)
    print("Supported edges:", len(supported_edges), file=sys.stderr)
    print("Info:", info, file=sys.stderr)

    # Iterate over new vertices and retrieve read support.
    # If not enough, put the vertex in the toremove set
    v_to_remove = set()
    for v, Cs in newv2cons.items():
        if v not in supported_vertices:
            # we do not have any consensus aligned to this vertex
            v_to_remove.add(v)
            continue
        reads = set()
        for c in Cs:
            if c in c_tokeep:
                reads |= np_reads_support[c]
        if len(reads) < minw:
            # vertex is not supported by minw reads (after consensus selection)
            v_to_remove.add(v)
    print(
        f"Need to remove {len(v_to_remove)} out of {len(newv2cons)} novel vertices",
        file=sys.stderr,
    )

    e_to_remove = set()
    for e, Cs in newe2cons.items():
        if e not in supported_edges:
            # we do not have any consensus aligned to this edge
            e_to_remove.add(e)
            continue
        reads = set()
        for c in Cs:
            if c in c_tokeep:
                reads |= np_reads_support[c]
        if len(reads) < minw:
            # edge is not supported by minw reads (after consensus selection)
            e_to_remove.add(e)

    print(
        f"Need to remove {len(e_to_remove)} out of {len(newe2cons)} novel edges",
        file=sys.stderr,
    )

    dump_gfa(gfa_fn, np_reads_support, v_to_remove, e_to_remove)

    # Get vertices information
    vseqs = {}
    for line in open(gfa_fn):
        if line.startswith("S"):
            _, v, seq, *_ = line.strip("\n").split("\t")
            v = int(v)
            if v in newv2cons and v not in v_to_remove:
                vseqs[v] = seq

    # Select and dump consensus supporting novel vertices and/or edges we kept
    c_tokeep = set()
    for v, Cs in newv2cons.items():
        if v not in v_to_remove:
            c_tokeep |= set(newv2cons[v])
    for e, Cs in newe2cons.items():
        if e not in e_to_remove:
            c_tokeep |= set(newe2cons[e])

    with open(resulting_c, "w") as rc:
        for line in open(palss_gaf_fn):
            fields = line.split("\t")
            name = fields[0]
            if name in c_tokeep:
                new_vertices = []
                new_edges = []
                if name in cons2newv:
                    for v in cons2newv[name]:
                        if v not in v_to_remove:
                            new_vertices.append(v)
                if name in cons2newe:
                    for e in cons2newe[name]:
                        if e not in e_to_remove:
                            new_edges.append(e)

                print(line.strip("\n"), end="", file=rc)
                if len(new_vertices) > 0:
                    print(
                        f'\tNV:Z:{"/".join([f"{v}|{vseqs[v]}" for v in new_vertices])}',
                        end="",
                        file=rc,
                    )
                if len(new_edges) > 0:
                    print(
                        f'\tNE:Z:{"/".join([f"{v1}|{v2}" for v1, v2 in new_edges])}',
                        end="",
                        file=rc,
                    )
                print("", file=rc)


if __name__ == "__main__":
    mode = sys.argv.pop(1)
    if mode == "1":
        support()
    elif mode == "2":
        post()
