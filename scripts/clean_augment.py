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
    knownV = set()
    new_paths = {}
    for line in open(gfa_fn):
        path = []
        skip = False
        if line.startswith("P"):
            line = line.strip("\n").split("\t")
            name = line[1]
            path = [int(x[:-1]) for x in line[2].split(",")]
            if name in np_reads_support:
                new_paths[name] = path
                skip = True
        elif line.startswith("W"):
            line = line.strip("\n").split("\t")
            path = [int(x) for x in re.split("[<>]", line[6][1:])]
        else:
            continue
        if skip:
            # path is a new path, do not add its vertices to knownV
            continue
        for v in path:
            knownV.add(v)
    return knownV, new_paths


def dump_gfa(gfa_fn, np_reads_support, to_remove):
    print("Outputting...", file=sys.stderr)
    for line in open(gfa_fn):
        if line.startswith("S"):
            _, v, seq = line.strip("\n").split("\t")
            v = int(v)
            if v in to_remove:
                continue
        elif line.startswith("L"):
            _, v1, _, v2, *_ = line.split("\t")
            v1, v2 = int(v1), int(v2)
            if v1 in to_remove or v2 in to_remove:
                continue
        elif line.startswith("P"):
            name = line.strip("\n").split("\t")[1]
            if name in np_reads_support or name.startswith(PREFIX):
                # do not print augmented paths referring to consensus
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
    # # mapping new vertex to consensus
    # txt_fn = sys.argv[4]

    np_reads_support = get_reads_support(gaf_fn)

    knownV, new_paths = parse_gfa(gfa_fn, np_reads_support)

    print("Extracting new segments...", file=sys.stderr)
    newV = {}
    for line in open(gfa_fn):
        if line.startswith("S"):
            _, v, seq = line.strip("\n").split("\t")
            v = int(v)
            if v not in knownV:
                newV[v] = set()

    print(
        f"Selecting segments to remove (due to low support, threshold: {minw})...",
        file=sys.stderr,
    )
    for name, path in new_paths.items():
        for v in path:
            if v in newV:
                newV[v].add(name)  # just store consensus

    to_remove = set()
    for v, Cs in newV.items():
        reads = set().union(*[np_reads_support[c] for c in Cs])
        if len(reads) < minw:
            to_remove.add(v)

    print(
        f"Need to remove {len(to_remove)} out of {len(newV)} novel vertices",
        file=sys.stderr,
    )
    dump_gfa(gfa_fn, np_reads_support, to_remove)
    # mapping = open(txt_fn, "w")
    for i, (v, Cs) in enumerate(newV.items()):
        if v not in to_remove:
            print(
                "P",
                f"{PREFIX}{i}",
                f"{v}+",
                "*",
                "CS:Z:" + ",".join([c for c in Cs]),
                sep="\t",
            )


def extract_mapping(gfa0_fn, gfa_fn):
    p2v = {}  # path to vertex
    p2cs = {}  # path to consensuses
    for line in open(gfa0_fn):
        if line[0] == "P":
            fields = line.strip("\n").split("\t")
            pname = fields[1]
            if pname.startswith(PREFIX):
                p2cs[pname] = fields[4].split(":")[-1].split(",")
                p2v[pname] = int(fields[2][:-1])

    old2new = {}
    new2cs = {}
    for line in open(gfa_fn):
        if line[0] == "P":
            fields = line.strip("\n").split("\t")
            pname = fields[1]
            if pname.startswith(PREFIX):
                v = int(fields[2][:-1])
                new2cs[v] = p2cs[pname]
                old2new[p2v[pname]] = v
    return old2new, new2cs


def post():
    gfa0_fn = sys.argv[1]  # augmented graph, original, with mapping
    gfa_fn = sys.argv[2]  # unchopped augmented graph (alignments refer to this graph)
    palss_gaf_fn = sys.argv[3]  # consensus from palss (to get reads support)
    original_gaf_fn = sys.argv[4]  # consensus aligned to original graph
    augmented_gaf_fn = sys.argv[5]  # consensus aligned tounchopped augmented graph
    minw = int(sys.argv[6])  # minimum support, remove all <
    resulting_c = sys.argv[7]

    np_reads_support = get_reads_support(palss_gaf_fn)
    old2new, newv2cons = extract_mapping(gfa0_fn, gfa_fn)

    cons2newv = {}
    for v, Cs in newv2cons.items():
        for c in Cs:
            if c not in cons2newv:
                cons2newv[c] = set()
            cons2newv[c].add(v)

    orig_alns = parse_gaf(original_gaf_fn)
    aug_alns = parse_gaf(augmented_gaf_fn)

    print("Original alignments:", len(orig_alns), file=sys.stderr)
    print("Augmented alignments:", len(aug_alns), file=sys.stderr)
    print("Using", len(set(orig_alns) & set(aug_alns)), "alignments", file=sys.stderr)

    # Save those consensus aligned better to augmented graph (c_tokeep) and get vertices (supported_vertices)
    supported_vertices = set()
    c_tokeep = set()
    info = [0, 0, 0, 0]
    for qname in set(orig_alns) & set(aug_alns):
        oc, onm, ocigar, opath = orig_alns[qname]
        ac, anm, acigar, apath = aug_alns[qname]
        if onm == 0 or onm == anm:
            i = 0
        elif onm < anm:
            print(qname, onm, oc, ocigar, anm, ac, acigar, file=sys.stderr)
            i = 1
        else:
            i = 2
            c_tokeep.add(qname)
            supported_vertices |= set([int(x) for x in re.split("[<>]", apath[1:])])
        info[i] += 1

    # Save all consensus not aligned to original graph
    for qname in set(aug_alns) - set(orig_alns):
        info[3] += 1
        c_tokeep.add(qname)
        supported_vertices |= set(
            [int(x) for x in re.split("[<>]", aug_alns[qname][3][1:])]
        )

    print("Supported vertices:", len(supported_vertices), file=sys.stderr)
    print("Info:", info, file=sys.stderr)

    to_remove = set()
    for v, Cs in newv2cons.items():
        if v not in supported_vertices:
            # we do not have any consensus aligned to this vertex
            to_remove.add(v)
            continue
        reads = set()
        for c in Cs:
            if c in c_tokeep:
                reads |= np_reads_support[c]
        if len(reads) < minw:
            # vertex is not supported by minw reads (after consensus selection)
            to_remove.add(v)

    print(
        f"Need to remove {len(to_remove)} out of {len(newv2cons)} novel vertices",
        file=sys.stderr,
    )

    dump_gfa(gfa_fn, np_reads_support, to_remove)

    # Get vertices information
    vseqs = {}
    for line in open(gfa_fn):
        if line.startswith("S"):
            _, v, seq, *_ = line.strip("\n").split("\t")
            v = int(v)
            if v in newv2cons and v not in to_remove:
                vseqs[v] = seq

    # Select and dump consensus supporting novel vertices we have kept
    c_tokeep = set()
    for v, cs in newv2cons.items():
        if v not in to_remove:
            c_tokeep |= set(newv2cons[v])
    with open(resulting_c, "w") as rc:
        for line in open(palss_gaf_fn):
            fields = line.split("\t")
            name = fields[0]
            if name in c_tokeep:
                new_vertices = []
                for v in cons2newv[name]:
                    if v not in to_remove:
                        new_vertices.append(v)
                print(
                    line.strip("\n"),
                    f'NV:Z:{"/".join([f"{v}|{vseqs[v]}" for v in new_vertices])}',
                    sep="\t",
                    file=rc,
                )


if __name__ == "__main__":
    mode = sys.argv.pop(1)
    if mode == "1":
        support()
    elif mode == "2":
        post()
