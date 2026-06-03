import sys
import re


def parse_gaf(gaf_fn, span_ratio=0.9):  # XXX: hardcoded
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


def main():
    gfa_fn = sys.argv[1] # augmented graph
    gaf_fn = sys.argv[2]  # consensus from palss (to get support)
    original_gaf_fn = sys.argv[3]  # consensus aligned to original graph
    augmented_gaf_fn = sys.argv[4]  # consensus aligned to augmented graph
    minw = int(sys.argv[5]) # minimum support, remove all <

    print("Getting read support for consensus...", file=sys.stderr)
    np_reads_support = {} # consensus to set of reads used to build it
    for line in open(gaf_fn):
        fields = line.strip("\n").split("\t")
        name = fields[0]
        np_reads_support[name] = set(fields[15].split(":")[-1].split("|"))

    print("Getting best alignments from original GAF file...", file=sys.stderr)
    orig_alns = parse_gaf(original_gaf_fn)

    print("Getting best alignments from augmented GAF file...", file=sys.stderr)
    aug_alns = parse_gaf(augmented_gaf_fn)

    print("Original alignments:", len(orig_alns), file=sys.stderr)
    print("Augmented alignments:", len(aug_alns), file=sys.stderr)
    print("Using", len(set(orig_alns) & set(aug_alns)), "alignments", file=sys.stderr)

    # Save those consensus aligned better to augmented graph (tokeep) and get vertices (supported_vertices)
    supported_vertices = set()
    tokeep = set()
    info = [0, 0, 0, 0]
    for qname in set(orig_alns) & set(aug_alns):
        oc, onm, ocigar, opath = orig_alns[qname]
        ac, anm, acigar, apath = aug_alns[qname]
        if onm == 0 or onm == anm:
            i = 0
        elif onm < anm:
            i = 1
        else:
            i = 2
            tokeep.add(qname)
            supported_vertices |= set([int(x) for x in re.split("[<>]", apath[1:])])
        info[i] += 1

    # Save all consensus not aligned to original graph
    for qname in set(aug_alns) - set(orig_alns):
        info[3] += 1
        tokeep.add(qname)
        supported_vertices |= set([int(x) for x in re.split("[<>]", aug_alns[qname][3][1:])])

    print("Supported vertices:", len(supported_vertices), file=sys.stderr)
    print("Info:", info, file=sys.stderr)

    # Iterate over graph paths. If path is old, store vertices as known. If path is new, store it in the dict
    print("Parsing paths...", file=sys.stderr)
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

    print("New paths:", len(new_paths), file=sys.stderr)

    # Iterate over new paths and get novel vertices
    nn = 0
    newV = {}
    to_remove = set()
    for name, path in new_paths.items():
        for v in path:
            if v not in knownV:
                # novel vertex
                nn += 1
                # keep the vertex if we have to keep the consensus and v is supported (we might have some vertices not supported by "realignment")
                if name in tokeep and v in supported_vertices:
                    if v not in newV:
                        newV[v] = set()
                    # for each novel vertex, store set of reads (union among all consensus supporting it)
                    newV[v] |= np_reads_support[name]
                else:
                    # otherwise, flag vertex to remove
                    to_remove.add(v)

    # flag all vertices supported by few reads to remove
    for v, w in newV.items():
        if len(w) < minw:
            to_remove.add(v)
    print(
        f"Vertices to remove: {len(to_remove)} out of {nn} novel vertices",
        file=sys.stderr,
    )

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
            line = line.strip("\n").split("\t")
            name = line[1]
            if name in np_reads_support:
                continue
        print(line, end="")


if __name__ == "__main__":
    main()
