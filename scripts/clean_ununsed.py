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
    gfa_fn = sys.argv[1]
    gaf_fn = sys.argv[2]  # selected consensus, from palss. For support
    re_gaf_fn = sys.argv[3]  # selected consensus, realigned (pass1)
    re2_gaf_fn = sys.argv[4]  # selected consensus, realigned (pass2)
    minw = int(sys.argv[5])

    print("Getting support from palss GAF file...", file=sys.stderr)
    np_reads_support = {}
    for line in open(gaf_fn):
        fields = line.strip("\n").split("\t")
        name = fields[0]
        np_reads_support[name] = set(fields[15].split(":")[-1].split("|"))

    print("Getting best alignments from pass1 GAF file...", file=sys.stderr)
    pass1_alns = parse_gaf(re_gaf_fn)
    print("Getting best alignments from pass1 GAF file...", file=sys.stderr)
    pass2_alns = parse_gaf(re2_gaf_fn)

    used_vertices = set()
    bad_alignments = set()
    for name, (c, nm, cigar, path) in pass2_alns.items():
        if name in pass1_alns:
            p1_c, p1_nm, p1_cigar, p1_path = pass1_alns[name]
            if p1_nm < nm:
                # this consensus was aligned better to pass1 graph, maybe due to vertex splitting (see comment in augment.sh)
                # XXX: it seems better to not do this
                bad_alignments.add(name)
                # continue
        used_vertices |= set([int(x) for x in re.split("[<>]", path[1:])])

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
                if name in bad_alignments:
                    # this consensus has been aligned bad to the new graph (pass2). So we keep all vertices anyway
                    # XXX: it seems better to not do this
                    # used_vertices |= set(path)
                    pass
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

    print("Parsing segments...", file=sys.stderr)
    toremove = set()
    newV = {}
    nn = 0
    for line in open(gfa_fn):
        if line.startswith("S"):
            _, v, seq = line.strip("\n").split("\t")
            v = int(v)
            if v in knownV:
                pass
            else:
                nn += 1
                if v not in used_vertices:
                    toremove.add(v)
                else:
                    newV[v] = set()
    print(
        f"Need to remove {len(toremove)} out of {nn} novel vertices",
        file=sys.stderr,
    )

    print(
        f"Selecting segments to remove (due to low support, threshold: {minw})...",
        file=sys.stderr,
    )
    for name, path in new_paths.items():
        for v in path:
            if v in newV:
                newV[v] |= np_reads_support[name]
    for v, w in newV.items():
        if len(w) < minw:
            toremove.add(v)
    print(
        f"Need to remove {len(toremove)} out of {nn} novel vertices",
        file=sys.stderr,
    )

    print("Outputting...", file=sys.stderr)
    for line in open(gfa_fn):
        if line.startswith("S"):
            _, v, seq = line.strip("\n").split("\t")
            v = int(v)
            if v in toremove:
                continue
        elif line.startswith("L"):
            _, v1, _, v2, *_ = line.split("\t")
            v1, v2 = int(v1), int(v2)
            if v1 in toremove or v2 in toremove:
                continue
        elif line.startswith("P"):
            line = line.strip("\n").split("\t")
            name = line[1]
            if name in np_reads_support:
                continue
        print(line, end="")


if __name__ == "__main__":
    main()
