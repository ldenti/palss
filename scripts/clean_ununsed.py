import sys
import re


def main():
    gfa_fn = sys.argv[1]
    gaf_fn = sys.argv[2]  # selected consensus, from palss. For support
    re_gaf_fn = sys.argv[3]  # selected consensus, realigned
    minw = int(sys.argv[4])

    print("Parsing GAF...", file=sys.stderr)
    np_reads_support = {}
    for line in open(gaf_fn):
        fields = line.strip("\n").split("\t")
        name = fields[0]
        np_reads_support[name] = set(fields[15].split(":")[-1].split("|"))

    print("Parsing GAF...", file=sys.stderr)
    used_vertices = set()
    for line in open(re_gaf_fn):
        fields = line.strip("\n").split("\t")
        name = fields[0]
        path = fields[5]
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
