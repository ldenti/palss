import sys


def main():
    gfa_fn = sys.argv[1]

    print("Parsing paths...", file=sys.stderr)
    V = {}
    for line in open(gfa_fn):
        path = []
        if line.startswith("P"):
            line = line.strip("\n").split("\t")
            path = [int(x[:-1]) for x in line[2].split(",")]
        elif line.startswith("W"):
            line = line.strip("\n").split("\t")
            path = [int(x) for x in line[6][1:].split(">")]
        else:
            continue
        for v in path:
            V[v] = ""

    print("Parsing segments...", file=sys.stderr)
    nV = {}
    for line in open(gfa_fn):
        if line.startswith("S"):
            _, v, seq = line.strip("\n").split("\t")
            v = int(v)
            if v in V:
                V[v] = seq
            else:
                nV[v] = seq

    print("Parsing edges...", file=sys.stderr)
    inE = {}
    outE = {}
    for line in open(gfa_fn):
        if line.startswith("L"):
            _, v1, _, v2, *_ = line.split("\t")
            v1, v2 = int(v1), int(v2)
            outE[v1] = outE[v1] + [v2] if v1 in outE else [v2]
            inE[v2] = inE[v2] + [v1] if v2 in inE else [v1]

    print("Selecting segments to remove...", file=sys.stderr)
    toremove = set()
    for v, oe in outE.items():
        for v1 in oe:
            if v1 not in nV:
                continue
            seq1 = nV[v1]
            for v2 in oe:
                if v1 == v2:
                    continue
                seq2 = V[v2] if v2 in V else nV[v2]
                if seq2 == seq1:
                    toremove.add(v1)
    for v, ie in inE.items():
        for v1 in ie:
            if v1 not in nV:
                continue
            seq1 = nV[v1]
            for v2 in ie:
                if v1 == v2:
                    continue
                seq2 = V[v2] if v2 in V else nV[v2]
                if seq2 == seq1:
                    toremove.add(v1)

    print(
        f"Need to remove {len(toremove)} out of {len(nV)} novel vertices",
        file=sys.stderr,
    )

    print(481183, 481183 in toremove, file=sys.stderr)

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
        print(line, end="")


if __name__ == "__main__":
    main()
