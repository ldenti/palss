import sys
import re

complement = {"A": "T", "T": "A", "C": "G", "G": "C"}


def rc(sequence):
    return "".join(complement[base] for base in reversed(sequence))


def main():
    gfa_fn = sys.argv[1]
    gaf_fn = sys.argv[2]
    minw = int(sys.argv[3])

    print("Parsing GAF...", file=sys.stderr)
    np_supports = {}
    for line in open(gaf_fn):
        fields = line.strip("\n").split("\t")
        name = fields[0]
        support = int(fields[15].split(":")[-1])
        np_supports[name] = support

    print("Parsing paths...", file=sys.stderr)
    new_paths = {}
    V = {}
    for line in open(gfa_fn):
        path = []
        skip = False
        if line.startswith("P"):
            line = line.strip("\n").split("\t")
            name = line[1]
            path = [int(x[:-1]) for x in line[2].split(",")]
            if name in np_supports:
                new_paths[name] = path
                skip = True
        elif line.startswith("W"):
            line = line.strip("\n").split("\t")
            path = [int(x) for x in re.split("[<>]", line[6][1:])]
        else:
            continue
        if skip:
            # path is a new path, do not add its vertices to V
            continue
        for v in path:
            V[v] = ""

    print("Parsing segments...", file=sys.stderr)
    nV = {}
    # nV_info = {}
    for line in open(gfa_fn):
        if line.startswith("S"):
            _, v, seq = line.strip("\n").split("\t")
            v = int(v)
            if v in V:
                pass  # V[v] = seq
            else:
                nV[v] = 0  # seq
                # nV_info[v] = []

    # print("Parsing links...", file=sys.stderr)
    # inE = {}
    # outE = {}
    # for line in open(gfa_fn):
    #     if line.startswith("L"):
    #         _, v1, _, v2, *_ = line.split("\t")
    #         v1, v2 = int(v1), int(v2)
    #         outE[v1] = outE[v1] + [v2] if v1 in outE else [v2]
    #         inE[v2] = inE[v2] + [v1] if v2 in inE else [v1]
    #
    # print("Selecting segments to remove...", file=sys.stderr)
    toremove = set()
    # for _, oe in outE.items():
    #     for v1 in oe:
    #         # for each outgoing vertex that is novel
    #         if v1 not in nV:
    #             continue
    #         seq1 = nV[v1]
    #         # check if we already have it
    #         for v2 in oe:
    #             if v1 == v2:
    #                 continue
    #             seq2 = V[v2] if v2 in V else nV[v2]
    #             if seq2 == seq1 or rc(seq2) == seq1:
    #                 toremove.add(v1)
    # for v, ie in inE.items():
    #     for v1 in ie:
    #         if v1 not in nV:
    #             continue
    #         seq1 = nV[v1]
    #         for v2 in ie:
    #             if v1 == v2:
    #                 continue
    #             seq2 = V[v2] if v2 in V else nV[v2]
    #             if seq2 == seq1:
    #                 toremove.add(v1)

    print("Selecting segments to remove (support)...", file=sys.stderr)
    for name, path in new_paths.items():
        for v in path:
            if v in nV:
                # if not isinstance(nV[v], int):
                #     nV[v] = 0
                nV[v] += np_supports[name]
                # nV_info[v].append(name)
    for v, w in nV.items():
        if w < minw:
            toremove.add(v)

    print(
        f"Need to remove {len(toremove)} out of {len(nV)} novel vertices",
        file=sys.stderr,
    )

    print("Outputting...", file=sys.stderr)
    for line in open(gfa_fn):
        if line.startswith("S"):
            _, v, seq = line.strip("\n").split("\t")
            v = int(v)
            if v in toremove:
                continue
            # if v in nV:
            #     line = line[:-1] + "\tNV:i:1\n"
        elif line.startswith("L"):
            _, v1, _, v2, *_ = line.split("\t")
            v1, v2 = int(v1), int(v2)
            if v1 in toremove or v2 in toremove:
                continue
        elif line.startswith("P"):
            line = line.strip("\n").split("\t")
            name = line[1]
            if name in np_supports:
                continue
        print(line, end="")


if __name__ == "__main__":
    main()
