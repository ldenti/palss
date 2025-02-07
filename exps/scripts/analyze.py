import argparse
import re
import pickle


def main(args):
    V = set()
    E = set()
    print("Parsing paths...")
    for line in open(args.GFA):
        path = []
        if line.startswith("P"):
            line = line.strip("\n").split("\t")
            path = [int(x[:-1]) for x in line[2].split(",")]
        elif line.startswith("W"):
            line = line.strip("\n").split("\t")
            path = [int(x) for x in re.split("[<>]", line[6][1:])]
        else:
            continue
        V.add(path[0])
        for v1, v2 in zip(path[:-1], path[1:]):
            V.add(v2)
            e = f"{v1}:{v2}"
            E.add(e)

    nV = {}
    nV_len = {}
    nE = {}
    print("Reiterating over graph...")
    for line in open(args.GFA):
        if line.startswith("S"):
            _, v, seq, *_ = line.strip("\n").split("\t")
            v = int(v)
            if v in V:
                continue
            nV_len[v] = len(seq)
            nV[v] = 0
        elif line.startswith("L"):
            _, v1, _, v2, _, _ = line.strip("\n").split("\t")
            e = f"{v1}:{v2}"
            if e in E:
                continue
            nE[e] = 0

    print(len(nV), len(nE))

    NM = []
    alignments = {}
    print("Iterating over GAF...")
    for line in open(args.GAF):
        line = line.strip("\n").split("\t")
        qidx = line[0]
        alignments[qidx] = alignments[qidx] + 1 if qidx in alignments else 1

        path = line[5]
        strand = path[0]
        path = [int(x) for x in re.split("[<>]", path[1:])]
        if strand == "<":
            path = path[::-1]

        v1 = path[0]
        if v1 in nV:
            nV[v1] = nV[v1] + 1

        for v1, v2 in zip(path[:-1], path[1:]):
            if v2 in nV:
                nV[v2] = nV[v2] + 1
            e = f"{v1}:{v2}"
            if e in nE:
                nE[e] = nE[e] + 1

        nm = line[13]
        nm = nm.split(":")
        assert nm[0] == "NM"
        nm = int(nm[2])
        NM.append(nm)

    with open(args.out, "wb") as of:
        pickle.dump(alignments, of, pickle.HIGHEST_PROTOCOL)
        pickle.dump(NM, of, pickle.HIGHEST_PROTOCOL)
        pickle.dump(nV, of, pickle.HIGHEST_PROTOCOL)
        pickle.dump(nE, of, pickle.HIGHEST_PROTOCOL)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="",
        description="",
    )
    parser.add_argument("GFA", help="")
    parser.add_argument("GAF", help="")
    parser.add_argument(
        "-o",
        dest="out",
        help="Output pickle file (default: out.pickle)",
        required=True,
        default="out.pickle",
    )
    args = parser.parse_args()
    main(args)
