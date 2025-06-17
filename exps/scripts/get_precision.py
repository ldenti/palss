import sys
import argparse
import re

import matplotlib.pyplot as plt
import seaborn as sns


def main():
    gfa_fn = sys.argv[1]
    gaf_fn = sys.argv[2]
    png_fn = sys.argv[3]

    Vseqs = {}
    V = set()
    print("Parsing paths...", file=sys.stderr)
    for line in open(gfa_fn):
        path = []
        if line.startswith("P"):
            line = line.strip("\n").split("\t")
            path = [int(x[:-1]) for x in line[2].split(",")]
        elif line.startswith("W"):
            line = line.strip("\n").split("\t")
            path = [int(x) for x in re.split("[<>]", line[6][1:])]
        else:
            continue
        for v in path:
            V.add(v)

    nV = {}
    # Vseqs = {}
    nV_len = {}
    print("Reiterating over graph...", file=sys.stderr)
    for line in open(gfa_fn):
        if line.startswith("S"):
            _, v, seq, *_ = line.strip("\n").split("\t")
            v = int(v)
            # Vseqs[v] = seq
            if v in V:
                continue
            nV_len[v] = len(seq)
            nV[v] = 0

    print(len(nV), file=sys.stderr)

    print("Iterating over GAF...", file=sys.stderr)
    for line in open(gaf_fn):
        line = line.strip("\n").split("\t")

        path = line[5]
        strand = path[0]
        path = [int(x) for x in re.split("[<>]", path[1:])]

        for v in path:
            if v in nV:
                nV[v] = nV[v] + 1
    zeros = [x for x, k in nV.items() if k == 0]
    ones = [x for x, k in nV.items() if k == 1]
    for zero in zeros:
        print(zero)
        # print(f">{zero}")
        # print(Vseqs[zero])
    print(0, len(zeros), len(nV), len(zeros) / len(nV), file=sys.stderr)
    print(1, len(ones), len(nV), len(ones) / len(nV), file=sys.stderr)

    ax = sns.histplot(
        nV.values(),
        # hue="Graph",
        discrete=True,
        element="step",
        legend=None,
        binrange=[0,30],
    )
    plt.title("Support for novel vertices")
    plt.xlabel("Support")
    plt.ylabel("#Vertices")
    plt.tight_layout()
    # plt.show()
    plt.savefig(png_fn)


if __name__ == "__main__":
    main()
