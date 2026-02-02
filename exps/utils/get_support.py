import sys
import os
import glob
import re


def get_novel_vertices(gfa_fn):
    V = set()
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
    for line in open(gfa_fn):
        if line.startswith("S"):
            _, v, seq, *_ = line.strip("\n").split("\t")
            v = int(v)
            if v in V:
                continue
            nV[v] = [0, len(seq)]
    return nV


def main():
    WD = sys.argv[1]

    print("graph,augmentation,n,w,d,v,supp,l")
    data = []
    for gfa_fn in glob.glob(
        os.path.join(WD, "n*", "palss-*", "pangenome-augmented-*.*.gfa")
    ):
        n = int(gfa_fn.split("/")[-3][1:])
        run = "oneout" if "oneout" in gfa_fn else "full"
        fn = gfa_fn.split("/")[-1]
        w = int(fn.split(".")[-2][1:])
        d = float(fn.split(".")[-4][1:] + "." + fn.split(".")[-3])
        augment = fn.split(".")[-5].split("-")[-1]

        if augment not in ["simple", "medium"]:
            continue

        novel_vertices = get_novel_vertices(gfa_fn)

        gaf_fn = os.path.join(
            WD, f"n{n}", "graphaligner", f"palss-{run}.{augment}.d{d}.w{w}.gaf"
        )

        for line in open(gaf_fn):
            line = line.strip("\n").split("\t")
            path = line[5]
            path = [int(x) for x in re.split("[<>]", path[1:])]
            for v in path:
                if v in novel_vertices:
                    novel_vertices[v][0] += 1
        for v, (supp, length) in novel_vertices.items():
            print(run, augment, n, w, d, v, supp, length, sep=",", flush=False)
        sys.stdout.flush()


if __name__ == "__main__":
    main()
