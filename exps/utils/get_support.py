import sys
import os
import glob
import re


def get_novel_vertices(gfa_fn, sample=""):
    otherV = set()
    sampleV = set()
    for line in open(gfa_fn):
        fields = line.strip("\n").split("\t")

        path = []
        if line.startswith("P"):
            path = [int(x[:-1]) for x in fields[2].split(",")]
        elif line.startswith("W"):
            path = [int(x) for x in re.split("[<>]", fields[6][1:])]
        else:
            continue

        is_sample = sample != "" and sample in fields[1]
        for v in path:
            if is_sample:
                sampleV.add(v)
            else:
                otherV.add(v)

    if sample == "":
        assert len(sampleV) == 0

    nV = {}
    for line in open(gfa_fn):
        if line.startswith("S"):
            _, v, seq, *_ = line.strip("\n").split("\t")
            v = int(v)

            if sample == "":
                if v not in otherV:
                    nV[v] = [0, len(seq)]
            else:
                if v in sampleV and v not in otherV:
                    nV[v] = [0, len(seq)]
    return nV


def main():
    WD = sys.argv[1]
    sample = sys.argv[2]

    print("graph,augmentation,n,w,d,v,supp,l")

    # PALSS pangenomes
    for gfa_fn in glob.glob(
        os.path.join(WD, "n*", "palss-*", "pangenome-augmented-*.*.gfa")
    ):
        print(gfa_fn, file=sys.stderr)
        n = int(gfa_fn.split("/")[-3][1:])
        run = "oneout" if "oneout" in gfa_fn else "full"
        fn = gfa_fn.split("/")[-1]
        w = int(fn.split(".")[-2][1:])
        d = float(fn.split(".")[-4][1:] + "." + fn.split(".")[-3])
        augment = fn.split(".")[-5].split("-")[-1]

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

    # Other pangenomes
    for gfa_fn in glob.glob(os.path.join(WD, "n*", "pangenome-*.gfa")):
        print(gfa_fn, file=sys.stderr)
        n = int(gfa_fn.split("/")[-2][1:])
        run = gfa_fn.split("/")[-1][:-4].split("-")[-1]
        if run not in ["full", "mgcactus"]:
            continue

        novel_vertices = get_novel_vertices(gfa_fn, sample)

        fn = f"{run}.gaf"
        if run != "mgcactus":
            fn = f"original-{run}.gaf"
        gaf_fn = os.path.join(WD, f"n{n}", "graphaligner", fn)
        for line in open(gaf_fn):
            line = line.strip("\n").split("\t")
            path = line[5]
            path = [int(x) for x in re.split("[<>]", path[1:])]
            for v in path:
                if v in novel_vertices:
                    novel_vertices[v][0] += 1
        for v, (supp, length) in novel_vertices.items():
            print(run, run, n, -1, -1, v, supp, length, sep=",", flush=False)
        sys.stdout.flush()


if __name__ == "__main__":
    main()
