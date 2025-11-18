import sys
import re
from collections import Counter


def main():
    gfa_fn = sys.argv[1]
    samples_fn = sys.argv[2]

    samples = set()
    for line in open(samples_fn):
        samples.add(line.strip("\n"))

    to_keep = {}
    # 0: to remove since in sample
    # 1: to keep since not in sample
    for line in open(gfa_fn):
        if line.startswith("S"):
            _, v, _ = line.split("\t")
            v = int(v)
            if v not in to_keep:
                # all vertices are by default "not-to-keep"
                to_keep[v] = 0
        elif line.startswith("P") or line.startswith("W"):
            name, path = "", []
            if line.startswith("P"):
                _, idx, path, _ = line.split("\t")
                name = idx.split("#")[0]
                path = [int(x[:-1]) for x in path.split(",")]
            else:
                _, name, haplotype, chrom, _, _, path = line.split("\t")
                path = [int(x) for x in re.split("[<>]", path[1:])]
            if name in samples:
                # we need to keep these vertices
                for v in path:
                    to_keep[v] = 1

    c = Counter(to_keep.values())
    print(c, file=sys.stderr)

    for line in open(gfa_fn):
        if line.startswith("S"):
            _, v, _ = line.split("\t")
            v = int(v)
            if to_keep[v] == 0:
                continue
        elif line.startswith("L"):
            _, v1, _, v2, _, _ = line.split("\t")
            v1, v2 = int(v1), int(v2)
            if to_keep[v1] == 0 or to_keep[v2] == 0:
                continue
        elif line.startswith("P"):
            name = line.split("\t")[1].split("#")[0]
            if name not in samples:
                continue
        elif line.startswith("W"):
            name = line.split("\t")[1]
            if name not in samples:
                continue
        print(line, end="")


if __name__ == "__main__":
    main()
