import sys
import re


def main():
    gaf_fn = sys.argv[1]
    gfa_fn = sys.argv[2]

    alns = set()
    for line in open(gaf_fn):
        alns.add(line.strip("\n").split("\t")[0])

    knownV = set()
    for line in open(gfa_fn):
        path = []
        skip = False
        if line.startswith("P"):
            line = line.strip("\n").split("\t")
            name = line[1]
            path = line[2]
            if name in alns:
                skip = True
            else:
                path = [int(x[:-1]) for x in line[2].split(",")]
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

    for line in open(gfa_fn):
        path = []
        skip = False
        if line[0] in ["H", "S", "L", "W"]:
            # XXX: assuming that "augmented" paths are P lines
            print(line, end="")
        elif line.startswith("P"):
            fields = line.strip("\n").split("\t")
            name = fields[1]
            path = fields[2]
            if name in alns:
                p = 0
                s = -1
                vertices = [(int(x[:-1]), x[-1]) for x in path.split(",")]
                # print(vertices, file=sys.stderr)
                # print([x[0] in knownV for x in vertices], file=sys.stderr)
                while p < len(vertices) and vertices[p][0] in knownV:
                    p += 1
                while s >= 0 and vertices[s][0] in knownV:
                    s -= 1
                if p <= s:
                    # if we have deletion, we just have reference vertices, so keep the path
                    vertices = vertices[p : s + 1]
                print(
                    "P",
                    name,
                    ",".join([str(x) + s for x, s in vertices]),
                    "*",
                    sep="\t",
                )
            else:
                print(line, end="")


if __name__ == "__main__":
    main()
