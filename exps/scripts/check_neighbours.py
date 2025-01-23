import sys
import pickle


def main():
    gaf_fn = sys.argv[1]
    gfa_fn = sys.argv[2]

    print("Parsing paths...")
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
            if v == 234:
                print(v)
            V[v] = 0

    nV = {}
    print("Reiterating over graph...")
    for line in open(gfa_fn):
        if line.startswith("S"):
            _, v, seq = line.strip("\n").split("\t")
            v = int(v)
            if v in V:
                continue
            nV[v] = [0, seq, [], []]
    print(len(V), len(nV))

    for line in open(gaf_fn):
        line = line.strip("\n").split("\t")
        path = line[5]
        strand = path[0]
        path = [int(x) for x in path[1:].split(strand)]
        if strand == "<":
            path = path[::-1]

        for v in path:
            if v in V:
                V[v] = V[v] + 1
            elif v in nV:
                nV[v][0] = nV[v][0] + 1
            else:
                assert False

    for line in open(gfa_fn):
        if line.startswith("L"):
            _, v1, _, v2, *_ = line.split("\t")
            v1, v2 = int(v1), int(v2)

            if v1 in nV and v2 in nV:
                continue

            if v1 in nV:
                nV[v1][3].append(v2)
            if v2 in nV:
                nV[v2][2].append(v1)

    for v in nV:
        w, seq, incoming, outgoing = nV[v]
        if w == 0:
            print(v, seq)
            for i in incoming:
                assert i in V
                print("<", i, V[i])
            for o in outgoing:
                assert o in V
                print(">", o, V[o])


if __name__ == "__main__":
    main()
