import sys
import argparse
import os
import glob
import re
from collections import Counter
from collections import defaultdict

regex = re.compile("([0-9]+[=XID])")


def parse_cigar(cigar):
    tokens = regex.split(cigar)
    tokens = [(int(x[:-1]), x[-1]) for x in tokens if x != ""]
    return tokens


"""
Ss: {v : (seq, seq_len, is_known)}
Ls: {(v,strand) : {(v, strand} : is_known}
"""


def load_graph(gfa_fn, sample=""):
    print("Parsing GFA...", file=sys.stderr)
    Ss, Ls = {}, {}

    for line in open(gfa_fn):
        if line.startswith("S"):
            _, idx, seq, *_ = line.strip("\n").split("\t")
            idx = int(idx)
            Ss[idx] = [seq, len(seq), False]
        elif line.startswith("L"):
            _, idx1, s1, idx2, s2, *_ = line.strip("\n").split("\t")
            idx1, idx2 = int(idx1), int(idx2)
            # L direction
            k = (idx1, s1 == "+")
            v = (idx2, s2 == "+")
            if k not in Ls:
                Ls[k] = {}
            Ls[k][v] = False
            # add also bidirection
            k = (idx2, s2 != "+")
            v = (idx1, s1 != "+")
            if k not in Ls:
                Ls[k] = {}
            Ls[k][v] = False
        elif line[0] in ["P", "W"]:
            # XXX: assuming P and W are after S and L lines
            fields = line.strip("\n").split("\t")
            is_novel = sample != "" and sample in fields[1]

            if is_novel:
                continue

            path = []
            if line.startswith("P"):
                path = [(int(x[:-1]), x[-1] == "+") for x in fields[2].split(",")]
            elif line.startswith("W"):
                path = list(re.findall(r"[<>]\d+|[^<>]+", fields[6]))
                path = [(int(x[1:]), x[0] == ">") for x in path]

            Ss[path[0][0]][2] = True
            for t1, t2 in zip(path[:-1], path[1:]):
                Ss[t2[0]][2] = True

                assert t1 in Ls
                Ls[t1][t2] = True
                # flag bidirection

                t1 = (t1[0], not t1[1])
                t2 = (t2[0], not t2[1])
                assert t2 in Ls
                Ls[t2][t1] = True

    print("Total vertices:", len(Ss), file=sys.stderr)
    print("Novel vertices:", sum([not Ss[x][2] for x in Ss]), file=sys.stderr)
    return Ss, Ls


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("GFA")
    parser.add_argument("GAF")
    parser.add_argument("-s", type=str, default="")

    parser.add_argument("-t", type=str, required=True)
    parser.add_argument("-n", type=str, required=True)
    parser.add_argument("-c", type=str, required=True)
    parser.add_argument("-l", type=str, required=True)

    args = parser.parse_args()

    vertices, edges = load_graph(args.GFA, args.s)

    Cs = [int(c) for c in args.c.split(",")]

    novel_vertices = {}
    for v, info in vertices.items():
        if not info[2]:
            novel_vertices[v] = []
    novel_edges = {}
    for v1 in edges:
        for v2 in edges[v1]:
            if not edges[v1][v2]:
                novel_edges[(v1, v2)] = 0

    print("tool,n,coverage,clen,kind,v,l,supp,qual")
    for line in open(args.GAF):
        line = line.strip("\n").split("\t")
        path = line[5]
        path = list(re.findall(r"[<>]\d+|[^<>]+", path))
        path = [(int(x[1:]), x[0] == ">") for x in path]

        ps, pe = int(line[7]), int(line[8])

        cigar = line[16].split(":")[-1]
        cigar = parse_cigar(cigar)
        ecigar = []
        for opl, op in cigar:
            ecigar += [op] * opl
        # print(ecigar)

        cp = 0  # where we are along the extended cigar
        plen = 0  # total path length

        # do first vertex since it might be cut
        v = path[0][0]
        plen += vertices[v][1]
        l = vertices[v][1] - ps
        if len(path) == 1:
            l -= vertices[v][1] - pe  # - 1 + 1, pe is not included
        local_ecigar = []
        vbases = 0
        while vbases < l:
            local_ecigar.append(ecigar[cp])
            if ecigar[cp] != "I":
                vbases += 1
            cp += 1

        if vertices[v][2] == False:
            # novel vertez
            op_counts = Counter(local_ecigar)
            novel_vertices[v].append(
                (op_counts["="], op_counts["X"] + op_counts["I"] + op_counts["D"])
            )

        # do internal vertices
        for v, _ in path[1:-1]:
            l = vertices[v][1]
            plen += l
            local_ecigar = []
            vbases = 0
            while vbases < l:
                local_ecigar.append(ecigar[cp])
                if ecigar[cp] != "I":
                    vbases += 1
                cp += 1

            if vertices[v][2] == False:
                op_counts = Counter(local_ecigar)
                novel_vertices[v].append(
                    (op_counts["="], op_counts["X"] + op_counts["I"] + op_counts["D"])
                )

        # do last vertex if path has more then 1 vertex. Last vertex might be cut
        if len(path) > 1:
            v = path[-1][0]
            l = pe - plen

            local_ecigar = []
            vbases = 0
            while vbases < l:
                local_ecigar.append(ecigar[cp])
                if ecigar[cp] != "I":
                    vbases += 1
                cp += 1

            if vertices[v][2] == False:
                op_counts = Counter(local_ecigar)
                novel_vertices[v].append(
                    (op_counts["="], op_counts["X"] + op_counts["I"] + op_counts["D"])
                )

        # do edges separately
        for (v1, strand1), (v2, strand2) in zip(path[:-1], path[1:]):
            if edges[(v1, strand1)][(v2, strand2)] == False:
                novel_edges[((v1, strand1), (v2, strand2))] += 1

    for v, info in novel_vertices.items():
        q = -1
        supp = len(info)
        if supp > 0:
            q = sum([x / (x + y) for x, y in info]) / len(info)
        for c in Cs:
            print(
                args.t,
                args.n,
                c,
                args.l,
                "vertex",
                v,
                vertices[v][1],
                supp,
                q,
                sep=",",
            )

    already_printed = set()
    for ((v1, s1), (v2, s2)), supp in novel_edges.items():
        if ((v1, s1), (v2, s2)) in already_printed:
            continue

        real_supp = supp + novel_edges[((v2, not s2), (v1, not s1))]
        already_printed.add(((v1, s1), (v2, s2)))
        already_printed.add(((v2, not s2), (v1, not s1)))

        for c in Cs:
            print(
                args.t,
                args.n,
                c,
                args.l,
                "edge",
                f"{v1}{'-+'[s1]}>{v2}{'-+'[s2]}",
                ".",
                real_supp,
                "1",
                sep=",",
            )


if __name__ == "__main__":
    main()
