import sys
import os
import glob
import re
from collections import Counter

regex = re.compile("([0-9]+[=XID])")


def parse_cigar(cigar):
    tokens = regex.split(cigar)
    tokens = [(int(x[:-1]), x[-1]) for x in tokens if x != ""]
    return tokens


def get_novel_vertices(gfa_fn, sample=""):
    allV = {}
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
                if fields[1].startswith("new.v"):
                    continue
                otherV.add(v)

    if sample == "":
        assert len(sampleV) == 0

    nV = {}
    for line in open(gfa_fn):
        if line.startswith("S"):
            _, v, seq, *_ = line.strip("\n").split("\t")
            v = int(v)

            allV[v] = len(seq)

            if sample == "":
                if v not in otherV:
                    nV[v] = []
            else:
                if v in sampleV and v not in otherV:
                    nV[v] = []
    return allV, nV


def get_novel_edges(gfa_fn, nV, sample=""):
    otherE = set()
    sampleE = set()
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
        for v1, v2 in zip(path[:-1], path[1:]):
            e = (min(v1, v2), max(v1, v2))
            if is_sample:
                sampleE.add(e)
            else:
                otherE.add(e)

    if sample == "":
        assert len(sampleE) == 0

    nE = {}
    for line in open(gfa_fn):
        if line.startswith("L"):
            _, v1, _, v2, _, *_ = line.strip("\n").split("\t")
            v1, v2 = int(v1), int(v2)
            e = (min(v1, v2), max(v1, v2))
            if sample == "":
                if e not in otherE:  # and v1 not in nV and v2 not in nV:
                    nE[e] = 0
            else:
                if (
                    e in sampleE and e not in otherE
                ):  # and v1 not in nV and v2 not in nV:
                    nE[e] = 0
    return nE


def main():
    gfa_fn = sys.argv[1]
    gaf_fn = sys.argv[2]
    sample = sys.argv[3]

    print("graph,augmentation,n,w,d,iden,kind,v,l,supp,qual")

    n = -1
    run = ""
    augmentation = ""
    i = -1
    w = -1
    d = -1

    fold = gfa_fn.split("/")[-2]
    vertices, novel_vertices = {}, {}
    if "palss" in fold:
        n = int(gfa_fn.split("/")[-3][1:])
        run = "oneout" if "oneout" in gfa_fn else "full"
        augmentation = "palss"
        fn = gfa_fn.split("/")[-1]
        i = float(fn.split(".")[-3][2:] + "." + fn.split(".")[-2])
        w = int(fn.split(".")[-4][1:])
        d = float(fn.split(".")[-6][1:] + "." + fn.split(".")[-5])

        vertices, novel_vertices = get_novel_vertices(gfa_fn)
        # novel_edges = get_novel_edges(gfa_fn, novel_vertices)
    else:
        n = int(gfa_fn.split("/")[-2][1:])
        run = "oneout"
        augmentation = gfa_fn.split("/")[-1][:-4].split("-")[-1]
        assert augmentation in ["full", "oneout", "mgcactus"]

        vertices, novel_vertices = get_novel_vertices(gfa_fn, sample)
        # novel_edges = get_novel_edges(gfa_fn, sample)

    for line in open(gaf_fn):
        line = line.strip("\n").split("\t")
        path = line[5]
        path = [int(x) for x in re.split("[<>]", path[1:])]
        ps, pe = int(line[7]), int(line[8])

        cigar = line[16].split(":")[-1]
        cigar = parse_cigar(cigar)
        ecigar = []
        for opl, op in cigar:
            ecigar += [op] * opl
        # print(ecigar)

        cp = 0  # where we are along the extended cigar
        plen = 0  # total path length

        v = path[0]
        plen += vertices[v]
        l = vertices[v] - ps  # first vertex might be cut
        if len(path) == 1:
            l -= vertices[v] - pe  # - 1 + 1, pe is not included
        local_ecigar = []
        vbases = 0
        while vbases < l:
            local_ecigar.append(ecigar[cp])
            if ecigar[cp] != "I":
                vbases += 1
            cp += 1

        if v in novel_vertices:
            # novel_vertices[v] += 1
            op_counts = Counter(local_ecigar)
            novel_vertices[v].append(
                (op_counts["="], op_counts["X"] + op_counts["I"] + op_counts["D"])
            )

        for v in path[1:-1]:
            l = vertices[v]
            plen += vertices[v]
            local_ecigar = []
            vbases = 0
            while vbases < l:
                local_ecigar.append(ecigar[cp])
                if ecigar[cp] != "I":
                    vbases += 1
                cp += 1

            if v in novel_vertices:
                op_counts = Counter(local_ecigar)
                novel_vertices[v].append(
                    (op_counts["="], op_counts["X"] + op_counts["I"] + op_counts["D"])
                )

        if len(path) > 1:
            v = path[-1]
            l = pe - plen  # last vertex might be cut

            local_ecigar = []
            vbases = 0
            while vbases < l:
                local_ecigar.append(ecigar[cp])
                if ecigar[cp] != "I":
                    vbases += 1
                cp += 1

            if v in novel_vertices:
                op_counts = Counter(local_ecigar)
                novel_vertices[v].append(
                    (op_counts["="], op_counts["X"] + op_counts["I"] + op_counts["D"])
                )

    for v, info in novel_vertices.items():
        q = -1
        supp = len(info)
        if supp > 0:
            q = sum([x / (x + y) for x, y in info]) / len(info)
        print(
            run,
            augmentation,
            n,
            w,
            d,
            i,
            "vertex",
            v,
            vertices[v],
            supp,
            q,
            sep=",",
        )
    # for (v1, v2), supp in novel_edges.items():
    #     print(
    #         run,
    #         f"palss{pv}",
    #         n,
    #         w,
    #         d,
    #         i,
    #         "edge",
    #         f"{v1}.{v2}",
    #         supp,
    #         -1,
    #         sep=",",
    #     )


if __name__ == "__main__":
    main()
