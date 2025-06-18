import sys
import argparse
import re

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def main():
    GAFs = sys.argv[1:-1]
    png_fn = sys.argv[-1]
    NM = []
    NAL = []
    COV = []
    for i, gaf in enumerate(GAFs):
        graph = "augmented"
        if i == 0:
            graph = "augmented"
        elif i == 1:
            graph = "full"
        elif i == 2:
            graph = "1out"
        nalignments = {}
        coverages = {}
        nms = {}
        print(f"Iterating over {gaf}...")
        for line in open(gaf):
            line = line.strip("\n").split("\t")
            qidx = line[0]
            ql = int(line[1])
            qs = int(line[2])  # closed
            qe = int(line[3])  # open

            if qidx not in coverages:
                coverages[qidx] = []
            coverages[qidx].append((qs, qe, ql))  # we can avoid repeating ql

            nalignments[qidx] = nalignments[qidx] + 1 if qidx in nalignments else 1

            nm = line[12]
            nm = nm.split(":")
            assert nm[0] == "NM"
            nm = int(nm[2])
            c = (qe - qs) / ql
            if qidx not in nms:
                nms[qidx] = (c, nm)
            else:
                old_c, old_nm = nms[qidx]
                if old_c == c:
                    nms[qidx] = (c, min(nm, old_nm))
                elif old_c < c:
                    nms[qidx] = (c, nm)
                else:
                    pass

        for qidx, (c, nm) in nms.items():
            NM.append([graph, nm])

        for qidx, fragments in coverages.items():
            fragments = list(set(fragments))
            l = fragments[0][-1]
            fragments.sort(key=lambda x: x[0])
            nonoverlapping = [[fragments[0][0], fragments[0][1]]]
            for s, e, _ in fragments[1:]:
                if s < nonoverlapping[-1][1]:
                    if e > nonoverlapping[-1][1]:
                        nonoverlapping[-1][1] = e
                else:
                    nonoverlapping.append([s,e])
            c = sum([e-s for s,e in nonoverlapping])
            COV.append([graph, float(c) / l])

        print(f"=== {graph}: ", sum(nalignments.values()))
        n = 5
        bp = {}
        for i in range(n):
            bp[str(i)] = 0
        bp[f"{n}+"] = 0
        for v in nalignments.values():
            # if v <= 10:
            k = str(v)
            if v >= n:
                k = f"{n}+"
            bp[k] += 1
        for k, v in bp.items():
            NAL.append([graph, k, v])

    NAL = pd.DataFrame(NAL, columns=["Graph", "#Al", "Count"])
    COV = pd.DataFrame(COV, columns=["Graph", "Coverage"])
    NM = pd.DataFrame(NM, columns=["Graph", "NM"])

    for graph in ["augmented", "full", "1out"]:
        #     print(f"=== {graph} ===")
        #     print(len(NM[NM["Graph"] == graph]))
        print(len(COV[(COV["Graph"] == graph) & (COV["Coverage"] == 1)]))
        #     print(len(NAL[NAL["Graph"] == graph]))

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(11, 5))

    sns.barplot(
        NAL,
        x="#Al",
        y="Count",
        hue="Graph",
        fill=True,
        alpha=0.5,
        linewidth=1,
        legend=False,
        gap=0.05,
        ax=ax1,
    )
    sns.barplot(
        NAL,
        x="#Al",
        y="Count",
        hue="Graph",
        fill=False,
        linewidth=1,
        legend=False,
        gap=0.05,
        ax=ax1,
    )
    ax1.set_title("Number of alignments per read")
    ax1.set_xlabel("#Alignments")

    y = sns.histplot(
        COV,
        x="Coverage",
        hue="Graph",
        element="step", # "bar"
        bins=100,
        cumulative=True, # False
        legend=False,
        binrange=[0, 0.99],
        ax=ax2,
    )
    ax2.set_ylabel("")
    ax2.set_title("Alignment coverage per read (no 1.0)")

    # print(".", sum(ax2.containers[0].datavalues))
    # print(sum(ax2.containers[1].datavalues))
    # print(sum(ax2.containers[2].datavalues))

    sns.histplot(
        NM,
        x="NM",
        hue="Graph",
        discrete=True,
        element="step",
        legend=True,
        ax=ax3,
        binrange=[0, 150],
    )
    ax3.set_ylabel("")
    ax3.set_title("NM (per best/longest alignment)")

    plt.tight_layout()
    # plt.show()
    plt.savefig(png_fn)


if __name__ == "__main__":
    main()
