import sys
import argparse
import re

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns

sns.set(style="whitegrid")


def main():
    GAFs = sys.argv[1:-1]
    png_fn = sys.argv[-1]
    NM = []
    NAL = []
    COV = []

    graphs = {
        0: "augmented",
        1: "full",
        2: "1out",
        # 3: "vg",
        # 4: "minigraph-cactus",
    }
    for i, gaf in enumerate(GAFs):
        graph = graphs[i]

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
                    nonoverlapping.append([s, e])
            c = sum([e - s for s, e in nonoverlapping])
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

    # STRAT = [0.25, 0.5, 0.75, 0.9, 0.99, 0.999, 0.9999, 1]
    STRAT = [0.999, 0.9999, 0.99999, 1]
    stratcov = []
    nreads = {}
    for graph in graphs.values():
        print(f"=== {graph} ===")
        #     print(len(NM[NM["Graph"] == graph]))
        full = len(COV[(COV["Graph"] == graph) & (COV["Coverage"] <= 1)])
        nreads[graph] = full
        for n in STRAT:
            partial = len(COV[(COV["Graph"] == graph) & (COV["Coverage"] <= n)])
            stratcov.append(
                [
                    graph,
                    n,
                    partial,
                    full,
                    round(partial / full * 100 if full > 0 else 0, 2),
                ]
            )

    print(nreads)
    stratcov = pd.DataFrame(
        stratcov, columns=["Graph", "n", "Partial", "Full", "Ratio"]
    )

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

    # sns.histplot(
    #     COV,
    #     x="Coverage",
    #     hue="Graph",
    #     element="step",  # "bar"
    #     bins=100,
    #     cumulative=True,  # False
    #     legend=False,
    #     binrange=[0, 1],
    #     ax=ax2,
    # )

    # sns.violinplot(COV[COV["Coverage"] > 0.9998], x="Graph", y="Coverage", ax=ax2)

    cmap = sns.color_palette()
    legend = []
    for i, n in enumerate(reversed(STRAT)):
        print(cmap[i])
        sns.barplot(data=stratcov[stratcov["n"] == n], x="Graph", y="Ratio", ax=ax2)
        legend.append(Patch(facecolor=cmap[i], edgecolor="w", label=n))
        if i != 0:
            ax2.bar_label(ax2.containers[i])
        else:
            ax2.bar_label(ax2.containers[0], labels=nreads.values())

    ax2.legend(handles=legend)
    ax2.set_ylabel("")
    ax2.set_title("Alignment coverage per read (%)")

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
        binrange=[0, 35],
        ax=ax3,
    )
    ax3.set_ylabel("")
    ax3.set_title("NM (per best/longest alignment)")

    plt.tight_layout()
    plt.show()
    # plt.savefig(png_fn)


if __name__ == "__main__":
    main()
