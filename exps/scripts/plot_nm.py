import sys
import os
import glob

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns

sns.set(style="whitegrid")

GRAPHS = ["1OUT", "FULL", "1OUT-AUG"]


def main():
    WD = sys.argv[1]
    png_fn = sys.argv[2]
    NM = []

    Ns = set()

    for gaf in glob.glob(os.path.join(WD, "*", "alignments*.gaf")):
        print(f"Analyzing {gaf}...")

        n = int(gaf.split("/")[-2])
        Ns.add(n)

        graph = ""
        if "w" in gaf:
            w = int(gaf.split(".")[-3][1:])
            if w != 3:
                continue
            graph = gaf.split("/")[-1].split("-")[1]
            graph = graph.upper() + "-AUG"
        else:
            graph = gaf.split("/")[-1].split("-")[1][:-4]
            graph = graph.upper()
        if graph not in GRAPHS:
            continue
        print(graph)

        nms = {}
        for line in open(gaf):
            line = line.strip("\n").split("\t")
            qidx = line[0]
            ql = int(line[1])
            qs = int(line[2])  # closed
            qe = int(line[3])  # open

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
            NM.append([n, graph, nm])

    Ns = list(Ns)
    Ns.sort()

    NM = pd.DataFrame(NM, columns=["n", "Graph", "NM"])

    fig, axes = plt.subplots(1, len(Ns), sharey=True)  # , figsize=(11, 5))
    for i, n in enumerate(Ns):
        sns.histplot(
            NM[NM["n"] == n],
            x="NM",
            hue="Graph",
            discrete=True,
            element="step",
            legend=True,
            binrange=[0, 35],
            ax=axes[i],
        )
        axes[i].set_ylabel("")
        axes[i].set_title(n)
    plt.suptitle("NM (per best/longest alignment)")
    plt.tight_layout()
    plt.show()
    # plt.savefig(png_fn)


if __name__ == "__main__":
    main()
