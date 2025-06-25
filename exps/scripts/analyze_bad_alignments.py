import sys
import argparse
import re

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns

sns.set(style="whitegrid")


def main():
    GAFs = sys.argv[1:]
    # png_fn = sys.argv[-1]

    NM = {}
    COV = {}

    graphs = {
        0: "augmented",
        1: "full",
        2: "1out",
        # 3: "vg",
        # 4: "minigraph-cactus",
    }
    for i, gaf in enumerate(GAFs):
        graph = i  # graphs[i]
        nms = {}
        coverages = {}

        NM[graph] = {}
        COV[graph] = {}

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
            NM[graph][qidx] = [c, nm]

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
            COV[graph][qidx] = float(c) / l

    QIDXS = set(NM[0].keys()) & set(NM[1].keys()) & set(NM[2].keys())
    total = 0
    partial = 0

    better_nm = []
    for qidx in QIDXS:
        total += 1
        c0, nm0 = NM[0][qidx]
        c1, nm1 = NM[1][qidx]
        c2, nm2 = NM[2][qidx]

        # if c0 != 1 or c1 != 1 or c2 != 1:
        #     partial += 1
        #     continue

        fc0, fc1, fc2 = COV[0][qidx], COV[1][qidx], COV[2][qidx]
        if nm1 < nm2:
            if nm0 < nm2:
                # print("<")
                # better_nm.append([qidx, 0, nm0])
                # better_nm.append([qidx, 1, nm1])
                # better_nm.append([qidx, 2, nm2])
                better_nm.append([qidx, nm2 - nm0, nm2 - nm1])
                if nm0 < nm1:
                    print(qidx, nm0, nm1, nm2)

            elif nm0 == nm2:
                pass  # print("=")
            else:
                pass  # print(">")
        elif nm1 == nm2:
            pass  # print("==")
        else:
            pass  # print(">>")

    print(partial, total, partial / total)

    df = pd.DataFrame(better_nm, columns=["qidx", "2-0", "2-1"])
    sns.scatterplot(df, x="2-0", y="2-1", alpha=0.5)

    xlim = plt.xlim()
    ylim = plt.ylim()

    max_lim = max(xlim[1], ylim[1])

    plt.xlim(xlim[0], max_lim)
    plt.ylim(ylim[0], max_lim)

    plt.xlabel("Δ(aug,1out)")
    plt.ylabel("Δ(full,1out)")
    plt.show()


if __name__ == "__main__":
    main()
