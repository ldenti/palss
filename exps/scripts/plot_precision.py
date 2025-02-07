import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
from matplotlib import pyplot as plt
import pickle
import glob
import os

from matplotlib.patches import Patch

sns.set(style="whitegrid")


def main():
    wd = sys.argv[1]

    NOVEL_1OUT = {}
    NOVEL_FULL = {}
    data = []
    runs = {}
    Gs = set()
    for f in glob.glob(os.path.join(wd, "*", "*.pkl")):
        n, fname = f.split("/")[-2:]
        n = int(n)
        Gs.add(n)

        t = os.path.splitext(fname)[0]
        if t.endswith("w3") or t.endswith("w5"):
            continue

        fb = open(f, "rb")
        AL = pickle.load(fb)
        NM = pickle.load(fb)
        nV = pickle.load(fb)
        nE = pickle.load(fb)
        fb.close()
        if "mgcactus" in t:
            continue
        if "augmented" not in t:
            assert len(nV) == 0 and len(nE) == 0
            continue

        if "1out" in t:
            runs[1] = "1OUT-AUG"
            NOVEL_1OUT[n] = nV
            t = 1  # "1OUT-AUG"
        else:
            runs[2] = "FULL-AUG"
            t = 2  # "FULL-AUG"
            NOVEL_FULL[n] = nV
        xv = len([x for x in nV.values() if x > 0])
        xe = len([x for x in nE.values() if x > 0])

        data.append([n, t, "V", len(nV), xv])
        print([n, t, "V", len(nV), xv])
        data.append([n, t, "E", len(nE), xe])

        print(n, t, "V", len(nV), xv, xv / len(nV) if len(nV) > 0 else 0)
        print(n, t, "E", len(nE), xe, xe / len(nE) if len(nE) > 0 else 0)

    df = pd.DataFrame(data, columns=["G", "Run", "Elem", "Added", "Used"])

    Gs = sorted(list(Gs))
    x = np.arange(len(Gs))  # the label locations
    width = 0.25  # the width of the bars
    multiplier = 0
    colors = [["royalblue", "lightsteelblue"], ["darkorange", "bisque"]]

    fig, (ax1, ax2, ax3) = plt.subplots(
        1, 3, layout="constrained", figsize=(11, 4), dpi=80
    )

    for r in [1, 2]:
        offset = width * multiplier

        d = df[(df["Run"] == r) & (df["Elem"] == "V")].sort_values(by="G")

        rects = ax1.bar(
            x + offset,
            d["Added"],
            width=width,
            label=runs[r],
            align="edge",
            color=colors[r - 1][0],
        )
        ax1.bar(
            x + offset,
            d["Used"],
            width=width,
            label=runs[r],
            align="edge",
            color=colors[r - 1][1],
        )
        ax1.bar_label(rects, labels=round(d["Used"] / d["Added"], 2), padding=3)

        d = df[(df["Run"] == r) & (df["Elem"] == "E")].sort_values(by="G")
        multiplier += 1

    ax1.set_title("(a) Novel Vertices")
    ax1.set_ylabel("Count")
    ax1.set_xlabel("n")
    ax1.set_xticks(x + width, Gs)

    legend_elements = [
        Patch(facecolor=colors[0][1], edgecolor=colors[0][0], label="1OUT-AUG"),
        Patch(facecolor=colors[1][1], edgecolor=colors[1][0], label="FULL-AUG"),
    ]
    ax1.legend(title="Graph", handles=legend_elements, loc="upper right")

    data = []
    for n in [1, 8]:  # NOVEL_1OUT
        for v, w in NOVEL_1OUT[n].items():
            data.append([n, "1OUT", w])
        for v, w in NOVEL_FULL[n].items():
            data.append([n, "FULL", w])
    df = pd.DataFrame(data, columns=["n", "Graph", "Support"])

    print(df[(df["n"] == 1) & (df["Graph"] == "1OUT") & (df["Support"] == 1)].shape[0])
    print(df[(df["n"] == 1) & (df["Graph"] == "1OUT")].shape[0])
    print(df[(df["n"] == 8) & (df["Graph"] == "1OUT") & (df["Support"] == 1)].shape[0])
    print(df[(df["n"] == 8) & (df["Graph"] == "1OUT")].shape[0])

    sns.histplot(
        data=df[df["n"] == 1],
        x="Support",
        hue="Graph",
        discrete=True,
        element="step",
        ax=ax2,
        legend=None,
    )
    ax2.set_ylabel("")
    ax2.set_title(f"(b) n=1")

    sns.histplot(
        data=df[df["n"] == 8],
        x="Support",
        hue="Graph",
        discrete=True,
        element="step",
        ax=ax3,
        legend=None,
    )
    ax3.set_ylim(ax2.get_ylim())
    ax3.set_ylabel("")
    ax3.set_title(f"(c) n=8")

    plt.tight_layout()
    # plt.savefig("precision.pdf")
    plt.show()


if __name__ == "__main__":
    main()
