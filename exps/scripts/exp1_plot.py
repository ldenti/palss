import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from intervaltree import IntervalTree
from Bio import SeqIO

import matplotlib.gridspec as gridspec
import SeabornFig2Grid as sfg

sns.set(style="whitegrid")


def main():
    fa_fn = sys.argv[1]
    segdup_fn = sys.argv[2]
    rmsk_fn = sys.argv[3]

    txt_fn = sys.argv[4]
    fai_fn = sys.argv[5]
    # l = int(sys.argv[6])

    n = 0
    sd_trees = {}
    for line in open(segdup_fn):
        chrom, s, e = line.split("\t")[:3]
        s, e = int(s), int(e)
        if chrom not in sd_trees:
            sd_trees[chrom] = IntervalTree()
        sd_trees[chrom][s:e] = 1
        n += 1
        if n % 10000 == 0:
            print(f"Loaded {n} segsups..", end="\r", file=sys.stderr)
    print(f"Loaded {n} segsups..", end="\n", file=sys.stderr)

    n = 0
    rm_trees = {}
    for line in open(rmsk_fn):
        chrom, s, e = line.split("\t")[:3]
        repclass = line.split("\t")[3]
        s, e = int(s), int(e)
        if chrom not in rm_trees:
            rm_trees[chrom] = IntervalTree()
        rm_trees[chrom][s:e] = repclass
        n += 1
        if n % 5000 == 0:
            print(f"Loaded {n} repeats..", end="\r", file=sys.stderr)
    print(f"Loaded {n} repeats..", end="\n", file=sys.stderr)

    d = []
    for record in SeqIO.parse(fa_fn, "fasta"):
        chrom, se = record.id.split(":")
        s, e = (int(x) for x in se.split("-"))
        d.append(e - s + 1)
        nn = record.count("N")
        if nn > 0.8:
            print(record.id, e - s + 1, f"N{nn/len(record)*100}", f"{nn}/{len(record)}")
        else:
            if sd_trees[chrom].overlaps(s, e):
                print(record.id, e - s + 1, "SEGDUP")
            else:
                overlaps = rm_trees[chrom].overlap(s, e)
                if len(overlaps) != 0:
                    repclasses = set()
                    for interval in overlaps:
                        repclasses.add(interval.data)
                    print(record.id, e - s + 1, ",".join(repclasses))
                else:
                    print(record.id, e - s + 1, "False")

    fig = plt.figure(figsize=(11, 5))

    df = pd.DataFrame(d)

    g1 = sns.displot(data=df, bins=100, legend=None)
    plt.title("(a)")
    g1.set_xlabels("Region Length (bp)")

    g1.set(
        xticks=[25000, 50000, 75000, 100000, 125000, 150000, 175000, 200000],
        xticklabels=["25k", "50k", "75k", "100k", "125k", "150k", "175k", "200k"],
    )

    data = {}
    for line in open(fai_fn):
        ridx, l = line.split("\t")[:2]
        l = int(l)
        data[ridx] = [0, l]

    for line in open(txt_fn):
        ridx, anchors, total = line.strip("\n").split("\t")
        anchors = int(anchors)
        if ridx in data:
            data[ridx][0] = anchors

    df = pd.DataFrame(list(data.values()), columns=["Anchors (#)", "Read Length (bp)"])
    print(df.describe())

    anchors_0 = len(df[df["Anchors (#)"] == 0])
    anchors_1 = len(df[df["Anchors (#)"] == 1])

    print("Reads   :", len(df))
    print("0 1 01  :", anchors_0, anchors_1, anchors_0 + anchors_1)
    print("01/reads:", (anchors_0 + anchors_1) / len(df))
    print("q0.01   :", df["Anchors (#)"].quantile(0.01))
    print("q0.02   :", df["Anchors (#)"].quantile(0.02))

    g2 = sns.jointplot(
        y="Anchors (#)",
        x="Read Length (bp)",
        data=df[df["Read Length (bp)"] <= 25000],
        color="seagreen",
        alpha=0.33,
        s=3,
    )

    g2.fig.axes[1].set_title("(b)")

    gs = gridspec.GridSpec(1, 2)
    mg0 = sfg.SeabornFig2Grid(g1, fig, gs[0])
    mg1 = sfg.SeabornFig2Grid(g2, fig, gs[1])
    gs.tight_layout(fig)
    plt.savefig("x1.pdf", bbox_inches="tight")
    plt.show()
    plt.close()

    sns.jointplot(
        y="Anchors (#)",
        x="Read Length (bp)",
        data=df,
        color="seagreen",
        alpha=0.33,
        s=3,
    )
    plt.tight_layout()
    # plt.savefig("x2.pdf", bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()
