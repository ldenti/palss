import sys
import os
import glob

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from intervaltree import IntervalTree


sns.set(style="whitegrid")


def main_fq():
    txt_fn = sys.argv[1]
    df = []
    for line in open(txt_fn):
        _, missed, length = line.strip("\n").split("\t")
        df.append([int(length) - int(missed), int(length)])

    df = pd.DataFrame(df, columns=["Anchors (#)", "Read Length (bp)"])
    print(df.describe())

    anchors_0 = len(df[df["Anchors (#)"] == 0])
    anchors_1 = len(df[df["Anchors (#)"] == 1])

    print("Reads   :", len(df))
    print("0 1 01  :", anchors_0, anchors_1, anchors_0 + anchors_1)
    print("01/reads:", (anchors_0 + anchors_1) / len(df))
    for q in [0.01, 0.02]:
        print(f"q{q}   :", df["Anchors (#)"].quantile(q))

    # sns.jointplot(
    #     y="Anchors (#)",
    #     x="Read Length (bp)",
    #     data=df[df["Read Length (bp)"] <= 35000],
    #     color="seagreen",
    #     alpha=0.33,
    #     s=3,
    # )
    # plt.tight_layout()
    # plt.show()


def main_fa():
    bed_fn = sys.argv[1]
    fai_fn = sys.argv[2]
    complex_fn = sys.argv[3]

    min_gap = 10000  # XXX: hardcoded

    print("Preparing chromosomes...")
    data = {}
    for line in open(fai_fn):
        cname, l = line.split("\t")[:2]
        if cname == "chrM":
            continue
        data[cname] = [0, int(l)]

    print("Loading complex regions...")
    complex_trees = {}
    for line in open(complex_fn):
        chrom, s, e, *_ = line.strip("\n").split("\t")
        s, e = int(s), int(e)
        if chrom not in complex_trees:
            complex_trees[chrom] = IntervalTree()
        complex_trees[chrom][s:e] = 1

    print("Loading missed regions...")
    missed = []
    for line in open(bed_fn):
        cname, s, e, *_ = line.strip("\n").split("\t")
        s, e = int(s), int(e)
        l = e - s
        if l > min_gap:
            flag = "Simple"
            if complex_trees[chrom].overlaps(s, e):
                flag = "Complex"
            missed.append([l, flag])
            data[cname][0] += l

    print("Plotting...")
    data = [[k] + v for k, v in data.items()]
    df = pd.DataFrame(data, columns=["Chromosome", "Missed", "Length"])
    df["%Missed"] = df["Missed"] / df["Length"]
    print(df)

    total_missed = df["Missed"].sum()
    total_length = df["Length"].sum()
    print(
        "Full genome:",
        total_missed,
        total_length,
        round(total_missed / total_length * 100, 3),
    )

    total_complex = sum([x[1] == "Complex" for x in missed])
    print("Complex:", total_complex, len(missed), total_complex / len(missed))
    missed = pd.DataFrame(missed, columns=["Length", "Type"])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 6))

    sns.histplot(
        data=missed,
        x="Length",
        hue="Type",
        multiple="stack",
        binrange=(0, 150000),
        ax=ax1,
    )
    #
    sns.barplot(data=df, x="Chromosome", y="Length", color="black", ax=ax2)
    sns.barplot(data=df, x="Chromosome", y="Missed", color="red", ax=ax2)

    xticks = ax2.get_xticks()
    xticklabels = [t.get_text() for t in ax2.get_xticklabels()]

    labels = {}
    for t in df.itertuples():
        labels[t.Chromosome] = round(t._4 * 100, 2)

    for p in ax2.patches[len(data) :]:
        x = p.get_x() + p.get_width() / 2 + 0.1
        height = p.get_height()

        closest_idx = (abs(xticks - x)).argmin()
        xticklabel = xticklabels[closest_idx]
        label = labels[xticklabel]
        ax2.annotate(
            str(label),
            xy=(x, height),
            xytext=(0, 6),
            textcoords="offset points",
            ha="center",
            va="bottom",
            color="red",
            rotation=90,
            fontsize=9,
        )

    ax1.set_title("(a)")
    ax2.set_title("(b)")
    for label in ax1.get_xticklabels():
        label.set_rotation(30)
    for label in ax2.get_xticklabels():
        label.set_rotation(45)

    plt.tight_layout()
    plt.show()


def main_merge():
    wd = sys.argv[1]
    fai_fn = sys.argv[2]

    data = {}
    for line in open(fai_fn):
        cname, l = line.split("\t")[:2]
        data[cname] = [int(l), {}]

    Ds = set()
    for bed_fn in glob.glob(os.path.join(wd, "d*", "missed_regions.g15000.bed")):
        d = float(bed_fn.split("/")[-2][1:])
        Ds.add(d)
        for line in open(bed_fn):
            cname, s, e, *_ = line.strip("\n").split("\t")
            if d not in data[cname][1]:
                data[cname][1][d] = 0
            data[cname][1][d] += int(e) - int(s)
    print(data)

    Ds = sorted(list(Ds))
    print(Ds)

    df = []
    for chrom, (chrom_l, missed) in data.items():
        if chrom == "chrM":
            continue
        row = [chrom, chrom_l] + [missed[d] for d in Ds]
        df.append(row)

    df = pd.DataFrame(
        df,
        columns=[
            "Chromosome",
            "Length",
            "d0.1 (%)",
            "d0.25 (%)",
            "d0.33 (%)",
            "d0.5 (%)",
        ],
    )

    tot = df.select_dtypes(include="number").sum()
    tot["Chromosome"] = "Genome"
    df.loc[len(df)] = tot

    for d in Ds:
        d = f"d{d} (%)"
        df[d] = (df[d] / df["Length"] * 100).round(2)

    print(
        df.to_latex(
            index=False,
            float_format="%.2f",
        )
    )


if __name__ == "__main__":
    mode = sys.argv.pop(1)
    if mode == "fa":
        main_fa()
    elif mode == "fq":
        main_fq()
    elif mode == "merge":
        main_merge()
