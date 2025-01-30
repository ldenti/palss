import sys
import os
import glob
import argparse

from matplotlib import pyplot as plt

import pandas as pd
import seaborn as sns


def main():
    fai = sys.argv[1]
    wd = sys.argv[2]
    dt = int(sys.argv[3])

    sizes = {}
    for line in open(fai):
        chrom, size, _, _, _ = line.split("\t")
        sizes[chrom] = int(size)

    data = []
    Gs = []
    for bed in glob.glob(os.path.join(wd, "*", "k*", "reference-anchors.bed")):
        info = bed.split("/")
        g = int(info[-3])
        # if g not in [1,16]:
        #    continue
        print("===", g, file=sys.stderr)
        Gs.append(g)

        last_chrom = ""
        regions = []
        for line in open(bed):
            chrom, s, e, _idx = line.split("\t")
            if "_" in chrom:
                continue
            if last_chrom != "" and chrom != last_chrom:
                print(f"Analyzing {last_chrom}, next chrom: {chrom}", file=sys.stderr)
                close = 0
                total = 0
                for (s1, e1), (s2, e2) in zip(regions[:-1], regions[1:]):
                    total += 1
                    d = s2 - e1
                    if d <= 100:
                        close += 1
                        continue
                    if d > dt:
                        print(
                            f"# {g} {last_chrom}:{e1+1}-{s2} ({d}) {last_chrom}:{s1}-{e1+1} {last_chrom}:{s2}-{e2+1}"
                        )
                    data.append([g, last_chrom, d])
                print(f"{close} / {total} = {close/total}", file=sys.stderr)
                regions = []
            last_chrom = chrom
            regions.append((int(s), int(e) - 1))  # 1-based, closed
        if len(regions) > 0:
            print(f"Analyzing {last_chrom}, no next chrom", file=sys.stderr)
            total += 1
            d = s2 - e1
            for (s1, e1), (s2, e2) in zip(regions[:-1], regions[1:]):
                total += 1
                d = max(0, s2 - e1)
                if d <= 100:
                    close += 1
                    continue
                if d > dt:
                    print(
                        f"# {g} {last_chrom}:{e1+1}-{s2} ({d}) {last_chrom}:{s1}-{e1+1} {last_chrom}:{s2}-{e2+1}"
                    )
                data.append([g, last_chrom, d])
            print(f"{close} / {total} = {close/total}", file=sys.stderr)
            regions = []

    Gs.sort()

    df = pd.DataFrame(data, columns=["G", "Chrom", "d"])
    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(
        2, 3, figsize=(9, 7), sharex=True, sharey=True
    )
    for g, ax in zip(Gs, [ax1, ax2, ax3, ax4, ax5, ax6]):
        sdf = df[df["G"] == g]
        print(sdf["d"].describe(), file=sys.stderr)
        sns.histplot(
            sdf,
            x="d",
            binrange=(1, sdf["d"].quantile(q=0.90)),
            # bins=max(1, int(df["d"].quantile(q=0.98) / 50)),
            legend=None,
            ax=ax,
        )
        ax.set_title(f"{g}")
    plt.tight_layout()
    plt.show()
    plt.savefig("d-hist.png")
    plt.close()

    data = []
    for g in Gs:
        for chrom in df["Chrom"].unique():
            sdf = df[(df["G"] == g) & (df["Chrom"] == chrom)]
            notok = sum(sdf[sdf["d"] > dt]["d"])
            data.append([g, chrom, notok, sizes[chrom]])

    df = pd.DataFrame(data, columns=["G", "Chrom", "NotOk", "Size"])
    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(
        2, 3, figsize=(9, 7), sharex=True, sharey=True
    )
    for g, ax in zip(Gs, [ax1, ax2, ax3, ax4, ax5, ax6]):
        sns.barplot(data=df[df["G"] == g], x="Chrom", y="Size", ax=ax)
        sns.barplot(data=df[df["G"] == g], x="Chrom", y="NotOk", ax=ax)
        ax.set_title(f"{g}")
    plt.tight_layout()
    plt.show()
    plt.savefig("d-bar.png")


def analyze(chrom, regions, ofp):
    of = open(f"{ofp}.{chrom}.txt", "w") if ofp != "" else sys.stdout

    data = []
    overlapping = 0
    consecutive = 0
    uncovered = 0
    regions.sort(key=lambda x: x[0])
    for (s1, e1), (s2, e2) in zip(regions[:-1], regions[1:]):
        if s2 < e1:
            overlapping += 1
        else:
            d = s2 - e1
            if d == 0:
                consecutive += 1
            else:
                if d > args.D:
                    uncovered += d
                    print(
                        f"# {chrom}:{e1+1}-{s2} ({d}) {chrom}:{s1}-{e1+1} {chrom}:{s2}-{e2+1}",
                        file=of,
                    )
                data.append(d)
    total = overlapping + consecutive + len(data)

    print(
        "Overlapping over total anchors:",
        overlapping,
        "/",
        total,
        overlapping / total,
        file=of,
    )
    print(
        "Consecutive over total anchors:",
        consecutive,
        "/",
        total,
        consecutive / total,
        file=of,
    )
    print("Data size:", len(data), file=of)

    df = pd.DataFrame(data, columns=["d"])
    print(f"=== {chrom} ===", file=of)
    print(
        df.describe(
            percentiles=[
                0.1,
                0.25,
                0.4,
                0.5,
                0.6,
                0.75,
                0.90,
                0.95,
                0.96,
                0.97,
                0.98,
                0.99,
            ]
        ),
        file=of,
    )
    if ofp != "":
        of.close()

    sns.histplot(
        df,
        x="d",
        binrange=(1, df["d"].quantile(q=0.98)),
        # bins=max(1, int(df["d"].quantile(q=0.98) / 50)),
        legend=None,
    )
    plt.title(chrom)
    if ofp == "":
        plt.show()
    else:
        plt.savefig(f"{ofp}.{chrom}.png")
    plt.close()

    return uncovered


def main2(args):
    chroms = []
    if args.chroms != "":
        chroms = args.chroms.split(",")

    sizes = {}
    for line in open(args.FAI):
        chrom, size, _, _, _ = line.split("\t")
        if "_" in chrom:
            continue
        sizes[chrom] = int(size)

    uncovered = {}
    last_chrom = ""
    regions = []
    for line in open(args.BED):
        chrom, s, e, _idx = line.split("\t")
        if "_" in chrom:
            continue
        if len(chroms) != 0 and chrom not in chroms:
            continue
        if chrom != last_chrom:
            if last_chrom != "":
                print(f"Analyzing {last_chrom}, next chrom: {chrom}", file=sys.stderr)
                uncovered[last_chrom] = analyze(last_chrom, regions, args.out)
            last_chrom = chrom
            regions = []
        regions.append((int(s), int(e) - 1))  # 1-based, closed

    if len(regions) != 0:
        uncovered[last_chrom] = analyze(last_chrom, regions, args.out)

    of = open(f"{args.out}.uncovered.txt", "w") if args.out != "" else sys.stdout

    uncovered_ratios = []
    for chrom in uncovered:
        ratio = uncovered[chrom] / sizes[chrom]
        print(chrom, uncovered[chrom], sizes[chrom], ratio, sep="\t", file=of)
        uncovered_ratios.append([chrom, uncovered[chrom], sizes[chrom], ratio])
    if args.out != "":
        of.close()

    # chrom_order = []
    # not_number = []
    # for c in uncovered_ratios:
    #     c = c[0]
    #     print("///", c)
    #     if c[0] != "c" and c.isdigit():
    #         chrom_order.append(int(c))
    #     elif c[0] == "c" and c[3:].isdigit():
    #         chrom_order.append(int(c[3:]))
    #     else:
    #         not_number.append(c[3:])
    # chrom_order.sort()
    # chrom_order = [f"chr{c}" for c in chrom_order] + [f"chr{c}" for c in not_number]
    # print(chrom_order)

    df = pd.DataFrame(uncovered_ratios, columns=["Chrom", "Unc.", "Size", "Ratio"])
    sns.barplot(df, x="Chrom", y="Size")  # , order=chrom_order)
    sns.barplot(df, x="Chrom", y="Unc.")  # , order=chrom_order)
    plt.tick_params("x", labelrotation=45)

    if args.out == "":
        plt.show()
    else:
        plt.savefig(f"{args.out}.png")


if __name__ == "__main__":
    # parser = argparse.ArgumentParser(
    #     prog="kan_hist",
    #     description="Analyze anchors distance from BED file",
    # )
    # parser.add_argument("BED", help="List of anchors (.bed)")
    # parser.add_argument("FAI", help="Index of reference file (.fai)")
    # parser.add_argument(
    #     "-c",
    #     help="Comma-separated list of chromosomes to consider (default: all)",
    #     dest="chroms",
    #     type=str,
    #     required=False,
    #     default="",
    # )
    # parser.add_argument(
    #     "-o",
    #     help="Comma-separated list of chromosomes to consider (default: stdout/show plot)",
    #     dest="out",
    #     type=str,
    #     required=False,
    #     default="",
    # )
    # parser.add_argument(
    #     "-D",
    #     help="Print all pairs whose distance is greater than this (default: 15000)",
    #     dest="D",
    #     type=int,
    #     required=False,
    #     default=15000,
    # )
    # parser.add_argument(
    #     "-b",
    #     help="Bin size (default: 50)",
    #     dest="bsize",
    #     type=int,
    #     required=False,
    #     default=50,
    # )
    # parser.add_argument(
    #     "-q",
    #     help="Quartile (default: 0.95)",
    #     dest="quart",
    #     type=float,
    #     required=False,
    #     default=0.95,
    # )

    # args = parser.parse_args()
    main()
