import sys
import argparse

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt


def main(args):
    chroms = []
    if args.chroms != "":
        chroms = args.chroms.split(",")

    of = open(args.out, "w") if args.out != "" else sys.stdout

    sizes = {}
    for line in open(args.FAI):
        chrom, size, _, _, _ = line.split("\t")
        sizes[chrom] = int(size)

    regions = {}
    for line in open(args.BED):
        chrom, s, e, _idx = line.split("\t")
        if len(chroms) != 0 and chrom not in chroms:
            continue
        if chrom not in regions:
            regions[chrom] = []
        regions[chrom].append((int(s), int(e) - 1))  # 1-based, closed

    data = []
    overlapping = 0
    consecutive = 0
    uncovered = {}
    for chrom, kmers in regions.items():
        uncovered[chrom] = 0
        kmers.sort(key=lambda x: x[0])
        for (s1, e1), (s2, e2) in zip(kmers[:-1], kmers[1:]):
            if s2 < e1:
                overlapping += 1
            else:
                d = s2 - e1
                if d == 0:
                    consecutive += 1
                else:
                    if d > args.D:
                        uncovered[chrom] += d
                        print(
                            f"# {chrom}:{e1+1}-{s2} ({d}) {chrom}:{s1}-{e1+1} {chrom}:{s2}-{e2+1}",
                            file=of,
                        )
                    data.append([chrom, d])
    total = overlapping + consecutive + len(data)

    for chrom, size in sizes.items():
        ratio = uncovered[chrom]/size if chrom in uncovered else -1
        print(chrom, ratio, sep="\t", file = of)

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

    df = pd.DataFrame(data, columns=["Path", "d"])
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

    sns.histplot(
        df,
        x="d",
        binrange=(1, df["d"].quantile(q=0.95)),
        bins=max(1, int(df["d"].quantile(q=0.95) / 50)),
        legend=None,
    )  # , hue="Path")
    if args.out == "":
        plt.show()
    else:
        plt.savefig(args.out + ".png")

    of.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="kan_hist",
        description="Analyze anchors distance from BED file",
    )
    parser.add_argument("BED", help="List of anchors (.bed)")
    parser.add_argument("FAI", help="Index of reference file (.fai)")
    parser.add_argument(
        "-c",
        help="Comma-separated list of chromosomes to consider (default: all)",
        dest="chroms",
        type=str,
        required=False,
        default="",
    )
    parser.add_argument(
        "-o",
        help="Comma-separated list of chromosomes to consider (default: stdout/show plot)",
        dest="out",
        type=str,
        required=False,
        default="",
    )
    parser.add_argument(
        "-D",
        help="Print all pairs whose distance is greater than this (default: 15000)",
        dest="D",
        type=int,
        required=False,
        default=15000,
    )
    parser.add_argument(
        "-b",
        help="Bin size (default: 50)",
        dest="bsize",
        type=int,
        required=False,
        default=50,
    )
    parser.add_argument(
        "-q",
        help="Quartile (default: 0.95)",
        dest="quart",
        type=float,
        required=False,
        default=0.95,
    )

    args = parser.parse_args()
    main(args)
