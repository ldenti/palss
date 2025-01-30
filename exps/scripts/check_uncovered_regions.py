import sys
from intervaltree import IntervalTree
from Bio import SeqIO
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


def main():
    fa_fn = sys.argv[1]
    segdup_fn = sys.argv[2]
    rmsk_fn = sys.argv[3]

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

    df = pd.DataFrame(d)
    sns.histplot(df, bins=100, legend=None)
    plt.xlim(0, max(d) + 25000)
    plt.xlabel("Length")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
