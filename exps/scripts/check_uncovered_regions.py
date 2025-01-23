import sys
from intervaltree import IntervalTree

from Bio import SeqIO


def main():
    fa_fn = sys.argv[1]
    segdup_fn = sys.argv[2]
    rmsk_fn = sys.argv[3]
    CHROM = "chr1"

    n = 0
    sd_trees = {}
    for line in open(segdup_fn):
        chrom, s, e = line.split("\t")[:3]
        # if chrom != CHROM:
        #     continue
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
        # if chrom != CHROM:
        #     continue
        repclass = "RM"  # line.split("\t")[11]
        s, e = int(s), int(e)
        if chrom not in rm_trees:
            rm_trees[chrom] = IntervalTree()
        rm_trees[chrom][s:e] = repclass
        n += 1
        if n % 5000 == 0:
            print(f"Loaded {n} repeats..", end="\r", file=sys.stderr)
    print(f"Loaded {n} repeats..", end="\n", file=sys.stderr)

    for record in SeqIO.parse(fa_fn, "fasta"):
        chrom, se = record.id.split(":")
        # if chrom != CHROM:
        #     continue
        s, e = (int(x) for x in se.split("-"))
        nn = record.count("N")
        if nn > 0.8:
            print(record.id, f"N{nn/len(record)*100}", f"{nn}/{len(record)}")
        else:
            if sd_trees[chrom].overlaps(s, e):
                print(record.id, "SEGDUP")
            else:
                overlaps = rm_trees[chrom].overlap(s, e)
                if len(overlaps) != 0:
                    repclasses = set()
                    for interval in overlaps:
                        repclasses.add(interval.data)
                    print(record.id, ",".join(repclasses))
                else:
                    print(record.id, "False")


if __name__ == "__main__":
    main()
