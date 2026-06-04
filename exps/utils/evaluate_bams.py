import sys
import os
import glob
import pysam


def main():
    WD = sys.argv[1]

    print(
        "Sample",
        "Coverage",
        "Aligner",
        "RName",
        "NM",
        # "=",
        # "X",
        # "I",
        # "D",
        "Length",
        "Clipped",
        sep=",",
    )
    for bam_fn in glob.glob(os.path.join(WD, "*", "c*", "*.bam")):
        print(bam_fn, file=sys.stderr)
        fields = bam_fn.split("/")
        sample = fields[-3]
        cov = fields[-2][1:]
        aligner = fields[-1].split(".")[-2]

        for aln in pysam.AlignmentFile(bam_fn, "rb"):
            qname = aln.query_name
            if aln.is_secondary or aln.is_supplementary or aln.is_unmapped:
                continue

            if not aln.has_tag("NM"):
                continue

            nm = aln.get_tag("NM")

            # m, x, i, d = 0, 0, 0, 0
            # for q, p, c in aln.get_aligned_pairs(with_seq=True):
            #     if q != None and p != None:
            #         if c.islower():
            #             x += 1
            #         else:
            #             m += 1
            #     elif q == None:
            #         d += 1
            #     else:
            #         i += 1
            print(
                sample,
                cov,
                aligner,
                qname,
                nm,
                # m,
                # x,
                # i,
                # d,
                aln.query_length,
                aln.get_cigar_stats()[0][4],
                sep=",",
                flush=False,
            )
        sys.stdout.flush()


if __name__ == "__main__":
    main()
