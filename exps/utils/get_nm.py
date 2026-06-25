import sys
import argparse
import os
import glob
import pysam


def parse_gaf(gaf_fn):
    nms = {}
    for line in open(gaf_fn):
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
    return nms


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("FN")  # "contigs" alignments to graph/reference
    parser.add_argument("TXT")  # "contigs" overlapping complex regions

    parser.add_argument("-t", type=str, required=True)
    parser.add_argument("-n", type=str, required=True)
    parser.add_argument("-c", type=str, required=True)
    parser.add_argument("-l", type=str, required=True)

    args = parser.parse_args()

    Cs = [int(c) for c in args.c.split(",")]

    mode = args.FN[-3:]
    assert mode in ["gaf", "bam"]

    Ns = [args.n]

    cpx_contigs = []
    for line in open(args.TXT):
        cpx_contigs.append(line.strip("\n"))

    labels = ["Simple", "Complex"]
    print("tool,n,coverage,clen,read,cov,nm,type")
    if mode == "bam":
        for aln in pysam.AlignmentFile(args.FN, "rb"):
            if aln.is_secondary or aln.is_supplementary or aln.is_unmapped:
                continue
            nm = aln.get_tag("NM") if aln.has_tag("NM") else -1
            for n in Ns:
                for cov in Cs:
                    print(
                        args.t,
                        n,
                        cov,
                        args.l,
                        aln.query_name,
                        1,
                        nm,
                        labels[int(aln.query_name in cpx_contigs)],
                        sep=",",
                    )
    else:
        nms = parse_gaf(args.FN)
        for qidx, (c, nm) in nms.items():
            # Cs is just a single value here
            for cov in Cs:
                print(
                    args.t,
                    args.n,
                    cov,
                    args.l,
                    qidx,
                    c,
                    nm,
                    labels[int(qidx in cpx_contigs)],
                    sep=",",
                    flush=False,
                )
        sys.stdout.flush()


if __name__ == "__main__":
    main()
