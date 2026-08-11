import sys
import argparse
import pysam


def parse_gaf(gaf_fn, min_c):
    nms = {}
    for line in open(gaf_fn):
        line = line.strip("\n").split("\t")

        qidx = line[0]

        ql = int(line[1])
        qs = int(line[2])  # closed
        qe = int(line[3])  # open

        c = (qe - qs) / ql

        if c < min_c:
            continue
        if len(line) < 13:
            print(gaf_fn, line[0], file=sys.stderr)
            continue
        nm = line[12]
        nm = nm.split(":")
        assert nm[0] == "NM"
        nm = int(nm[2])

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

    parser.add_argument("--gaf", nargs="+", required=True)
    parser.add_argument("--bam", nargs="*")
    parser.add_argument("-c", type=float, default=0.9)

    args = parser.parse_args()

    assert len(args.gaf) % 2 == 0 and len(args.bam) % 2 == 0

    gafs = [args.gaf[i : i + 2] for i in range(0, len(args.gaf), 2)]
    bams = [args.bam[i : i + 2] for i in range(0, len(args.bam), 2)]

    print("Reference", "Read", "Coverage", "NM", sep=",")
    for name, gaf in gafs:
        nms = parse_gaf(gaf, args.c)
        for qidx, (c, nm) in nms.items():
            print(name, qidx, c, nm, sep=",")

    for name, bam in bams:
        for aln in pysam.AlignmentFile(bam, "rb"):
            if aln.is_secondary or aln.is_supplementary or aln.is_unmapped:
                continue
            nm = aln.get_tag("NM") if aln.has_tag("NM") else -1
            c = aln.query_alignment_length / aln.infer_read_length()
            if c < args.c:
                continue
            print(
                name,
                aln.query_name,
                c,
                nm,
                sep=",",
            )


if __name__ == "__main__":
    main()
