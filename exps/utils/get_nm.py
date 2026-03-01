import sys
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


def parse_bam(bam_fn, Ns=[-1]):
    fn = bam_fn.split("/")[-1]
    for aln in pysam.AlignmentFile(bam_fn, "rb"):
        if aln.is_secondary or aln.is_supplementary or aln.is_unmapped:
            continue
        nm = aln.get_tag("NM") if aln.has_tag("NM") else -1
        for n in Ns:
            print(fn, "reference", n, "-1", "-1", aln.query_name, 1, nm, sep=",")


def main():
    WD = sys.argv[1]
    sample = sys.argv[2]

    Ns = set()
    print("fn,graph,n,w,d,read,cov,nm")
    for gaf_fn in glob.glob(os.path.join(WD, "n*", "graphaligner", "*.gaf")):
        print(gaf_fn, file=sys.stderr)
        n = int(gaf_fn.split("/")[-3][1:])
        Ns.add(n)
        fn = gaf_fn.split("/")[-1]
        w, d = -1, -1
        graph = ""
        run = "full" if "full" in fn else "oneout"
        augmentation = ""
        if "original" in fn:
            graph = "original"
        elif "mgcactus" in fn:
            graph = "mgcactus"
        else:
            graph = "palss"
            w = int(fn.split(".")[-2][1:])
            d = float(fn.split(".")[-4][1:] + "." + fn.split(".")[-3])

        nms = parse_gaf(gaf_fn)
        for qidx, (c, nm) in nms.items():
            print(fn, graph, n, w, d, qidx, c, nm, sep=",", flush=False)
        sys.stdout.flush()

    parse_bam(os.path.join(WD, f"{sample}-reads.bam"), Ns)
    parse_bam(os.path.join(WD, f"{sample}-reads.tohaps.bam"), Ns)


if __name__ == "__main__":
    main()
