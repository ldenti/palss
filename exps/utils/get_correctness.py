import sys
import os
import glob
import pysam


def parse_gaf(gaf_fn):
    newv = {}
    v2c = {}
    c2v = {}
    for line in open(gaf_fn):
        fields = line.strip("\n").split("\t")
        name = fields[0]
        NV = set([tuple(x.split("|")) for x in fields[-1].split(":")[-1].split("/")])
        if name not in c2v:
            c2v[name] = set()
        for v, seq in NV:
            v = int(v)
            newv[v] = seq
            c2v[name].add(v)
            if v not in v2c:
                v2c[v] = set()
            v2c[v].add(name)
    return newv, v2c, c2v


def main():
    wd = sys.argv[1]

    print(
        "n",
        "Graph",
        "d",
        "w",
        "Consensus",
        "Length",
        "AlType",
        "NM",
        "M",
        "X",
        "I",
        "D",
        "NV",
        "NB",
        sep=",",
    )
    for gaf_fn in glob.glob(
        os.path.join(wd, "n*", "palss-*", "resulting-consensus.*.gaf")
    ):
        fields = gaf_fn.split("/")

        supp = int(fields[-1].split(".")[-2][1:])
        density = float(fields[-1].split(".")[-4][1:] + "." + fields[-1].split(".")[-3])
        graph = fields[-2].split("-")[1]
        n = int(fields[-3][1:])

        bam_fn = gaf_fn[:-4] + ".bam"

        newv, v2c, c2v = parse_gaf(gaf_fn)

        for aln in pysam.AlignmentFile(bam_fn, "rb"):
            qname = aln.query_name
            ql = aln.query_length
            at = "1"
            if aln.is_secondary:
                at = "2"
            elif aln.is_supplementary:
                at = "S"
            elif aln.is_unmapped:
                at = "M"

            nm = aln.get_tag("NM") if aln.has_tag("NM") else -1

            m, x, i, d = 0, 0, 0, 0
            for q, p, c in aln.get_aligned_pairs(with_seq=True):
                if q != None and p != None:
                    if c.islower():
                        x += 1
                    else:
                        m += 1
                elif q == None:
                    d += 1
                else:
                    i += 1

            nv = len(c2v[qname])
            badded = sum([len(newv[v]) for v in c2v[qname]])
            print(
                n,
                graph,
                density,
                supp,
                qname,
                ql,
                at,
                nm,
                m,
                x,
                i,
                d,
                nv,
                badded,
                sep=",",
            )


if __name__ == "__main__":
    main()
