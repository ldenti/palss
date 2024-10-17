import sys


def main():
    gfa_fn = sys.argv[1]
    pgcf_fn = sys.argv[2]
    # PP = 0 # sys.argv[3]  # which path to use as reference

    vertices = {}
    name = ""
    pathl = 0
    path = []
    # p = 0
    for line in open(gfa_fn):
        if line.startswith("S"):
            _, idx, seq, *_ = line.strip("\n").split("\t")
            vertices[idx] = len(seq)
        elif line[0] in ["W", "P"]:
            if True:  # p == PP:
                if line[0] == "W":
                    pass
                    # _, name, hap, chrom, _, _, path, *_ = line.strip("\n").split("\t")
                    # name = name + "#" + hap
                    # s = path[0]
                    # path = path[1:].split(s)
                    # CHECKME: do we need to reverse the path if < ?
                else:
                    _, name, path, *_ = line.strip("\n").split("\t")
                    s = path[-1] == "+"
                    path = [x[:-1] for x in path.split(",")]
                    if not s:
                        path = path[::-1]
                pathl = sum([vertices[x] for x in path])
            # p += 1

    print("##fileformat=VCFv4.2")
    print('##FILTER=<ID=PASS,Description="All filters passed">')
    print(f"##contig=<ID={name},length={pathl}>")
    print('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Variant type">')
    print('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Variant size">')
    print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    print(
        "#CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "INFO",
        "FORMAT",
        "sample",
        sep="\t",
    )
    for line in open(pgcf_fn):
        idx, hap, svtype, svlen, pname, subpath, qual, ref, alt, offset, cigar = line.strip(
            "\n"
        ).split("\t")
        if pname != name:
            continue
        svlen = "-" + svlen if svtype == "DEL" else svlen
        sv = subpath.split(">")[0]
        pos = 0
        for v in path:
            if v == sv:
                break
            pos += vertices[v]
        print(
            pname,
            pos + int(offset),
            idx,
            ref,
            alt,
            qual,
            "PASS",
            f"SVTYPE={svtype};SVLEN={svlen}",
            "GT",
            "./.",
            sep="\t",
        )
        # print(line)


if __name__ == "__main__":
    main()
