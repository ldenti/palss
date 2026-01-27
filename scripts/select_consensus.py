import sys


def parse_gaf(fn, span_ratio=0.9):
    alns = {}
    reads = set()
    for line in open(fn):
        fields = line.strip("\n").split("\t")
        name, ql, qs, qe = fields[:4]
        reads.add(name)

        ql, qs, qe = int(ql), int(qs), int(qe)
        # if qs != 0 or qe != ql - 1:
        if (qe - qs) / ql < span_ratio:
            continue

        nm = int(fields[12].split(":")[-1])
        cigar = fields[16].split(":")[-1]
        alns[name] = (nm, cigar)

    return alns, reads


def main():
    original_gaf = sys.argv[1]
    augmented_gaf = sys.argv[2]
    palss_gaf = sys.argv[3]
    span_ratio = 0.9  # float(sys.argv[4])

    original, original_reads = parse_gaf(original_gaf, span_ratio)
    augmented, augmented_reads = parse_gaf(augmented_gaf, span_ratio)

    print("Original:", len(original), len(original_reads), file=sys.stderr)
    print("Augmented:", len(augmented), len(augmented_reads), file=sys.stderr)

    print(set(original) - set(augmented), file=sys.stderr)

    print("Using", len(set(original) & set(augmented)), "alignments", file=sys.stderr)
    tokeep = set()
    info = [0, 0, 0]
    for qname in set(original) & set(augmented):
        onm, ocigar = original[qname]
        anm, acigar = augmented[qname]
        if onm == 0 or onm == anm:
            i = 0
        elif onm < anm:
            # print(qname, file=sys.stderr)
            # print(original[qname], file=sys.stderr)
            # print(augmented[qname], file=sys.stderr)
            # print("", file=sys.stderr)
            i = 1
        else:
            i = 2
            tokeep.add(qname)
        info[i] += 1
    print(info, file=sys.stderr)

    rescued = 0
    for line in open(palss_gaf):
        name = line.strip("\n").split("\t")[0]
        if name in tokeep:
            print(line, end="")
        if name in augmented and name not in original:
            print(line, end="")
            rescued += 1
    print(rescued, file=sys.stderr)


if __name__ == "__main__":
    main()
