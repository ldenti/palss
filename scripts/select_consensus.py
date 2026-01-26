import sys


def parse_gaf(fn):
    alns = {}
    reads = set()
    for line in open(fn):
        fields = line.strip("\n").split("\t")
        name, ql, qs, qe = fields[:4]
        reads.add(name)

        ql, qs, qe = int(ql), int(qs), int(qe)
        if qs != 0 or qe != ql - 1:
            continue

        nm = int(fields[12].split(":")[-1])
        cigar = fields[16].split(":")[-1]
        alns[name] = (nm, cigar)

    return alns, reads


def main():
    original_gaf = sys.argv[1]
    augmented_gaf = sys.argv[2]
    palss_gaf = sys.argv[3]

    original, original_reads = parse_gaf(original_gaf)
    augmented, augmented_reads = parse_gaf(augmented_gaf)

    print("Original:", len(original), len(original_reads), file=sys.stderr)
    print("Augmented:", len(augmented), len(augmented_reads), file=sys.stderr)

    print("Using", len(set(original) & set(augmented)), "alignments", file=sys.stderr)
    tokeep = set()
    for qname in set(original) & set(augmented):
        onm, ocigar = original[qname]
        anm, acigar = augmented[qname]
        if onm == 0 or onm == anm:
            continue
        if onm < anm:
            print(qname, file=sys.stderr)
            print(original[qname], file=sys.stderr)
            print(augmented[qname], file=sys.stderr)
            print("", file=sys.stderr)
            continue
        tokeep.add(qname)

    for line in open(palss_gaf):
        name = line.strip("\n").split("\t")[0]
        if name in tokeep:
            print(line, end="")


if __name__ == "__main__":
    main()
