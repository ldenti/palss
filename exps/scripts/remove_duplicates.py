import sys
from pysam import VariantFile


def get_var(record):
    chrom = record.contig
    pos = record.pos
    alleles = "/".join(record.alleles)
    alleles = alleles.upper()
    return f"{chrom}:{pos}:{alleles}"


def main():
    vcf_fn = sys.argv[1] if len(sys.argv) > 1 else sys.stdin
    vcf = VariantFile(vcf_fn)

    print(vcf.header, end="")
    kept = set()
    tot=0
    filt=0
    for record in vcf:
        v = get_var(record)
        tot += 1
        if v in kept:
            v = ":".join(v.split(":")[:-1])
            print(f"Skipping {v}", file=sys.stderr)
            filt += 1
            continue
        record.ref = record.ref.upper()
        record.alts = (x.upper() for x in record.alts)
        print(record, end="")
        kept.add(v)
    print(filt, "/", tot, file=sys.stderr)

if __name__ == "__main__":
    main()
