import sys
from pysam import VariantFile


def v2uidx(v):
    return f"{v.contig}.{v.pos}.{'-'.join(v.alleles)}"

def main():
    vcf1_fn = sys.argv[1]
    vcf2_fn = sys.argv[2]
    sample = sys.argv[3]

    vcf2 = []
    for record in VariantFile(vcf2_fn):
        vcf2.append(v2uidx(record))
        assert sample == record.samples[0].name

    vcf1 = VariantFile(vcf1_fn)
    print(vcf1.header, end="")
    for record in vcf1:
        uidx = v2uidx(record)
        if uidx not in vcf2:
            print(record, end="")

if __name__ == "__main__":
    main()
