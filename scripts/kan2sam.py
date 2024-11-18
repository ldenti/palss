import sys

from Bio import SeqIO


def main():
    fa_fn = sys.argv[1]
    txt_fn = sys.argv[2]

    print("@HD", "VN:1.6", "SO:coordinate", sep="\t")
    for record in SeqIO.parse(fa_fn, "fasta"):
        print("@SQ", f"SN:{record.id}", f"LN:{len(record.seq)}", sep="\t")

    regions = {}
    i = 0
    for line in open(txt_fn):
        region = line.split(" ")[0]
        kmer = line.strip("\n").split(" ")[-1]
        chrom = region.split(":")[0]
        s, e = region.split(":")[1].split("-")
        print(i, 0, chrom, s, 60, f"{len(kmer)}M", "*", 0, 0, kmer, "*", sep="\t")
        i += 1


if __name__ == "__main__":
    main()
