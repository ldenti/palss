import sys
from Bio import SeqIO


def main():
    fx_path = sys.argv[1]
    sfs_path = sys.argv[2]

    SFS = {}
    for line in open(sfs_path):
        idx, l, strand, kept, vertices = line.strip("\n").split(" ")
        if kept == "0":
            continue
        qname = idx.split(":")[0]
        s, e = idx.split(":")[1].split("-")
        s, e = int(s), int(e)
        if qname not in SFS:
            SFS[qname] = []
        SFS[qname].append((s, e, strand == "1"))

    fq = False
    for record in SeqIO.parse(fx_path, "fasta"):
        idx = record.id
        if idx not in SFS:
            continue
        sfs = SFS[idx]
        sfs.sort(key=lambda x: x[0])
        for s, e, plus in sfs:
            if fq:
                pass
                # print(f"@{idx}:{s}-{e}")
                # print(record.seq[s:e])
                # print("+")
                # print(
                #     "".join(
                #         [
                #             chr(q + 33)
                #             for q in record.letter_annotations["phred_quality"][s:e]
                #         ]
                #     )
                # )
            else:
                print(f">{idx}:{s}-{e} {plus}")
                if plus:
                    print(record.seq[s:e])
                else:
                    print(record.seq[s:e].reverse_complement())


if __name__ == "__main__":
    main()
