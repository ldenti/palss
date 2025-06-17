import sys
from Bio import SeqIO


gaf_fn = sys.argv[1]
fa_fn = sys.argv[2]

alignments = {}
for line in open(gaf_fn):
    idx, ql, qs, qe, strand, path, pl, ps, pe, residues, blocks, mapq, nm, *rest = (
        line.strip("\n").split("\t")
    )

    cidx, gt, l, w, a1, a2 = idx.split(".")
    v1, off1 = [int(x) for x in a1.split(":")]
    v2, off2 = [int(x) for x in a2.split(":")]

    nm = int(nm.split(":")[2])

    if path[0] == "<":
        continue
    path = [int(x) for x in path[1:].split(">")]
    if path[0] != v1 or path[-1] != v2:
        continue
    if idx in alignments:
        alignments[idx] = [
            min(nm, alignments[idx][0]),
            ",".join([str(x) for x in path]),
        ]
    else:
        alignments[idx] = [nm, ",".join([str(x) for x in path])]

for record in SeqIO.parse(fa_fn, "fasta"):
    if record.id in alignments:
        if alignments[record.id][0] > 0:
            print(record.id, alignments[record.id][1], str(record.seq))
