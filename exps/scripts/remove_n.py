#!/usr/bin/env python3

import sys
from Bio import SeqIO

def main():
    fq_path = sys.argv[1]
    for record in SeqIO.parse(fq_path, "fastq"):
        # print(record.id, sys.stderr)
        if "N" in record.seq:
            continue
        SeqIO.write(record, sys.stdout, "fastq")

if __name__ == "__main__":
    main()
