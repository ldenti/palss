#!/usr/bin/env python3

import sys
import gzip

def main():
    fq = sys.stdin
    if len(sys.argv) > 1:
        fq = sys.argv[1]
        if fq[-3:] == ".gz":
            fq = gzip.open(fq, "rt")
        else:
            fq = open(fq)

    record = ["", "", "", ""]  # header, seq, sep, quality = "", "", "", ""
    for i, line in enumerate(fq, 1):
        ii = i % 4
        record[ii - 1] = line.strip("\n")
        if ii == 0:
            if "N" in record[1]:
                continue
            for x in record:
                print(x)


if __name__ == "__main__":
    main()
