#!/usr/bin/env python3

import sys
import gzip


def main():
    fa_fn = sys.argv[1]
    txt_fn = sys.argv[2]

    names = []
    for line in open(txt_fn):
        names.append(line.strip("\n"))

    if fa_fn[-3:] == ".gz":
        fa = gzip.open(fa_fn, "rt")
    else:
        fa = open(fa_fn)

    record = ["", ""]  # header, seq
    for i, line in enumerate(fa, 1):
        ii = i % 2
        record[ii - 1] = line.strip("\n")
        if ii == 0:
            if record[0][1:] in names:
                for x in record:
                    print(x)
    if record[0][1:] in names:
        for x in record:
            print(x)


if __name__ == "__main__":
    main()
