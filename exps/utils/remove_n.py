#!/usr/bin/env python3

import sys


def main():
    fq = open(sys.argv[1]) if len(sys.argv) > 1 else sys.stdin
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
