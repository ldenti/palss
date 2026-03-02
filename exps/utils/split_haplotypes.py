import sys
import gzip


def split(header, seq, size=50000):
    l = len(seq)
    if l < size:
        print(f">{header}.0")
        print(seq)
    else:
        i = 0
        p = 0
        while p < l:
            print(f">{header}.{i}")
            print(seq[p : p + size])
            p += int(size / 2)
            i += 1


def main():
    fa_fn = sys.argv[1]

    if fa_fn[-3:] == ".gz":
        fa = gzip.open(fa_fn, "rt")
    else:
        fa = open(fa_fn)

    header = ""
    seq = ""
    for line in fa:
        line = line.strip("\n")
        if line[0] == ">":
            if seq != "":
                split(header, seq)
            header = line[1:]
            seq = ""
        else:
            seq += line
    #
    split(header, seq)


if __name__ == "__main__":
    main()
