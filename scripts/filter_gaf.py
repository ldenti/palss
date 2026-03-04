import sys

gaf_fn = sys.argv[1]
min_identity = float(sys.argv[2])

for line in open(gaf_fn):
    identity = float(line.split("\t")[15].split(":")[2])
    if identity >= min_identity:
        print(line, end="")
