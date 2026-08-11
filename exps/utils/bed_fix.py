import sys

map_fn = sys.argv[1]
bed_fn = sys.argv[2]

mapping = {}
for line in open(map_fn):
    x, y = line.strip("\n").split("\t")
    mapping[y] = f"chr{x}"

for line in open(bed_fn):
    fields = line.strip("\n").split("\t")
    if fields[0] not in mapping:
        continue
    fields[0] = mapping[fields[0]]
    print("\t".join(fields))
