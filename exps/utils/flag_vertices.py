import sys
import re
import random
from collections import deque
from intervaltree import IntervalTree, Interval


def parse_bed(fn, flank=0):
    # XXX: assuming merged BED
    trees = {}
    for line in open(fn):
        chrom, s, e, *info = line.strip("\n").split("\t")
        s, e = int(s), int(e)
        s, e = s - flank, e + flank
        if chrom not in trees:
            trees[chrom] = IntervalTree()
        trees[chrom][s:e] = True
    return trees


def load_graph(gfa_fn, ref_name):
    segments_lengths = {}
    in_links = {}
    out_links = {}
    ref_offsets = {}

    for line in open(gfa_fn):
        fields = line.strip("\n").split("\t")
        if fields[0] == "S":
            _, idx, seq, *rest = fields
            idx = int(idx)
            segments_lengths[idx] = len(seq)
            if idx not in in_links:
                in_links[idx] = []
            if idx not in out_links:
                out_links[idx] = []
        elif fields[0] == "L":
            _, idx1, _, idx2, _, _ = fields
            idx1, idx2 = int(idx1), int(idx2)
            if idx2 not in in_links:
                in_links[idx2] = []
            if idx1 not in out_links:
                out_links[idx1] = []
            in_links[idx2].append(idx1)
            out_links[idx1].append(idx2)

    segments = set()
    for line in open(gfa_fn):
        fields = line.strip("\n").split("\t")
        if fields[0] not in ["P", "W"]:
            continue
        if ref_name not in fields[1]:
            continue
        path = []
        if line.startswith("P"):
            assert False
            path = [int(x[:-1]) for x in fields[2].split(",")]
        else:
            path = [int(x) for x in re.split("[<>]", fields[6][1:])]

        contig = fields[3]
        ref_offsets[contig] = {}
        offset = 0
        for idx in path:
            segments.add(idx)
            ref_offsets[contig][idx] = offset
            offset += segments_lengths[idx]

    return segments, in_links, out_links, ref_offsets


def bfs_find_first_flagged(v0, segments, links, max_v=500):
    """
    BFS from v; returns the first flagged vertex, or None if none reachable. Visits each node at most once.
    """
    seen = {v0}
    q = deque([v0])
    while q and len(seen) < max_v:
        u = q.popleft()
        for v in links.get(u, ()):
            if v in seen:
                continue
            if v in segments:
                return v
            seen.add(v)
            q.append(v)
    return None


def get_region(segments, in_links, out_links, ref_offsets, idx, size=5000):
    contig, start, end = "", -1, -1
    rsv, rev = None, None

    rsv = bfs_find_first_flagged(idx, segments, in_links, size)
    if rsv != None:
        rev = bfs_find_first_flagged(idx, segments, out_links, size)

    if rsv != None and rev != None:
        for c, offsets in ref_offsets.items():
            if rsv in offsets and rev in offsets:
                contig = c
                start = offsets[rsv]
                end = offsets[rev]
                break
        assert contig != ""
    return contig, start, end


def flag(tree, start, end):
    label = ""
    if tree.overlaps(start, end):
        return "Complex"
    else:
        return "Simple"
    #     print(
    #         line, f"{contig}:{start}-{end}", end - start, f"{rsv}>{rev}", t, sep=","
    #     )


def main():
    gfa_fn = sys.argv[1]
    support_table = sys.argv[2]
    bed_fn = sys.argv[3]

    print("Loading regions..", file=sys.stderr)
    trees = parse_bed(bed_fn)
    print("Regions loaded!", file=sys.stderr)

    print("Loading graph..", file=sys.stderr)
    segments, in_links, out_links, ref_offsets = load_graph(
        gfa_fn, "CHM13"
    )  # XXX: hardcoded
    print("Graph loaded!", file=sys.stderr)

    print("Iterating over novel vertices..", file=sys.stderr)

    regions_memoization = {}
    for line in open(support_table):
        if line.startswith("tool"):
            print(line.strip("\n"), "region", "type", sep=",")
            continue

        line = line.strip("\n")
        fields = line.split(",")
        contig, start, end = "", -1, -1
        label, region = "", ""

        if fields[4] == "vertex":
            idx = int(fields[5])
            if idx in regions_memoization:
                label, region = regions_memoization[idx]
            else:
                assert idx not in segments
                contig, start, end = get_region(
                    segments, in_links, out_links, ref_offsets, idx
                )
                if contig != "":
                    label = flag(trees[contig], start, end)
                    region = f"{contig}:{start}-{end}"
                else:
                    label = "Unknown"
                    region = "."
                regions_memoization[idx] = (label, region)
        else:
            v1, v2 = fields[5].split(">")
            v1, v2 = int(v1[:-1]), int(v2[:-1])

            contig1, start1, end1 = "", -1, -1
            contig2, start2, end2 = "", -1, -1
            label1, label2 = "", ""
            region1, region2 = "", ""

            # Source
            if v1 in regions_memoization:
                label1, region1 = regions_memoization[v1]
            else:
                contig1, start1, end1 = get_region(
                    segments, in_links, out_links, ref_offsets, v1
                )
                if contig1 != "":
                    label1 = flag(trees[contig1], start1, end1)
                    region1 = f"{contig1}:{start1}-{end1}"
                else:
                    label1 = "Unknown"
                    region1 = "."
                regions_memoization[v1] = (label1, region1)

            # Sink
            if v2 in regions_memoization:
                label2, region2 = regions_memoization[v2]
            else:
                contig2, start2, end2 = get_region(
                    segments, in_links, out_links, ref_offsets, v2
                )
                if contig2 != "":
                    label2 = flag(trees[contig2], start2, end2)
                    region2 = f"{contig2}:{start2}-{end2}"
                else:
                    label2 = "Unknown"
                    region2 = "."
                regions_memoization[v2] = (label2, region2)

            if region1.split(":")[0] != region2.split(":")[0]:
                label = "TRA"
                region = region1 + ">" + region2
            else:
                if label1 == "Unknown" or label2 == "Unknown":
                    label = "Unknown"
                    region = "."
                else:
                    if label1 != label2:
                        label = "Mixed"
                    else:
                        label = label1
                    region = (
                        region1.split(":")[0]
                        + region1.split(":")[1].split("-")[0]
                        + ":"
                        + region2.split(":")[1].split("-")[1]
                    )
        print(line, region, label, sep=",")


if __name__ == "__main__":
    main()
