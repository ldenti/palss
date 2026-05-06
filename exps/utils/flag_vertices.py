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
    print("fn,n,kind,v,l,supp,qual,region,refspan,refvertices,type")
    for line in open(support_table):
        if line.startswith("fn"):
            continue
        line = line.strip("\n")
        fields = line.split(",")
        idx = int(fields[3])

        contig, start, end = "", -1, -1
        rsv, rev = None, None

        assert idx not in segments

        rsv = bfs_find_first_flagged(idx, segments, in_links, 5000)  # XXX: hardcoded
        if rsv != None:
            rev = bfs_find_first_flagged(
                idx, segments, out_links, 5000
            )  # XXX: hardcoded

        if rsv == None or rev == None:
            print(line, ".", ".", ".", "Unkwnow", sep=",")
        else:
            contig, start, end = "", -1, -1
            for c, offsets in ref_offsets.items():
                if rsv in offsets and rev in offsets:
                    contig = c
                    start = offsets[rsv]
                    end = offsets[rev]
                    break

            assert contig != ""

            t = ""
            if trees[contig].overlaps(start, end):
                t = "Complex"
            else:
                t = "Simple"
            print(
                line, f"{contig}:{start}-{end}", end - start, f"{rsv}>{rev}", t, sep=","
            )


if __name__ == "__main__":
    main()
