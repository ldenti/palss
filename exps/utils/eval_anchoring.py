import sys
import os
import glob
import re

regex = re.compile("([0-9]+[=XID])")


def parse_gaf(gaf_fn):
    alignments = {}
    for line in open(gaf_fn):
        qname, ql, qs, qe, strand, path, pl, ps, pe, *_rest, cigar = line.strip(
            "\n"
        ).split("\t")

        ps = int(ps)
        qs, qe, ql = int(qs), int(qe), int(ql)
        c = (qe - qs) / ql

        old_c = 0
        if qname in alignments:
            old_c = alignments[qname][-1]

        if c > old_c:
            cigar = cigar[5:]
            alignments[qname] = [path, cigar, qs, qe, ql, ps, c]
    return alignments


def parse_cigar(cigar):
    tokens = regex.split(cigar)
    tokens = [(int(x[:-1]), x[-1]) for x in tokens if x != ""]
    return tokens


# Get alignment pairs {position : (vertex, operation)}
def process_alignment(segments, alignment):
    path, cigar, qs, qe, ql, ps, c = alignment
    path = re.split("[<>]", path[1:])
    ext_path = []
    for v in path:
        ext_path += [v] * segments[v]

    cigar = parse_cigar(cigar)
    ext_cigar = []
    for opl, op in cigar:
        ext_cigar += [op] * opl

    alignment = {}
    # Fill clips
    for apos in range(0, qs):
        alignment[apos] = ("~", "~")
    for apos in range(qe, ql + 1):
        alignment[apos] = ("~", "~")

    ppos = ps
    apos = qs
    cpos = 0
    while cpos < len(ext_cigar):
        op = ext_cigar[cpos]
        if op == "=" or op == "X":
            assert ppos < len(ext_path), f"{len(ext_path)}, {ppos}, {qname}"
            alignment[apos] = (ext_path[ppos], op)
            apos += 1
            ppos += 1
        elif op == "I":
            alignment[apos] = (ext_path[ppos], op)
            apos += 1
        elif op == "D":
            ppos += 1
        else:
            print("!!!", op, file=sys.stderr)
            return 1
        cpos += 1
    return alignment


def main():
    gfa_fn = sys.argv[1]
    gaf_fn = sys.argv[2]
    txt_fn = sys.argv[3]

    alignments = parse_gaf(gaf_fn)

    segments = {}
    for line in open(gfa_fn):
        if not line.startswith("S"):
            continue
        _, s, seq, *_ = line.strip("\n").split("\t")
        segments[s] = len(seq)

    last_qname = ""
    alignment = None
    for line in open(txt_fn):
        if not line.startswith("0"):
            continue
        _, qname, s, l, e, _, _, _, _, v1, v2, _, _, _, _ = line.strip("\n").split("\t")

        s, e = int(s), int(e)
        if qname != last_qname:
            alignment = process_alignment(segments, alignments[qname])

        cat = ""
        if s not in alignment or e - 1 not in alignment:
            print(
                qname,
                s,
                s in alignment,
                e - 1,
                e - 1 in alignment,
                file=sys.stderr,
            )
            # return 1
        true_vertices = [alignment[s][0], alignment[e - 1][0]]
        if set([v1, v2]) == set(true_vertices):
            cat = "GOOD"
        else:
            if "~" in true_vertices:
                # CLIP
                if len(set(true_vertices)) == 1:
                    cat = "FULL_CLIP"
                elif len(set([v1, v2]) & set(true_vertices)) == 1:
                    cat = "HALF_CLIP"
                else:
                    cat = "BAD_CLIP"
            else:
                if len(set([v1, v2]) & set(true_vertices)) == 1:
                    cat = "HALF_BAD"
                else:
                    cat = "BAD"
        print(
            qname,
            s,
            e,
            v1,
            v2,
            true_vertices[0] if true_vertices[0] != "~" else -1,
            true_vertices[1] if true_vertices[1] != "~" else -1,
            cat,
            sep=",",
        )

        last_qname = qname


if __name__ == "__main__":
    main()
