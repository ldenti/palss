import sys
import re
from Bio import SeqIO

rc_map = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def split(cigar):
    pairs = []
    p = re.compile(r"([0-9]+[=XID])")  # TODO check if these are enough
    # TODO check why we have empty strings
    for x in p.split(cigar):
        if x == "":
            continue
        pairs.append((int(x[:-1]), x[-1]))
    return pairs


def compact(cigar):
    ccigar = [cigar[0]]
    for l, op in cigar[1:]:
        if op == ccigar[-1][1]:
            ccigar[-1] = (ccigar[-1][0] + l, op)
        else:
            ccigar.append((l, op))
    return ccigar


def cigar2str(cigar):
    return "cg:Z:" + "".join(f"{l}{op}" for l, op in cigar)


def rc(S):
    return "".join([rc_map[c] for c in S[::-1]])


def cut(alignment, interval, gfa_s):
    #   0     1   2   3       4     5     6   7   8         9     10    11   12-
    qname, qlen, qs, qe, strand, path, plen, ps, pe, residues, ablen, mapq, *opt = (
        alignment
    )
    qlen = int(qlen)
    cigar = split((opt[-1].split(":")[2]))

    path_seq = "".join(
        [
            gfa_s[node] if path[0] == ">" else rc(gfa_s[node])
            for node in path[1:].split(path[0])
        ]
    )[
        int(ps) : int(pe) + 1
    ]  # CHECKME pe+1 ?
    pp = 0  # position on path
    qp = int(qs)  # position on query
    ccigar = []
    addflag = False
    stopflag = False
    cut = 0
    pskipped1 = 0
    pkept = 0
    pskipped2 = 0
    for opl, op in cigar:
        if (
            not addflag
            and op != "D"
            and (qp == interval[0] or (qp < interval[0] and qp + opl > interval[0]))
        ):
            addflag = True
            cut = interval[0] - qp
        if addflag and op != "D" and (qp < interval[1] and qp + opl >= interval[1]):
            stopflag = True
            cut = qp + opl - interval[1]

        edit = ""
        if op == "=":
            qp += opl
            pp += opl
        elif op == "X":
            edit = path_seq[pp : pp + opl]
            qp += opl
            pp += opl
        elif op == "D":
            edit = path_seq[pp : pp + opl]
            pp += opl
        elif op == "I":
            qp += opl
        else:
            assert False, f"Unknown operation in CIGAR string {op}"

        if addflag:
            if cut != 0 and op != "I" and not stopflag:
                pskipped1 += cut
            ccigar.append((opl - cut, op, edit))
            if op != "I":
                pkept += opl - cut
            cut = 0
        else:
            if op != "I":
                pskipped1 += opl
        if stopflag:
            break

    # print("Skipped prefix:", pskipped1, file=sys.stderr)
    # print("Path length:", pkept, file=sys.stderr)
    alignment[7] = int(alignment[7]) + pskipped1
    alignment[8] = alignment[7] + pkept
    alignment[2] = interval[0]
    alignment[3] = interval[1]
    alignment[-1] = ccigar

    p = 0
    s, e = alignment[7], alignment[8]
    new_path = []
    new_path_l = 0
    new_s, new_e = 0, -1
    for node in path[1:].split(path[0]):
        l = len(gfa_s[node])
        if p + l - 1 < s:
            pass  # shift += l
        elif p < s and p + l - 1 >= s:
            # print(f"Starting {p} {s} {p+l-1}", file=sys.stderr)
            new_path.append(str(node))
            new_path_l += l
            new_s = s - p
        elif p >= s and p + l - 1 <= e:
            new_path.append(str(node))
            new_path_l += l
        elif p < e and p + l >= e:
            # print(f"Ending {p} {p+l} {e}", file=sys.stderr)
            new_path.append(str(node))
            new_path_l += l
            new_e = new_path_l - (p + l - e)
            break
        else:
            # print(f"Ending {p} {p+l} {e}", file=sys.stderr)
            new_e = new_path_l
            break
        p += l
    # print("New path lenght:", len(new_path), len(new_path)*32)
    sep = path[0]
    alignment[5] = sep + sep.join(new_path)
    alignment[6] = new_path_l
    alignment[7] = new_s
    alignment[8] = new_e


def get_split(alignments):
    intervals = []
    for al in alignments:
        qname, qlen, qs, qe, strand, path, plen, ps, pe, residues, ablen, mapq, *opt = (
            al
        )
        qs, qe = int(qs), int(qe)
        if len(intervals) == 0 or qs > intervals[-1][1]:
            intervals.append((qs, qe))
        else:
            intervals.append((intervals[-1][1], qe))
            intervals[-2] = (intervals[-2][0], qs)
    return intervals


def smooth(alignments, read, ofq, ogaf, minl=10):
    seq = ""
    p = 0
    full_cigar = []
    poffset = 0
    for al in alignments:
        int_s, int_e, cigar = al[2], al[3], al[-1]
        if int_s > p:
            op = "N"
            if p == 0:
                op = "H"  # TODO S/H if hard clipping
            # TODO decide if we want to add clips to the sequence
            if op == "N":
                seq += read[p:int_s]
            # full_cigar.append((int_s - p, op))
            p += int_s - p

        insertions = []
        deletions = []
        rp = 0
        for l, op, _ in cigar:
            if op == "I":
                insertions.append((rp, rp+1))
            elif op == "D":
                deletions.append((rp, rp+l))
            if op != "I":
                rp += l

        itosmooth = []
        if len(insertions) > 0:
            i_clusters = [[insertions[0]]]
            for (s,e) in insertions[1:]:
                if s - i_clusters[-1][-1][1] < 20: # FIXME hardcoded
                    i_clusters[-1].append((s,e))
                else:
                    i_clusters.append([(s,e)])
            for c in i_clusters:
                f = sum([e - s for s,e in c]) < minl
                itosmooth += [f for _ in c]
        dtosmooth = []
        if len(deletions) > 0:
            d_clusters = [[deletions[0]]]
            for (s,e) in deletions[1:]:
                if s - d_clusters[-1][-1][1] < 20: # FIXME hardcoded
                    d_clusters[-1].append((s,e))
                else:
                    d_clusters.append([(s,e)])
            for c in d_clusters:
                f = sum([e - s for s,e in c]) < minl
                dtosmooth += [f for _ in c]

        local_cigar = []
        plen = 0
        iidx, didx = 0, 0
        al[2] += poffset  # move starting position
        print(cigar)
        for l, op, c in cigar:
            assert (op in ["=", "I"] and c == "") or (
                op in ["X", "D"] and c != ""
            ), f"{l} {op} {c}"
            if op == "=":
                print(p, l, op, read[p : p + l])
                seq += read[p : p + l]
                local_cigar.append((l, op))
                p += l
            elif op == "X":
                seq += c
                print(p, l, op, c)
                local_cigar.append((l, "="))
                p += l
            elif op == "I":
                if itosmooth[iidx]:
                    p += l
                    poffset -= l
                    print("Smoothing", l, op)
                else:
                    seq += read[p : p + l]
                    print(p, l, op, read[p : p + l])
                    local_cigar.append((l, op))
                    p += l
                iidx += 1
            elif op == "D":
                if dtosmooth[didx]:
                    seq += c
                    print("Smoothing", l, op, c)
                    local_cigar.append((l, "="))
                    poffset += l
                else:
                    print(p, l, op, c)
                    local_cigar.append((l, op))
                didx += 1
            else:
                assert False, f"Unkown CIGAR operation, {op}, why?"
        al[3] += poffset  # move ending position
        al[-1] = cigar2str(compact(local_cigar))
        full_cigar += local_cigar
    if p < int(alignments[-1][1]):
        # TODO decide if we want to add clips to the sequence
        pass
        # seq+=read[p:]
        # full_cigar.append((len(read) - p, "H")) # TODO S/H
    full_cigar = compact(full_cigar)
    print(f"@{alignments[0][0]} {cigar2str(full_cigar)}", file=ofq)
    print(seq, file=ofq)
    print("+", file=ofq)
    print(seq, file=ofq)

    for al in alignments:
        # qname, qlen, qs, qe, strand, path, plen, ps, pe, residues, ablen, mapq, *opt
        if al[3] - al[2] > 100: # FIXME: hardcoded
            print(
                al[0],
                len(seq),
                al[2],
                al[3],
                al[4],
                al[5],
                al[6],
                al[7],
                al[8],
                "*",
                "*",
                al[11],
                al[-1],
                sep="\t",
                file=ogaf,
            )  # al[9], al[10],
    alignments.append(len(seq))


def main():
    gfa_fn = sys.argv[1]
    gaf_fn = sys.argv[2]
    fq_fn = sys.argv[3]
    oprefix = sys.argv[4]

    out_fq = open(oprefix + ".smooth.fq", "w")
    out_gaf = open(oprefix + ".smooth.gaf", "w")

    # TODO do not load in memory the full FASTQ
    # (this should be easier when reimplementing this into graphaligner)
    reads = {}
    for read in SeqIO.parse(fq_fn, "fastq"):
        reads[read.id] = str(read.seq)

    gfa_s = {}
    for line in open(gfa_fn):
        if line.startswith("S"):
            _, idx, seq, *_ = line.strip("\n").split("\t")
            gfa_s[idx] = seq

    alignments = []
    for line in open(gaf_fn):
        tokens = line.strip("\n").split("\t")
        if int(tokens[11]) < 20:
            # MAPQ
            # print("Skipping", tokens[0], "due to low MAPQ", file=sys.stderr)
            continue

        if len(alignments) > 0 and alignments[-1][0] != tokens[0]:
            # we got a new read, so we analyze the current read
            if len(alignments) > 3:
                # TODO hardcoded
                print(
                    "Skipping",
                    alignments[0][0],
                    "due to multiple alignments",
                    file=sys.stderr,
                )
                alignments = []
            else:
                alignments.sort(key=lambda x: int(x[2]))
                intervals = get_split(alignments)
                for i in range(len(alignments)):
                    cut(alignments[i], intervals[i], gfa_s)
                smooth(alignments, reads[alignments[0][0]], out_fq, out_gaf)
                alignments = []
                # sys.exit(1)
        alignments.append(tokens)
    if len(alignments) <= 3:
        alignments.sort(key=lambda x: int(x[2]))
        intervals = get_split(alignments)
        for i in range(len(alignments)):
            cut(alignments[i], intervals[i], gfa_s)
        smooth(alignments, reads[alignments[0][0]], out_fq, out_gaf)
    else:
        # TODO
        print(
            "Skipping", alignments[0][0], "due to multiple alignments", file=sys.stderr
        )

    out_fq.close()
    out_gaf.close()


if __name__ == "__main__":
    # TODO argparse
    # clips, verbose, minlen, numsplits, mapq
    main()
