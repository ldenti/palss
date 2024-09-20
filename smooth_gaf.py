import sys
import re
from Bio import SeqIO

rc_map = {"A":"T", "C":"G", "G":"C", "T":"A", "N":"N"}

def split(cigar):
    pairs = []
    p = re.compile(r"([0-9]+[=XID])") # TODO check if these are enough
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
    qname, qlen, qs, qe, strand, path, plen, ps, pe, residues, ablen, mapq, *opt = alignment

    qlen = int(qlen)
    cigar = split((opt[-1].split(":")[2]))
    
    path_seq = "".join([gfa_s[node] if path[0] == ">" else rc(gfa_s[node]) for node in path[1:].split(path[0])])[int(ps):int(pe)+1] # CHECKME pe+1 ?
    pp = 0 # position on path
    qp = int(qs) # position on query
    ccigar = []
    for opl, op in cigar:
        if op == "=":
            ccigar.append((opl, op, ""))
            qp += opl
            pp+= opl
        elif op == "X":
            ccigar.append((opl, op, path_seq[pp:pp+opl]))
            qp += opl
            pp+=opl
        elif op == "D":
            ccigar.append((opl,op, path_seq[pp:pp+opl]))
            pp+=opl
        elif op == "I":
            ccigar.append((opl, op, ""))
            qp += opl
        else:
            assert False, f"Unkown operation in CIGAR string {op}"
    alignment.append(interval[0])
    alignment.append(interval[1])
    alignment.append(ccigar)


def get_split(alignments):
    intervals = []
    for al in alignments:
        qname, qlen, qs, qe, strand, path, plen, ps, pe, residues, ablen, mapq, *opt = al
        qs, qe= int(qs), int(qe)
        if len(intervals) == 0 or qs > intervals[-1][1]:
            intervals.append((qs, qe))
        else:
            intervals.append((intervals[-1][1], qe))
            intervals[-2] = (intervals[-2][0], qs)
    return intervals


def smooth(alignments, read, minl=10):
    seq = ""
    p = 0
    full_cigar = []
    for al in alignments:
        int_s, int_e, cigar = al[-3:]
        if int_s > p:
            # TODO decide if we want to add clips to the sequence
            # seq += read[p:int_s]
            op = "N"
            if p == 0:
                op = "S"
            full_cigar.append((int_s - p, op))
            p += int_s - p
        for l, op, c in cigar:
            assert (op in ["=", "I"] and c == "") or (op in ["X", "D"] and c != ""), f"{l} {op} {c}"
            if op == "=":
                seq += read[p:p+l]
                full_cigar.append((l, op))
                p += l
            elif op == "X":
                seq += c
                full_cigar.append((l, "="))
                p += l
            elif op == "I":
                if l < minl:
                    p+=l
                else:
                    seq += seq[p:p+l]
                    full_cigar.append((l, op))
                    p += l
            elif op == "D":
                if l < minl:
                    # print(".", ">", c)
                    seq += c
                    full_cigar.append((l, "="))
                else:
                    full_cigar.append((l, op))
        assert p == int_e
    if p < int(alignments[-1][1]):
        # TODO decide if we want to add clips to the sequence
        # seq+=read[p:]
        full_cigar.append((len(read) - p, "S"))
    full_cigar = compact(full_cigar)
    print(f"@{alignments[0][0]} {cigar2str(full_cigar)}")
    print(seq)
    print("+")
    print(seq)

def main():
    gfa_fn = sys.argv[1]
    gaf_fn = sys.argv[2]
    fq_fn = sys.argv[3]
    
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
        if len(alignments) == 0:
            alignments.append(tokens)
            continue

        if alignments[-1][0] != tokens[0]:
            # we got a new read, so we analyze the current read
            if len(alignments) > 3:
                # TODO hardcoded
                print("Skipping", alignments[0][0], "due to multiple alignments", file=sys.stderr)
                alignments = []
            else:
                alignments.sort(key=lambda x: x[2])
                intervals = get_split(alignments)
                for i in range(len(alignments)):
                    cut(alignments[i], intervals[i], gfa_s)
                smooth(alignments, reads[alignments[0][0]])
                alignments = []
        alignments.append(tokens)
    if len(alignments) <= 3:
        alignments.sort(key=lambda x: x[2])
        intervals = get_split(alignments)
        for i in range(len(alignments)):
            cut(alignments[i], intervals[i], gfa_s)
        smooth(alignments, reads[alignments[0][0]])
    else:
        # TODO
        print("Skipping", alignments[0][0], "due to multiple alignments", file=sys.stderr)
    

if __name__ == "__main__":
    # TODO argparse
    # clips, verbose, minlen, numsplits, mapq
    main()