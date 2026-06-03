import sys
import argparse
import re
from collections import namedtuple
from Bio.Align import PairwiseAligner

Alignment = namedtuple(
    "Alignment",
    ["idx", "ql", "qs", "qe", "pl", "ps", "pe", "cov", "NM", "path", "cigar"],
)

regex = re.compile("([0-9]+[=XID])")

complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N":"N"}


def rc(sequence):
    return "".join(complement[base] for base in reversed(sequence))


def parse_cigar(cigar):
    tokens = regex.split(cigar)
    tokens = [(int(x[:-1]), x[-1]) for x in tokens if x != ""]
    return tokens


def split_path(path):
    split = re.split(r"([<>])", path)
    split = split[1:]
    return [split[i] + split[i + 1] for i in range(0, len(split), 2)]


def load_segments(gfa_fn, segments):
    for line in open(gfa_fn):
        if line.startswith("S"):
            _, idx, seq, *_ = line.strip("\n").split("\t")
            if idx in segments:
                segments[idx] = seq


def parse_gaf(gaf_fn, min_cov):
    alns = {}
    for line in open(gaf_fn):
        fields = line.strip("\n").split("\t")
        qidx = fields[0]
        ql, qs, qe = (int(x) for x in fields[1:4])
        # [qs, qe)
        pl, ps, pe = (int(x) for x in fields[6:9])

        c = (qe - qs) / ql
        if c < min_cov:
            continue

        nm = fields[13].split(":")
        assert nm[0] == "NM"
        nm = int(nm[2])

        cg = fields[18].split(":")
        assert cg[0] == "cg"
        cg = cg[2]

        aln = Alignment(
            idx=qidx,
            ql=ql,
            qs=qs,
            qe=qe,
            cov=c,
            NM=nm,
            path=fields[5],
            cigar=cg,
            ps=ps,
            pe=pe,
            pl=pl,
        )

        if qidx not in alns:
            alns[qidx] = aln
        else:
            print(f"[W] Contig {qidx} has more alignments", file=sys.stderr)
            old_nm = alns[qidx].NM
            if old_nm >= nm:
                alns[qidx] = aln
    return alns


def generate_cigar_string(alignment):
    matches = 0
    nm = 0
    ext_cigar = []
    seq1, seq2 = alignment[0], alignment[1]
    for a, b in zip(seq1, seq2):
        if a == b:
            ext_cigar.append(("=", a, b))
            matches += 1
        else:
            if a != "-" and b != "-":
                ext_cigar.append(("X", a, b))
                nm += 1
            else:
                op = "D" if b == "-" else "I"
                ext_cigar.append((op, a, b))

    if ext_cigar[0][0] != "=" or ext_cigar[-1][0] != "=":
        return "", "", -1, -1, -1

    total_cigar_length = len(ext_cigar)

    cigar = []
    last_op = ext_cigar[0][0]
    opl = 1
    for op, a, b in ext_cigar[1:]:
        if op != last_op:
            cigar.append(f"{opl}{last_op}")
            last_op = op
            opl = 1
        else:
            opl += 1
    cigar.append(f"{opl}{last_op}")

    cs = []
    p1, p2 = 0, 0
    for o in cigar:
        opl, op = int(o[:-1]), o[-1]
        if op == "=":
            cs.append(f":{opl}")
            p1 += opl
            p2 += opl
        elif op == "X":
            for i in range(opl):
                cs.append(f"*{seq1[p1]}{seq2[p2]}")
                p1 += 1
                p2 += 1
        elif op == "I":
            cs.append("+" + seq2[p2 : p2 + opl])
            p1 += opl  # we need this since in seq1/seq2 we have -
            p2 += opl
        elif op == "D":
            cs.append("-" + seq1[p1 : p1 + opl])
            p1 += opl
            p2 += opl  # we need this since in seq1/seq2 we have -

    return "".join(cigar), "".join(cs), matches, total_cigar_length, nm


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("GFA")
    parser.add_argument("FA")
    parser.add_argument("GAF")
    parser.add_argument("SFS")
    parser.add_argument("-c", "--cov", type=float, default=0.95)
    parser.add_argument("-m", "--matches", type=int, default=7)
    args = parser.parse_args()

    sequences = {}
    idx, seq = "", ""
    for line in open(args.FA):
        if line.startswith(">"):
            if idx != "":
                sequences[idx] = seq
            idx = line[1:-1]
            seq = ""
        else:
            seq += line.strip("\n")
    sequences[idx] = seq

    # ===

    alns = parse_gaf(args.GAF, args.cov)
    alns = {k: v for k, v in alns.items() if v.NM > 0}

    # ===

    segments = {}
    for aln in alns.values():
        for v in split_path(aln.path):
            v = v[1:]
            segments[v] = ""
    load_segments(args.GFA, segments)

    ### Extract alignment pairs
    ext_alns = {}
    for c, aln in alns.items():
        ext_path = []
        path = split_path(aln.path)
        path_seq = ""
        for v in path:
            reverse = v[0] == "<"
            v = v[1:]
            ext_path += [(v, reverse)] * len(segments[v])
            path_seq += segments[v] if not reverse else rc(segments[v])
        assert len(ext_path) == aln.pl
        # ext_path = ext_path[aln.ps :]

        pp = aln.ps
        qp = aln.qs

        cigar = parse_cigar(aln.cigar)
        ext_cigar = []
        for opl, op in cigar:
            ext_cigar += [op] * opl

        assert sum([x != "I" for x in ext_cigar]) == aln.pe - aln.ps

        apairs = []
        # Fill clips
        for _ in range(0, aln.qs):
            apairs.append("~")

        last_v = ext_path[pp][0]
        offset = 0
        for op in ext_cigar:
            if pp - offset == len(segments[last_v]):
                offset += len(segments[last_v])
                last_v = ext_path[pp][0]

            if op in ["X", "="]:
                apairs.append(
                    (
                        ext_path[pp],
                        pp,
                        qp,
                        pp - offset,
                        path_seq[pp],
                        sequences[c][qp],
                        op,
                    )
                )
                pp += 1
                qp += 1
            elif op == "D":
                apairs.append((ext_path[pp], pp, qp, pp - offset, path_seq[pp], "", op))
                pp += 1
            elif op == "I":
                apairs.append(
                    (ext_path[pp], pp, qp, pp - offset, "", sequences[c][qp], op)
                )
                qp += 1
            else:
                print(f"[E] Unknown cigar ({op}) operation for {c}", file=sys.stderr)
                return 1
        for _ in range(aln.qe, aln.ql + 1):
            apairs.append("~")
        ext_alns[c] = (apairs, path_seq)

    for line in open(args.SFS):
        if not line.startswith("0"):
            continue
        fields = line.strip("\n").split("\t")
        qname = fields[1]
        qs, ql, qe = (int(x) for x in fields[2:5])
        seq = fields[-1]
        if qname not in alns:
            continue
        aln = alns[qname]
        if qs >= aln.qs and qe < aln.qe:
            apairs, pseq = ext_alns[qname]
            # print(qs, qe, aln.qs, aln.qe)

            n = 0
            while qs >= 0 and n < args.matches:
                if apairs[qs][-1] != "=":
                    n = 0
                else:
                    n += 1
                qs -= 1
            if qs < aln.qs:
                continue

            n = 0
            while qe < len(apairs) and n < args.matches:
                if apairs[qe][-1] != "=":
                    n = 0
                else:
                    n += 1
                qe += 1
            if qe >= aln.qe:
                continue

            # print(qs, qe, aln.qs, aln.qe)
            # for i in range(qs, qe + 1):
            #     print(apairs[i])
            # print(apairs[qs], apairs[qe])
            p_subseq = pseq[apairs[qs][1] : apairs[qe][1] + 1]
            q_subseq = sequences[qname][apairs[qs][2] : apairs[qe][2] + 1]

            if p_subseq == q_subseq:
                continue

            aligner = PairwiseAligner()
            aligner.open_gap_score = -16
            aligner.extend_gap_score = -2
            aligner.match_score = 1
            aligner.mismatch_score = -9

            alignments = aligner.align(p_subseq, q_subseq)
            cigar, cs, matches, total_cigar_length, nm = generate_cigar_string(
                alignments[0]
            )

            if len(cigar) == 0:
                continue

            path = []
            for q in range(qs, qe + 1):
                v, reverse = apairs[q][0]
                if len(path) == 0 or path[-1] != v:
                    path.append("><"[reverse])
                    path.append(v)

            # print(apairs[qs])
            pl = sum([len(segments[v]) for v in path[1::2]])
            ps = apairs[qs][3]
            pe = ps + len(p_subseq)

            print(
                f"{qname}:{qs}-{qe}",
                apairs[qe][2] - apairs[qs][2] + 1,
                0,
                apairs[qe][2] - apairs[qs][2] + 1,
                "+",
                "".join(path),
                pl,
                ps,
                pe,
                matches,
                total_cigar_length,
                60,
                f"NM:i:{nm}",
                f"cg:Z:{cigar}",
                f"cs:Z:{cs}",
                f"ps:Z:{p_subseq}",
                f"qs:Z:{q_subseq}",
                sep="\t",
            )
        else:
            print(
                f"[W] Skipping {qname}:{qs}-{qe} due to alignment clipping",
                file=sys.stderr,
            )


def gaf_heatmap():
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns

    parser = argparse.ArgumentParser()
    parser.add_argument("GAF")
    # parser.add_argument("SFS")
    # parser.add_argument("-p", "--palss", action="store_true")
    args = parser.parse_args()

    contigs = {}
    intervals = {}
    for line in open(args.GAF):
        fields = line.strip("\n").split("\t")
        qname = fields[0]
        ql, qs, qe = (int(x) for x in fields[1:4])
        if qname not in contigs:
            contigs[qname] = [0] * ql
            intervals[qname] = []
        intervals[qname].append((qs, qe, (qe - qs + 1) / ql))
        for q in range(qs, qe):
            contigs[qname][q] += 1

    for c in intervals:
        intervals[c].sort()
        print(c, ", ".join([f"({a},{b},{x})" for a, b, x in intervals[c]]))

    max_ql = max([len(x) for x in contigs.values()])
    for c in contigs:
        contigs[c] += [-1] * (max_ql - len(contigs[c]))

    data = np.array(list(contigs.values()))
    sns.heatmap(
        data,
        annot=False,
        # fmt="d",
        cmap="coolwarm",
        square=False,
        cbar=True,
        # linewidths=0.5,
        # linecolor="black",
        # alpha=0.5,
        xticklabels=False,
        yticklabels=False,
    )
    plt.title(f"n: {len(contigs)}, l: {max_ql}")
    # plt.show()
    plt.savefig("1.png")


if __name__ == "__main__":
    main()
