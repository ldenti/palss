import sys
import argparse
import re

regex = re.compile("([0-9]+[=XID])")


def parse_cs(cs):
    patterns = {
        "=": r"[ACGTN]+",
        ":": r"[0-9]+",
        "*": r"[acgtn][acgtn]",
        "+": r"[acgtn]+",
        "-": r"[acgtn]+",
        "~": r"[acgtn]{2}[0-9]+[acgtn]{2}",
    }
    results = []
    regex = re.compile(r"(=|:|\*|\+|\-|\~)")
    tokens = regex.split(cs)[1:]
    return [(tokens[i], tokens[i + 1]) for i in range(0, len(tokens) - 1, 2)]


def get_run_length(seq):
    c = seq[0]
    run_len = 1
    while run_len < len(seq) and seq[run_len] == c:
        run_len += 1
    return run_len


def merge_cs(cs):
    new_cs = [cs[0]]
    if new_cs[0][0] == ":":
        new_cs[0] = (":", int(new_cs[0][1]))
    for op, d in cs[1:]:
        if op == new_cs[-1][0] and op == ":":
            new_cs[-1] = (":", int(new_cs[-1][1]) + int(d))
        else:
            new_cs.append((op, d))
    return new_cs


def c_to_str(c):
    return "".join([f"{x}{y}" for x, y in c])


def build_cg(cs):
    cg = []
    for op, d in cs:
        if op == ":":
            cg.append((int(d), "="))
        elif op == "*":
            opl = len(d) // 2
            cg.append((opl, "X"))
        elif op == "-":
            opl = len(d)
            cg.append((opl, "D"))
        elif op == "+":
            opl = len(d)
            cg.append((opl, "I"))
        else:
            assert False
    cg1 = [cg[0]]
    for l, op in cg[1:]:
        if op == cg1[-1][1]:
            cg1[-1] = (cg1[-1][0] + l, op)
        else:
            cg1.append((l, op))
    return cg1


def main():
    parser = argparse.ArgumentParser()

    # parser.add_argument("GFA")
    parser.add_argument("GAF")
    parser.add_argument("-s", type=int, default=3, help="Minimum I/D size")
    parser.add_argument("-r", type=int, default=7, help="Run length")

    args = parser.parse_args()

    # Ss = {}
    # for line in open(gaf_fn):
    #     path = line.split("\t")[5]
    #     for x in re.split("[<>]", path[1:]):
    #         Ss[x] = ""
    # for line in open(gfa_fn):
    #     if not line.startswith("S"):
    #         continue
    #     _, x, seq, *_ = line.strip("\n").split("\t")
    #     if x in Ss:
    #         Ss[x] = seq

    for line in open(args.GAF):
        (
            name,
            ql,
            qs,
            qe,
            strand,
            path,
            tpl,
            ps,
            pe,
            resm,
            cl,
            mapq,
            ascore,
            cg,
            cs,
            rs,
            qseq,
            pseq,
        ) = line.strip("\n").split("\t")

        assert int(qs) == 0

        cs = parse_cs(cs[5:])

        qseq = qseq[5:]
        pseq = pseq[5:]

        new_cs = []
        new_ql = 0
        new_qseq = ""

        qofx, pofx = 0, 0
        drm, irm = 0, 0
        for op, d in cs:
            if op == ":":
                opl = int(d)
                #
                new_cs.append((op, d))
                new_ql += opl
                new_qseq += qseq[qofx : qofx + opl]
                #
                qofx += opl
                pofx += opl
            elif op == "*":
                opl = len(d) // 2
                #
                new_cs.append((op, d))
                new_ql += opl
                new_qseq += qseq[qofx : qofx + opl]
                #
                qofx += opl
                pofx += opl
            elif op == "+":
                opl = len(d)
                assert d == qseq[qofx : qofx + opl]

                keep = True
                if opl <= args.s:
                    run_l = get_run_length(pseq[pofx:])
                    if run_l > args.r:
                        keep = False
                if keep:
                    new_cs.append((op, d))
                    new_qseq += d
                    new_ql += opl
                else:
                    # we do not need to do anything here, we just remove
                    irm += opl
                qofx += opl
            elif op == "-":
                opl = len(d)
                assert d == pseq[pofx : pofx + opl]

                keep = True
                if opl <= args.s:
                    run_l = get_run_length(pseq[pofx:])
                    if run_l > args.r:
                        keep = False

                if keep:
                    new_cs.append((op, d))
                else:
                    drm += opl
                    # fill deletion
                    new_qseq += d
                    new_ql += opl
                    new_cs.append((":", opl))
                pofx += opl

        assert len(new_qseq) == new_ql

        new_cs = merge_cs(new_cs)
        new_cg = build_cg(new_cs)
        new_cs_s = "".join([f"{x}{y}" for x, y in new_cs])
        new_cg_s = "".join([f"{x}{y}" for x, y in new_cg])

        print(
            name,
            drm,
            irm,
            c_to_str(cs),
            cg[5:],
            ql,
            new_ql,
            new_cs_s,
            new_cg_s,
            file=sys.stderr,
        )

        print(
            name,
            new_ql,
            qs,
            new_ql,
            strand,
            path,
            tpl,
            ps,
            pe,
            resm,
            cl,
            mapq,
            ascore,
            f"cg:Z:{new_cg_s}",
            f"cs:Z:{new_cs_s}",
            rs,
            f"qs:Z:{new_qseq}",
            f"ps:Z:{pseq}",
            sep="\t",
        )


if __name__ == "__main__":
    main()
