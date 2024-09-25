import sys
import re


class SFS:
    def __init__(self, qname, s, e):
        self.qname = qname
        self.s = s  # 0-based [
        self.e = e  # 0-based )
        self.nodes = []

    def place(self, nodes, aprx=False):
        self.nodes = nodes
        self.approximate = aprx

    def __repr__(self):
        return (
            f"{self.qname}:{self.s}-{self.e}\t{'>'.join([str(x) for x in self.nodes])}"
        )


def split(cigar):
    pairs = []
    p = re.compile(r"([0-9]+[=XID])")  # TODO check if these are enough
    # TODO check why we have empty strings
    for x in p.split(cigar):
        if x == "":
            continue
        pairs.append((int(x[:-1]), x[-1]))
    return pairs


def place(specifics, alignments, gfa_s):
    specifics_it = iter(specifics)
    alignments_it = iter(alignments)

    sfs = next(specifics_it, None)
    al = next(alignments_it, None)
    al_i = 0

    # sfs_cat = []
    updateal = True
    while sfs != None and al != None:
        s, e = sfs.s, sfs.e
        if updateal:
            (
                qname,
                qlen,
                qs,
                qe,
                strand,
                path,
                plen,
                ps,
                pe,
                residues,
                ablen,
                mapq,
                *opt,
            ) = al
            cigar = split((opt[-1].split(":")[2]))
            path = [int(x) for x in path[1:].split(path[0])]
            qs, qe = int(qs), int(qe)  # [qs,qe)
            ps, pe = int(ps), int(pe)  # [qs,qe)
            qp = qs  # current position on query
            pp = ps  # current position on query
            updateal = False
        if e <= qs:
            if al_i == 0:
                # TODO Clipping at begin
                sfs.place(["SCLIP"], aprx=True)
            else:
                if alignments[al_i - 1][5][0] == alignments[al_i][5][0]:
                    x, y = int(
                        alignments[al_i - 1][5][1:].split(alignments[al_i - 1][5][0])[
                            -1
                        ]
                    ), int(alignments[al_i][5][1:].split(alignments[al_i][5][0])[0])
                    sfs.place([x, "...", y], aprx=True)
                else:
                    # TODO complex event?
                    sfs.place(["COMPLEX"], aprx=True)
            sfs = next(specifics_it, None)
        elif s < qs and e > qs:
            if al_i == 0:
                # TODO Clipping at begin
                sfs.place(["SCLIP"], aprx=True)
            else:
                if alignments[al_i - 1][5][0] == alignments[al_i][5][0]:
                    x = int(
                        alignments[al_i - 1][5][1:].split(alignments[al_i - 1][5][0])[
                            -1
                        ]
                    )

                    subpath = []
                    nx = 0  # node idx
                    cx = 0
                    while True:
                        node = path[nx]
                        node_l = gfa_s[node]
                        l, op = cigar[cx]
                        if l >= node_l:
                            nx += 1
                            l = node_l
                            cigar[cx] = (cigar[cx][0] - node_l, cigar[cx][1])
                            if cigar[cx][0] <= 0:
                                cx += 1
                        else:
                            cx += 1
                        if op == "=" or op == "I":
                            qp += l
                        elif op == "D":
                            pass
                        else:
                            assert False
                        if qp >= s:
                            if len(subpath) == 0 or node != subpath[-1]:
                                subpath.append(node)
                        if qp >= e:
                            break
                    assert len(subpath) > 0
                    sfs.place([x, "..."] + subpath, aprx=True)
                else:
                    # TODO complex event?
                    sfs.place(["COMPLEX"], aprx=True)
            sfs = next(specifics_it, None)
        elif s >= qs and e <= qe:
            subpath = []
            nx = 0  # node idx
            cx = 0
            while True:
                node = path[nx]
                node_l = gfa_s[node]
                l, op = cigar[cx]
                if l >= node_l:
                    nx += 1
                    l = node_l
                    cigar[cx] = (cigar[cx][0] - node_l, cigar[cx][1])
                    if cigar[cx][0] <= 0:
                        cx += 1
                else:
                    cx += 1
                if op == "=" or op == "I":
                    qp += l
                elif op == "D":
                    pass
                else:
                    assert False
                if qp >= s:
                    if len(subpath) == 0 or node != subpath[-1]:
                        subpath.append(node)
                if qp >= e:
                    break
            assert len(subpath) > 0
            sfs.place(sorted(subpath))
            sfs = next(specifics_it, None)
        elif s < qe and e > qe:
            if al_i == len(alignments) - 1:
                # TODO Clipping at end
                sfs.place(["ECLIP"], aprx=True)
            else:
                y = int(
                    alignments[al_i + 1][5][1:].split(alignments[al_i + 1][5][0])[0]
                )
                if alignments[al_i][5][0] == alignments[al_i + 1][5][0]:
                    subpath = []
                    nx = 0  # node idx
                    cx = 0
                    while True:
                        node = path[nx]
                        node_l = gfa_s[node]
                        l, op = cigar[cx]
                        if l >= node_l:
                            nx += 1
                            l = node_l
                            cigar[cx] = (cigar[cx][0] - node_l, cigar[cx][1])
                            if cigar[cx][0] <= 0:
                                cx += 1
                        else:
                            cx += 1
                        if op == "=" or op == "I":
                            qp += l
                        elif op == "D":
                            pass
                        else:
                            assert False
                        if qp >= s:
                            if len(subpath) == 0 or node != subpath[-1]:
                                subpath.append(node)
                        if cx == len(cigar) or qp >= e:
                            break
                    assert len(subpath) > 0
                    sfs.place(subpath + ["...", y], aprx=True)
                else:
                    # TODO complex event?
                    sfs.place(["COMPLEX"], aprx=True)
            sfs = next(specifics_it, None)
        elif s >= qe:
            # move to next alignment
            al = next(alignments_it, None)
            al_i += 1
            updateal = True
    # for x in sfs_cat:
    #     print(alignments[0][0], x[0], x[1])
    # print("")
    #     # [qs,qe)

    #     cigar = split((opt[-1].split(":")[2]))


def main():
    gfa_fn = sys.argv[1]
    sfs_fn = sys.argv[2]
    gaf_fn = sys.argv[3]

    gfa_s = {}
    for line in open(gfa_fn):
        if line.startswith("S"):
            _, idx, seq, *_ = line.strip("\n").split("\t")
            gfa_s[int(idx)] = len(seq)

    sfs = {}
    for line in open(sfs_fn):
        idx, st, le = line.strip("\n").split("\t")
        if idx not in sfs:
            sfs[idx] = []
        sfs[idx].append(SFS(idx, int(st), int(st) + int(le)))

    print(len(sfs))

    alignments = []
    for line in open(gaf_fn):
        tokens = line.strip("\n").split("\t")
        if tokens[0] not in sfs:
            continue
        if int(tokens[11]) < 20:
            continue

        if len(alignments) > 0 and alignments[-1][0] != tokens[0]:
            # we got a new read, so we analyze the current read
            assert (
                len(alignments) <= 3
            ), "Too many alignments for read {alignments[0][0]}, why?"
            place(sfs[alignments[0][0]], alignments, gfa_s)
            for _ in sfs[alignments[0][0]]:
                print(_)

            alignments = []
        alignments.append(tokens)
    if len(alignments) > 0:
        assert (
            len(alignments) <= 3
        ), "Too many alignments for read {alignments[0][0]}, why?"

    placed = 0
    total = 0
    for qname, specifics in sfs.items():
        for s in specifics:
            total += 1
            placed += len(s.nodes) > 0
    print(placed, total, placed / total)


if __name__ == "__main__":
    main()
