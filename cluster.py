import sys
import re
import networkx as nx


class SFS:
    def __init__(self, qname, s, e):
        self.qname = qname
        self.s = s  # 0-based [
        self.e = e  # 0-based )
        self.nodes = []
        self.path = False

    def place(self, nodes, path=False):
        self.nodes = nodes
        self.path = path

    def __repr__(self):
        return f"{self.qname}:{self.s}-{self.e}\t{self.path}\t{'>'.join([str(x) for x in self.nodes])}"


def split(cigar):
    pairs = []
    p = re.compile(r"([0-9]+[=XID])")  # TODO check if these are enough
    # TODO check why we have empty strings
    for x in p.split(cigar):
        if x == "":
            continue
        pairs.append((int(x[:-1]), x[-1]))
    return pairs


# Strong assumption: no SFS cover more than 1 smoothed alignment. For this reason, we filter out short alignments while smoothing (we smooth the read though)
def place(specifics, alignments, graph):
    specifics_it = iter(specifics)
    alignments_it = iter(alignments)

    sfs = next(specifics_it, None)
    al = next(alignments_it, None)
    al_i = 0

    updateal = True
    while sfs != None and al != None:
        s, e = sfs.s, sfs.e
        if True:  # updateal:
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
            nidx = 0  # node idx
            cigar = split((opt[-1].split(":")[2]))
            cx = 0  # current cigar operation
            path = [int(x) for x in path[1:].split(path[0])]
            qs, qe = int(qs), int(qe)  # [qs,qe)
            ps, pe = int(ps), int(pe)  # [ps,pe)
            qp = qs  # current position on query
            updateal = False
        if e <= qs:
            # Case 1: sfs precedes alignment
            if al_i == 0:
                # TODO Clipping at begin
                pass
            else:
                if alignments[al_i - 1][5][0] == alignments[al_i][5][0]:
                    x, y = int(
                        alignments[al_i - 1][5][1:].split(alignments[al_i - 1][5][0])[
                            -1
                        ]
                    ), int(alignments[al_i][5][1:].split(alignments[al_i][5][0])[0])
                    if alignments[al_i - 1][5][0] == "<":
                        x, y = y, x
                    if x > y:
                        # CHECKME: we may just need to invert x and y
                        print(
                            f"(1) Skipping subgraph from {x} to {y} - insertion on {alignments[al_i][0]} {alignments[al_i][5][0]} ?",
                            file=sys.stderr,
                        )
                    else:
                        if (
                            y - x < 100
                        ):  # hardcoded + approximation assuming topological sorting
                            print(
                                f"Computing subgraph from {x} to {y}", file=sys.stderr
                            )
                            allpaths = nx.all_simple_paths(graph, source=x, target=y)
                            nodes = set(node for path in allpaths for node in path)
                            sfs.place(nodes)
                        else:
                            print(f"Skipping subgraph from {x} to {y}", file=sys.stderr)
                else:
                    # TODO complex event?
                    pass
            sfs = next(specifics_it, None)
        elif s < qs and e > qs:
            # Case 2: sfs starts before alignment and ends after
            if al_i == 0:
                # TODO Clipping at begin
                pass
            else:
                if alignments[al_i - 1][5][0] == alignments[al_i][5][0]:
                    x = int(
                        alignments[al_i - 1][5][1:].split(alignments[al_i - 1][5][0])[
                            -1
                        ]
                    )

                    subpath = []
                    used = ps
                    while True:
                        node = path[nidx]
                        node_l = graph.nodes[node]["l"] - used
                        l, op = cigar[cx]
                        if op == "I":
                            cigar[cx] = (0, cigar[cx][1])
                            cx += 1
                        else:
                            if l >= node_l:
                                nidx += 1
                                used = 0
                                l = node_l
                                cigar[cx] = (cigar[cx][0] - node_l, cigar[cx][1])
                                if cigar[cx][0] <= 0:
                                    cx += 1
                            else:
                                used += l
                                cigar[cx] = (0, cigar[cx][1])
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
                    y = subpath[0]
                    if alignments[al_i - 1][5][0] == "<":
                        x, y = y, x
                    if x > y:
                        # CHECKME: we may just need to invert x and y
                        print(
                            f"(2) Skipping subgraph from {x} to {y} - insertion on {alignments[al_i][0]} {alignments[al_i][5][0]} ?",
                            file=sys.stderr,
                        )
                    else:
                        if (
                            y - x < 100
                        ):  # hardcoded + approximation assuming topological sorting
                            print(
                                f"Computing subgraph from {x} to {y}", file=sys.stderr
                            )
                            allpaths = nx.all_simple_paths(graph, source=x, target=y)
                            nodes = set(node for path in allpaths for node in path)
                            sfs.place(nodes | set(subpath))
                        else:
                            print(f"Skipping subgraph from {x} to {y}", file=sys.stderr)
                else:
                    # TODO complex event?
                    pass
            sfs = next(specifics_it, None)
        elif s >= qs and e <= qe:
            # Case 3: sfs is inside alignment
            subpath = []
            used = ps
            while True:
                node = path[nidx]
                node_l = graph.nodes[node]["l"] - used
                l, op = cigar[cx]
                if op == "I":
                    cigar[cx] = (0, cigar[cx][1])
                    cx += 1
                else:
                    if l >= node_l:
                        nidx += 1
                        used = 0
                        l = node_l
                        cigar[cx] = (cigar[cx][0] - node_l, cigar[cx][1])
                        if cigar[cx][0] <= 0:
                            cx += 1
                    else:
                        used += l
                        cigar[cx] = (0, cigar[cx][1])
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
            sfs.place(sorted(subpath), path=True)
            sfs = next(specifics_it, None)
        elif s < qe and e > qe:
            # Case 3: sfs overlaps alignment end
            if al_i == len(alignments) - 1:
                # TODO Clipping at end
                pass
            else:
                y = int(
                    alignments[al_i + 1][5][1:].split(alignments[al_i + 1][5][0])[0]
                )
                if alignments[al_i][5][0] == alignments[al_i + 1][5][0]:
                    subpath = []
                    used = ps
                    while True:
                        node = path[nidx]
                        node_l = graph.nodes[node]["l"] - used
                        l, op = cigar[cx]
                        if op == "I":
                            cigar[cx] = (0, cigar[cx][1])
                            cx += 1
                        else:
                            if l >= node_l:
                                nidx += 1
                                used = 0
                                l = node_l
                                cigar[cx] = (cigar[cx][0] - node_l, cigar[cx][1])
                                if cigar[cx][0] <= 0:
                                    cx += 1
                            else:
                                used += l
                                cigar[cx] = (0, cigar[cx][1])
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
                                # in this case, subpath is path[nidx:] for sure
                        if cx == len(cigar) or qp >= e:
                            break

                    assert len(subpath) > 0
                    assert subpath[-1] == path[-1]
                    x = subpath[0]
                    if alignments[al_i][5][0] == "<":
                        x, y = y, x

                    if x > y:
                        print(
                            f"(3) Skipping subgraph from {x} to {y} - insertion on {alignments[al_i][0]} {alignments[al_i][5][0]} ?",
                            file=sys.stderr,
                        )
                    else:
                        if (
                            y - x < 100
                        ):  # hardcoded + approximation assuming topological sorting
                            print(
                                f"Computing subgraph from {x} to {y}", file=sys.stderr
                            )
                            allpaths = nx.all_simple_paths(graph, source=x, target=y)
                            nodes = set(node for path in allpaths for node in path)
                            sfs.place(nodes | set(subpath))
                        else:
                            print(
                                f"(3) Skipping subgraph from {x} to {y}",
                                file=sys.stderr,
                            )
                else:
                    # TODO complex event?
                    pass
            al = next(alignments_it, None)
            al_i += 1
            updateal = True
            sfs = next(specifics_it, None)
        elif s >= qe:
            # Case 4: sfs is after current alignment
            # move to next alignment
            al = next(alignments_it, None)
            al_i += 1
            updateal = True


def main():
    gfa_fn = sys.argv[1]
    sfs_fn = sys.argv[2]
    gaf_fn = sys.argv[3]

    print("Iterating over nodes", file=sys.stderr)
    graph = nx.DiGraph()
    for line in open(gfa_fn):
        if line.startswith("S"):
            _, idx, seq, *_ = line.strip("\n").split("\t")
            graph.add_node(int(idx), l=len(seq))
    print("Iterating over edges", file=sys.stderr)
    for line in open(gfa_fn):
        if line.startswith("L"):
            _, idx1, _, idx2, _, *_ = line.strip("\n").split("\t")
            graph.add_edge(int(idx1), int(idx2))

    print("Loading specific strings", file=sys.stderr)
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
            orientations = set(al[5][0] for al in alignments)
            if len(orientations) > 1:
                print(
                    "Alignments with different orientations. Skipping for now",
                    file=sys.stderr,
                )
            else:
                orientation = list(orientations)[0]
                # alignments.sort(key=lambda x: int(x[2]), reverse=orientation == "<")
                alignments.sort(key=lambda x: int(x[2]))
                place(sfs[alignments[0][0]], alignments, graph)
                for _ in sfs[alignments[0][0]]:
                    print(_)

            alignments = []
        alignments.append(tokens)
    if len(alignments) > 0:
        assert (
            len(alignments) <= 3
        ), "Too many alignments for read {alignments[0][0]}, why?"
        orientations = set(al[5][0] for al in alignments)
        if len(orientations) > 1:
            print(
                "Alignments with different orientations. Skipping for now",
                file=sys.stderr,
            )
        else:
            orientation = list(orientations)[0]
            # alignments.sort(key=lambda x: int(x[2]), reverse=orientation == "<")
            alignments.sort(key=lambda x: int(x[2]))
            place(sfs[alignments[0][0]], alignments, graph)

    placed = 0
    total = 0
    for qname, specifics in sfs.items():
        for s in specifics:
            total += 1
            placed += len(s.nodes) > 0
    print(placed, total, placed / total)


if __name__ == "__main__":
    main()
