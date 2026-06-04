import sys
import argparse
import re

# import time

trans = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(trans)[::-1]


"""
Ss: {v : (seq, seq_len, is_known)}
Ls: {(v,strand) : {(v, strand} : is_known}
"""


def load_graph(gfa_fn):
    print("Parsing GFA...", file=sys.stderr)
    Ss, Ls = {}, {}
    kSs = set()
    kLs = {}
    new_paths = {}

    for line in open(gfa_fn):
        if line.startswith("S"):
            _, idx, seq, *_ = line.strip("\n").split("\t")
            idx = int(idx)
            Ss[idx] = (seq, len(seq), False)
        elif line.startswith("L"):
            _, idx1, s1, idx2, s2, *_ = line.strip("\n").split("\t")
            idx1, idx2 = int(idx1), int(idx2)
            # L direction
            k = (idx1, s1 == "+")
            v = (idx2, s2 == "+")
            if k not in Ls:
                Ls[k] = {}
            Ls[k][v] = False
            # add also bidirection
            k = (idx2, s2 != "+")
            v = (idx1, s1 != "+")
            if k not in Ls:
                Ls[k] = {}
            Ls[k][v] = False
        elif line.startswith("P"):
            # XXX: assuming P lines to be new paths
            line = line.strip("\n").split("\t")
            name = line[1]
            new_paths[name] = [(int(x[:-1]), x[-1] == "+") for x in line[2].split(",")]
        elif line.startswith("W"):
            line = line.strip("\n").split("\t")
            # XXX: do we need strand here? If v+ is known and we have v-, what happens?
            path = list(re.findall(r"[<>]\d+|[^<>]+", line[6]))
            path = [(int(x[1:]), x[0] == ">") for x in path]

            Ss[path[0][0]] = (Ss[path[0][0]][0], Ss[path[0][0]][1], True)
            for t1, t2 in zip(path[:-1], path[1:]):
                Ss[t2[0]] = (Ss[t2[0]][0], Ss[t2[0]][1], True)

                assert t1 in Ls
                Ls[t1][t2] = True
                # flag bidirection

                t1 = (t1[0], not t1[1])
                t2 = (t2[0], not t2[1])
                assert t2 in Ls
                Ls[t2][t1] = True

    for k, x in Ss.items():
        if x[2]:
            print("NOVEL", k, file=sys.stderr)
    print("Known vertices:", sum([x[2] for k, x in Ss.items()]), file=sys.stderr)
    print("New paths:", len(new_paths), file=sys.stderr)

    return Ss, Ls, new_paths


def dp(Ss, Ls, P, source, source_strand, sink, sink_strand):
    from collections import defaultdict

    # Define DP[i] = set of vertices v such that there exists a path from source to v that exactly matches S[:i] (so i matches)

    m = len(P)

    DP = [set() for _ in range(m + 1)]
    DP[Ss[source][1]].add((source, source_strand))

    # store predecessors for backtracking
    # predecessors: pred[(v, j)] -> list of (u, i) where i + len(label[v]) == j
    predecessors = defaultdict(list)

    for i in range(m + 1):
        for u, u_strand in DP[i]:
            # print(i, u, u_strand, Ls.get((u, u_strand), ()))
            for v, v_strand in Ls.get((u, u_strand), ()):
                v_seq, v_len, _ = Ss[v]
                if not v_strand:
                    v_seq = revcomp(v_seq)
                j = i + v_len
                if j <= m and P[i:j] == v_seq:
                    DP[j].add((v, v_strand))
                    predecessors[(v, v_strand, j)].append((u, u_strand, i))
                # TODO: store idx in frame j, so we have there vertex and strand

    # we should have sink in DP[m]
    # CHECKME: but in some cases we also have other vertices
    if (sink, sink_strand) not in DP[m]:
        yield []

    # for i, dp in enumerate(DP):
    #     print(f"DP[{i}] ({len(dp)}) [ ", end="")
    #     for x in dp:
    #         print(x[0], end=" ")
    #     print("]")

    for k in predecessors:
        predecessors[k].sort(key=lambda x: x[0], reverse=True)

    # for k, Vs in predecessors.items():
    #     print(f"PRED[({k[0]},{k[2]})] ({len(Vs)}) [ ", end="", file=sys.stderr)
    #     for v in Vs:
    #         print(f"({v[0]},{v[2]})", end=" ", file=sys.stderr)
    #     print("]", file=sys.stderr)

    start_pos = Ss[source][1]
    assert len(DP[start_pos]) == 1

    stack = [(sink, sink_strand, m, [(sink, sink_strand)])]

    while stack:
        # print("STACK:", len(stack), " ".join([f"{x[0]},{x[1]}" for x in stack]))
        # print(len(stack), "paths in the queue")
        vertex, vertex_strand, pos, path = stack.pop()
        if pos == start_pos:
            path.reverse()
            yield path
        else:
            for new_vertex, new_vertex_strand, new_pos in predecessors[
                (vertex, vertex_strand, pos)
            ]:
                # print(vertex, pos, ":", new_vertex, new_pos)
                stack.append(
                    (
                        new_vertex,
                        new_vertex_strand,
                        new_pos,
                        path + [(new_vertex, new_vertex_strand)],
                    )
                )
            # print("")
    # elapsed_time = time.perf_counter() - start_time
    # print(f'BT: {elapsed_time}s', file=sys.stderr)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("GFA")
    parser.add_argument("GAF")
    parser.add_argument("-s", "--supp", default=2, type=int)
    parser.add_argument("-m", "--maxp", default=1024, type=int)
    parser.add_argument("--gaf", default="", type=str)

    args = parser.parse_args()

    print("Parsing GAF...", file=sys.stderr)
    np_reads_support = {}
    for line in open(args.GAF):
        fields = line.strip("\n").split("\t")
        name = fields[0]
        np_reads_support[name] = set(fields[15].split(":")[-1].split("|"))

    Ss, Ls, new_paths = load_graph(args.GFA)

    real_novel = {}
    for cname, path in new_paths.items():
        seq = ""
        for v, s in path:
            seq += Ss[v][0] if s else revcomp(Ss[v][0])

        # if cname != "palss-100.0":
        #     continue

        if len(path) == 1:
            print(f"Skipping path {cname}: single vertex path", file=sys.stderr)
            continue

        # FIXME: actually if we have all known vertices, we might have a deletion. So we also have to keep track of novel edges. Here and later
        if all([Ss[v][2] for v, _ in path]):
            # print("".join([f"{'<>'[s]}{v}" for v, s in path]))
            # print(",".join([f"{v}{'-+'[s]}" for v, s in path]))
            # print([Ls[v1][v2] for v1, v2 in zip(path[:-1], path[1:])])
            # for v1, v2 in zip(path[:-1], path[1:]):
            #     print(v1, v2, v2 in Ls[v1], Ls[v1])
            #     print(f'"L\\t{v1[0]}\\t\\{'-+'[v1[1]]}\\t{v2[0]}\\t\\{'-+'[v2[1]]}"')
            # print(v1, v2, Ss[v1[0]][2], Ss[v2[0]][2], Ls[v1][v2])
            if all([Ls[v1][v2] for v1, v2 in zip(path[:-1], path[1:])]):
                # print(
                #     f"Skipping path {cname}: all known. Seq: {len(seq)}",
                #     file=sys.stderr,
                # )
                continue

        # real_novel[cname] = 1
        # continue
        # assert len(path) > 1

        source, source_strand = path[0]
        sink, sink_strand = path[-1]

        # print(",".join([f"{v}{'-+'[s]}[{'NK'[Ss[v][2]]}]" for v, s in path]))

        best_novel_path_count = float("inf")
        best_novel_path = []
        for n, p in enumerate(
            dp(Ss, Ls, seq, source, source_strand, sink, sink_strand), 1
        ):
            if n - 1 == args.maxp:
                break

            for (v1, s1), (v2, s2) in zip(path[:-1], path[1:]):
                assert (v2, s2) in Ls[
                    (v1, s1)
                ], f"{cname}: L error in path: {v1},{s1}, {v2},{s2}"

            # print(",".join([f"{v}{'-+'[s]}[{'NK'[Ss[v][2]]}]" for v, s in p]))
            for (v1, s1), (v2, s2) in zip(p[:-1], p[1:]):
                assert (v2, s2) in Ls[
                    (v1, s1)
                ], f"{cname}: L error in p: {v1},{s1}, {v2},{s2}"

            rseq = ""
            for v, s in p:
                rseq += Ss[v][0] if s else revcomp(Ss[v][0])
            # print(seq)
            # print(rseq)
            assert seq == rseq, f"{cname} : wrong path from dp"

            # print([x[0] for x in p], file=sys.stderr)
            # print(
            #     cname,
            #     sum([not Ss[v][2] for v, _ in p]),
            #     sum([not Ls[v1][v2] for v1, v2 in zip(p[:-1], p[1:])]),
            #     file=sys.stderr,
            # )
            # for v1, v2 in zip(p[:-1], p[1:]):
            #     print(v1[0], v2[0], not Ls[v1][v2], file=sys.stderr)

            if all([Ss[v][2] for v, _ in p]) and all(
                [Ls[v1][v2] for v1, v2 in zip(p[:-1], p[1:])]
            ):
                best_novel_path = []  # to force skipping this consensus
                # print(
                #     "K",
                #     cname,
                #     len(path),
                #     len(seq),
                #     ",".join([f"{v}{'-+'[s]}[{'NK'[Ss[v][2]]}]" for v, s in p]),
                #     file=sys.stderr,
                # )
                break
            novel_path_count = sum([not Ss[v][2] for v, _ in p])
            if novel_path_count < best_novel_path_count:
                best_novel_path_count = novel_path_count
                best_novel_path = p
            # print(
            #      "N",
            #      cname,
            #     novel_path_count,
            #     ",".join([f"{v}{'-+'[s]}[{'NK'[Ss[v][2]]}]" for v, s in p]),
            #     file=sys.stderr,
            # )
        if best_novel_path != []:
            print(cname, best_novel_path_count, len(best_novel_path), file=sys.stderr)
            real_novel[cname] = best_novel_path
            # print(
            #     "X",
            #     cname,
            #     ",".join([f"{v}{'-+'[s]}[{'NK'[Ss[v][2]]}]" for v, s in path]),
            #     file=sys.stderr,
            # )
            # print(
            #     "X",
            #     cname,
            #     ",".join([f"{v}{'-+'[s]}[{'NK'[Ss[v][2]]}]" for v, s in best_novel_path]),
            #     file=sys.stderr,
            # )
        assert n != 0, f"{cname}"

    print(
        f"{len(real_novel)}/{len(new_paths)} novel paths will be retained.",
        file=sys.stderr,
    )

    print(
        f"Selecting segments/links to remove...",
        file=sys.stderr,
    )

    # XXX: we need 2 passes since if a vertex is supported by at least one path, we want to keep it even if it's not supported by others
    newV = {}
    newL = {}
    for cname, path in real_novel.items():
        print(cname, ">", end="", file=sys.stderr)
        for v, strand in path:
            print("", v, end="", file=sys.stderr)
            if not Ss[v][2]:
                if v not in newV:
                    newV[v] = set()
                newV[v].add(cname)
        print("", file=sys.stderr)
        for (v1, strand1), (v2, strand2) in zip(path[:-1], path[1:]):
            if not Ls[(v1, strand1)][(v2, strand2)]:
                if (v1, strand1) not in newL:
                    newL[(v1, strand1)] = {}
                if (v2, strand2) not in newL[(v1, strand1)]:
                    newL[(v1, strand1)][(v2, strand2)] = set()
                newL[(v1, strand1)][(v2, strand2)].add(cname)

    # print(len(newV), len(newL), file=sys.stderr)

    v_to_remove = set()
    for cname, path in new_paths.items():
        for v, strand in path:
            if not Ss[v][2] and v not in newV:
                v_to_remove.add(v)
    print(
        f"Need to remove {len(v_to_remove)} novel vertices",
        file=sys.stderr,
    )

    l_to_remove = set()
    for cname, path in new_paths.items():
        for (v1, strand1), (v2, strand2) in zip(path[:-1], path[1:]):
            if not Ls[(v1, strand1)][(v2, strand2)]:
                if (v1, strand1) not in newL or (v2, strand2) not in newL[
                    (v1, strand1)
                ]:
                    l_to_remove.add(((v1, strand1), (v2, strand2)))
                    l_to_remove.add(((v2, not strand2), (v1, not strand1)))

    print(
        f"Need to remove {len(l_to_remove)} novel edges",
        file=sys.stderr,
    )

    for v, Cs in newV.items():
        reads = set().union(*[np_reads_support[c] for c in Cs])
        if len(reads) < args.supp:
            v_to_remove.add(v)
    print(
        f"Need to remove {len(v_to_remove)} novel vertices",
        file=sys.stderr,
    )

    for e1, E2 in newL.items():
        for e2, Cs in E2.items():
            reads = set().union(*[np_reads_support[c] for c in Cs])
            if len(reads) < args.supp:
                v1, strand1 = e1
                v2, strand2 = e2

                l_to_remove.add(((v1, strand1), (v2, strand2)))
                l_to_remove.add(((v2, not strand2), (v1, not strand1)))
    print(
        f"Need to remove {len(l_to_remove)} novel edges",
        file=sys.stderr,
    )

    if args.gaf != "":
        ogaf = open(args.gaf, "w")
        for line in open(args.GAF):
            fields = line.strip("\n").split("\t")
            name = fields[0]
            if name in real_novel:
                print(line, end="", file=ogaf)
        ogaf.close()

    # dump_gfa(gfa_fn, np_reads_support, v_to_remove, e_to_remove)

    print("Outputting...", file=sys.stderr)
    for line in open(args.GFA):
        if line.startswith("S"):
            _, v, seq = line.strip("\n").split("\t")
            v = int(v)
            if v in v_to_remove:
                continue
        elif line.startswith("L"):
            _, v1, strand1, v2, strand2, *_ = line.split("\t")
            v1, v2 = int(v1), int(v2)
            strand1, strand2 = strand1 == "+", strand2 == "+"
            if v1 in v_to_remove or v2 in v_to_remove:
                continue
            if ((v1, strand1), (v2, strand2)) in l_to_remove or (
                (v2, not strand2),
                (v1, not strand1),
            ) in l_to_remove:
                continue
            # e = (min(v1, v2), max(v1, v2))
            # if e in e_to_remove:
            #     continue
        elif line.startswith("P"):
            name = line.strip("\n").split("\t")[1]
            if name in new_paths:
                # do not print augmented paths referring to consensus or mapping paths
                continue
        print(line, end="")


if __name__ == "__main__":
    main()
