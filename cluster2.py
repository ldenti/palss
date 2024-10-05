import sys
import networkx as nx


class SFS:
    def __init__(self, qname, s, l, nodes):
        self.qname = qname
        self.s = s
        self.l = l
        self.nodes = nodes

    def __repr__(self):
        return (
            f"{self.qname}:{self.s}-{self.l}\t{'>'.join([str(x) for x in self.nodes])}"
        )


def main():
    gfa_fn = sys.argv[1]
    sfs_fn = sys.argv[2]

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

    print("Populating specifics", file=sys.stderr)
    SS = []
    for line in open(sfs_fn):
        idx, s, l, nodes, strand, rl = line.strip("\n").split("\t")
        n1, n2 = nodes.split(">")
        n1, n2 = int(n1), int(n2)
        if n1 == -1 or n2 == -1:
            continue
        if n2 - n1 > 100:
            # assuming topological order, remove nodes that are too far away
            # we could do a bounded two-way visit
            # print(n1, n2, nx.shortest_path(graph, n1, n2))
            continue
        assert n1 <= n2
        nodes = set([n1])
        if n1 != n2:
            allpaths = nx.all_simple_paths(graph, source=n1, target=n2)
            nodes = set(node for path in allpaths for node in path)
        assert len(nodes) > 0
        SS.append(SFS(idx, int(s), int(l), nodes))

    clusters = [[SS[0]]]
    for s in SS[1:]:
        ii = -1
        for i, c in enumerate(clusters):
            cg = set.union(*[set(s2.nodes) for s2 in c])
            if len(set(s.nodes) & cg) > 0:
                ii = i
                break
        if ii != -1:
            clusters[i].append(s)
        else:
            clusters.append([s])

    for i, c in enumerate(clusters):
        # if i != 1718:
        #     continue
        if len(c) < 2:
            continue
        cgraph = set.union(*[set(s.nodes) for s in c])
        coord = f"chr19:{(min(cgraph)-1)*32}-{(max(cgraph)-1)*32}"

        specifics = {}
        for s in c:
            if s.qname not in specifics:
                specifics[s.qname] = []
            specifics[s.qname].append(s)

        merged = []
        for qn, ss in specifics.items():
            s = min([x.s for x in ss])
            e = max([x.s+x.l for x in ss])
            merged.append(SFS(qn, s, e-s, set()))
        # for qn, s in specifics.items():
        #     print(s)

        print(
            i,
            len(c),
            coord,
            ",".join([f"{s.qname}:{s.s}-{s.s+s.l}" for s in merged]),
            ",".join([str(x) for x in cgraph]),
            sep="\t",
        )


if __name__ == "__main__":
    main()
