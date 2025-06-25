import sys
import re
import random


def main():
    gfa_fn = sys.argv[1]
    gaf_fn = sys.argv[2]
    min_l = 1000

    V = set()
    print("Parsing paths...", file=sys.stderr)
    for line in open(gfa_fn):
        path = []
        if line.startswith("P"):
            line = line.strip("\n").split("\t")
            path = [int(x[:-1]) for x in line[2].split(",")]
        elif line.startswith("W"):
            line = line.strip("\n").split("\t")
            path = [int(x) for x in re.split("[<>]", line[6][1:])]
        else:
            continue
        for v in path:
            V.add(v)

    nV = {}
    Vseqs = {}
    INs = {}
    OUTs = {}
    refs = set()

    print("Reiterating over graph...", file=sys.stderr)
    for line in open(gfa_fn):
        if line.startswith("S"):
            _, v, seq, *_ = line.strip("\n").split("\t")
            v = int(v)
            Vseqs[v] = seq
            if v in V:
                continue
            nV[v] = 0
        elif line.startswith("L"):
            _, idx1, _, idx2, _, _ = line.strip("\n").split("\t")
            idx1, idx2 = int(idx1), int(idx2)
            if idx1 not in OUTs:
                OUTs[idx1] = []
            OUTs[idx1].append(idx2)
            if idx2 not in INs:
                INs[idx2] = []
            INs[idx2].append(idx1)
        elif line.startswith("P"):
            # assuming reference on P line
            _, _, path, _ = line.strip("\n").split("\t")
            path = [int(x[:-1]) for x in path.split(",")]
            for idx in path:
                refs.add(idx)

    print(f"We have {len(nV)} novel vertices", file=sys.stderr)

    print("Iterating over GAF...", file=sys.stderr)
    for line in open(gaf_fn):
        line = line.strip("\n").split("\t")

        path = line[5]
        strand = path[0]
        path = [int(x) for x in re.split("[<>]", path[1:])]

        for v in path:
            if v in nV:
                nV[v] = nV[v] + 1
    for v, w in nV.items():
        precontig = [v]
        l = 0
        while l < min_l:
            if precontig[-1] not in INs:
                break
            Ps = INs[precontig[-1]]
            pred = None
            for p in Ps:
                if p in refs:
                    pred = p
            if pred == None:
                pred = random.choice(Ps)
            precontig.append(pred)
            l += len(Vseqs[pred])
        precontig.reverse()

        postcontig = [v]
        l = 0
        while l < min_l:
            if postcontig[-1] not in OUTs:
                break
            Ss = OUTs[postcontig[-1]]
            succ = None
            for s in Ss:
                if s in refs:
                    succ = s
            if succ == None:
                succ = random.choice(Ss)
            postcontig.append(succ)
            l += len(Vseqs[succ])

        seq = ""
        for v in precontig:
            seq += Vseqs[v]
        for v in postcontig[1:]:
            seq += Vseqs[v]

        print(
            f">{postcontig[0]}.{w}",
            f"{len(seq) > min_l*2}",
            ">".join([str(x) for x in precontig + postcontig[1:]]),
        )
        print(seq)


if __name__ == "__main__":
    main()
