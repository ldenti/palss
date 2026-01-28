import sys
import os
import glob


def parse_gaf(gaf_fn):
    nms = {}
    for line in open(gaf_fn):
        line = line.strip("\n").split("\t")
        qidx = line[0]
        ql = int(line[1])
        qs = int(line[2])  # closed
        qe = int(line[3])  # open

        nm = line[12]
        nm = nm.split(":")
        assert nm[0] == "NM"
        nm = int(nm[2])
        c = (qe - qs) / ql
        if qidx not in nms:
            nms[qidx] = (c, nm)
        else:
            old_c, old_nm = nms[qidx]
            if old_c == c:
                nms[qidx] = (c, min(nm, old_nm))
            elif old_c < c:
                nms[qidx] = (c, nm)
            else:
                pass
    return nms


def main():
    WD = sys.argv[1]

    print("fn,graph,n,w,d,read,cov,nm")
    data = []
    for gaf_fn in glob.glob(os.path.join(WD, "n*", "graphaligner", "*.gaf")):
        n = int(gaf_fn.split("/")[-3][1:])
        fn = gaf_fn.split("/")[-1]
        w, d = -1, -1
        graph = ""
        run = "full" if "full" in fn else "oneout"
        if "original" in fn:
            graph = "original"
        elif "mgcactus" in fn:
            graph = "mgcactus"
        else:
            graph = "palss"
            w = int(fn.split(".")[-2][1:])
            d = float(fn.split(".")[-4][1:] + "." + fn.split(".")[-3])

        nms = parse_gaf(gaf_fn)
        for qidx, (c, nm) in nms.items():
            print(fn, graph, n, w, d, qidx, c, nm, sep=",", flush=False)
        sys.stdout.flush()


if __name__ == "__main__":
    main()
