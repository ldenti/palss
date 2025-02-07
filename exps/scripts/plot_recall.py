import sys
import glob
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# from scipy.stats import pearsonr
# from itertools import chain, combinations
# from upsetplot import from_memberships
# from upsetplot import plot

sns.set(style="whitegrid")


# # https://stackoverflow.com/a/40986475
# def powerset(iterable):
#     "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
#     s = list(iterable)  # allows duplicate elements
#     return chain.from_iterable(combinations(s, r) for r in range(len(s) + 1))


def parse_gaf(gaf, l):
    data = {}
    for line in open(gaf):
        line = line.strip("\n").split("\t")
        qidx = line[0]
        matches = int(line[9])
        ql = int(line[1])
        qlen = int(line[3]) - int(line[2])
        if qlen / ql < l:
            continue
        tot = int(line[10])
        cr = matches / tot

        nm = line[13]
        nm = nm.split(":")
        assert nm[0] == "NM"
        nm = int(nm[2])

        if qidx not in data or cr > data[qidx][0] / data[qidx][1]:
            data[qidx] = [matches, tot, qlen / ql, nm]
    return data


def main():
    wd = sys.argv[1]
    l = float(sys.argv[2])

    gafs = {}
    for f in glob.glob(os.path.join(wd, "*", "alignments*.gaf")):
        n, fname = f.split("/")[-2:]
        n = int(n)

        t = os.path.splitext(fname)[0]

        if "mgcactus" in t:
            t = "MGC"
        else:
            if "augmented" not in t:
                if "1out" in t:
                    t = "1OUT"
                else:
                    t = "FULL"
            else:
                if "1out" in t:
                    t = "1OUT-AUG"
                else:
                    t = "FULL-AUG"

        if t == "FULL-AUG":
            continue

        if n not in gafs:
            gafs[n] = {}
        gafs[n][t] = parse_gaf(f, l)

    graphs = ["1OUT", "1OUT-AUG", "FULL", "MGC"]  # "FULL-AUG"
    data = []
    for n in gafs:
        for g in graphs:
            for d in gafs[n][g].values():
                data.append([n, g, d[-1], d[-2]])

    df = pd.DataFrame(data, columns=["G", "Graph", "NM", "%Aligned"])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5), dpi=80)

    x = sns.boxplot(
        data=df,
        x="G",
        y="NM",
        hue="Graph",
        hue_order=graphs,
        showfliers=False,
        palette="Set1",
        ax=ax1,
    )
    ax1.legend(loc=1, ncols=2)
    ax1.set_xlabel("n")
    yl = x.get_ylim()
    ax1.set_yticks(np.arange(-5, int(yl[1]) + (int(yl[1]) % 5), 5))
    ax1.set_title("(a)")

    PW_GRAPHS = {
        1: "FULL/1OUT",
        # 2: "FULL/FULL-AUG",
        3: "FULL/1OUT-AUG",
        # 4: "FULL/MGC",
        5: "1OUT-AUG/1OUT",
        6: "1OUT-AUG/MGC",
        7: "MGC/1OUT",
    }
    data = []
    for n in gafs:
        gaf1 = gafs[n]["FULL"]
        gaf2 = gafs[n]["1OUT"]
        # gaf3 = gafs[n]["FULL-AUG"]
        gaf4 = gafs[n]["1OUT-AUG"]
        gaf5 = gafs[n]["MGC"]

        # === upset plot =======
        # keys = ["FULL", "1OUT", "FULL-AUG", "1OUT-AUG"]
        # union = set.union(*[set(gafs[n][k]) for k in keys])
        # print(len(union))
        # Ks = []
        # Vs = []
        # for comb in powerset(keys):
        #     if len(comb) == 0:
        #         continue
        #     inside = set.intersection(*[set(gafs[n][k]) for k in comb])
        #     if len(comb) == 4:
        #         print(len(inside))
        #         continue
        #     outside = set.union(*[set(gafs[n][k]) for k in keys if k not in comb])
        #     print(comb, inside - outside)
        #     Ks.append(comb)
        #     Vs.append(len(inside - outside))
        #
        # df_membership = from_memberships(Ks, data=Vs)
        # plot(df_membership)
        # plt.show()
        # plt.close()

        for k in (
            set(gaf1.keys())
            | set(gaf2.keys())
            # | set(gaf3.keys())
            | set(gaf4.keys())
            | set(gaf5.keys())
        ):
            m1 = gaf1[k][0] if k in gaf1 else None
            m2 = gaf2[k][0] if k in gaf2 else None
            # m3 = gaf3[k][0] if k in gaf3 else None
            m4 = gaf4[k][0] if k in gaf4 else None
            m5 = gaf5[k][0] if k in gaf5 else None

            if m1 != None and m2 != None:
                data.append([n, PW_GRAPHS[1], m1, m2, m1 - m2])
            # if m1 != None and m3 != None:
            #     data.append([n, PW_GRAPHS[2], m1, m3, m1 - m3])
            if m1 != None and m4 != None:
                data.append([n, PW_GRAPHS[3], m1, m4, m1 - m4])
            # if m1 != None and m5 != None:
            #     data.append([n, PW_GRAPHS[4], m1, m5, m1 - m5])
            if m2 != None and m4 != None:
                data.append([n, PW_GRAPHS[5], m4, m2, m4 - m2])
            if m4 != None and m5 != None:
                data.append([n, PW_GRAPHS[6], m4, m5, m4 - m5])
            # if m2 != None and m5 != None:
            #     data.append([n, PW_GRAPHS[7], m5, m2, m5 - m2])

            if n == 1 and m2 != None and m4 != None and m4 - m2 < 0:
                print("===", k, m4 - m2)
                print(gaf2[k])
                print(gaf4[k])
                print("")

    print(len(data))
    df = pd.DataFrame(data, columns=["G", "Graph", "GAF1", "GAF2", "D"])
    G = df["G"].unique()
    G.sort()
    sns.boxplot(
        data=df, x="G", y="D", hue="Graph", showfliers=False, palette="Set2", ax=ax2
    )

    ax2.legend(loc=1, ncols=1)

    ax2.set_xlabel("n")
    ax2.set_ylabel("ΔMatches")
    ax2.set_title("(b)")

    plt.tight_layout()
    plt.savefig("x.png", bbox_inches="tight")
    plt.savefig("x.pdf", bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()
