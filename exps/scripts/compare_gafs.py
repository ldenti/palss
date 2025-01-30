import sys
import glob
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from itertools import chain, combinations
from upsetplot import from_memberships
from upsetplot import plot

sns.set(style="whitegrid")


# https://stackoverflow.com/a/40986475
def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)  # allows duplicate elements
    return chain.from_iterable(combinations(s, r) for r in range(len(s) + 1))


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
        if qidx not in data or cr > data[qidx][0] / data[qidx][1]:
            data[qidx] = [matches, tot, qlen / ql]
    return data


def main():
    wd = sys.argv[1]
    l = float(sys.argv[2])

    gafs = {}
    for f in glob.glob(os.path.join(wd, "*", "alignments*.gaf")):
        n, fname = f.split("/")[-2:]
        n = int(n)
        t = os.path.splitext(fname)[0]
        if t.endswith("w3") or t.endswith("w5"):
            continue

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

        if n not in gafs:
            gafs[n] = {}
        gafs[n][t] = parse_gaf(f, l)

    PW_GRAPHS = {
        1: "FULL/1OUT",
        2: "FULL/FULL-AUG",
        3: "FULL/1OUT-AUG",
        4: "1OUT-AUG/1OUT",
    }
    data = []
    for n in gafs:
        gaf1 = gafs[n]["FULL"]
        gaf2 = gafs[n]["1OUT"]
        gaf3 = gafs[n]["FULL-AUG"]
        gaf4 = gafs[n]["1OUT-AUG"]

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

        # df_membership = from_memberships(Ks, data=Vs)
        # plot(df_membership)
        # plt.show()
        # plt.close()

        # print(
        #     len(set(gaf2) - set(gaf4)),
        #     len(set(gaf2) & set(gaf4)),
        #     len(set(gaf4) - set(gaf2)),
        # )
        # print(set(gaf2) - set(gaf4), set(gaf4) - set(gaf2))

        # print("N", "FULL", "1OUT", "FULL-AUG", "1OUT-AUG")
        # print(n, len(gaf1), len(gaf2), len(gaf3), len(gaf4))
        # print(
        #     n,
        #     len(
        #         set(gaf1.keys())
        #         & set(gaf2.keys())
        #         & set(gaf3.keys())
        #         & set(gaf4.keys())
        #     ),
        # )
        for k in (
            set(gaf1.keys()) & set(gaf2.keys()) & set(gaf3.keys()) & set(gaf4.keys())
        ):
            m1 = gaf1[k][0]
            m2 = gaf2[k][0]
            m3 = gaf3[k][0]
            m4 = gaf4[k][0]

            # if min(m1[1], m2[1]) / max(m1[1], m2[1]) < 0.95:
            #     continue
            # if m1[0] > m2[0]:
            #     if m1[0] - m2[0] > 1000:
            #         print(k, m1, m2)
            data.append([n, PW_GRAPHS[1], m1, m2, m1 - m2])
            data.append([n, PW_GRAPHS[2], m1, m3, m1 - m3])
            data.append([n, PW_GRAPHS[3], m1, m4, m1 - m4])
            data.append([n, PW_GRAPHS[4], m4, m2, m4 - m2])
            if m4 - m2 > 50:
                if gaf2[k][2] > 0.95 and gaf4[k][2] > 0.95:
                    print("===", k, m4 - m2)
                    print(gaf2[k])
                    print(gaf4[k])
                    print("")

    print(len(data))
    df = pd.DataFrame(data, columns=["G", "Graph", "GAF1", "GAF2", "D"])
    G = df["G"].unique()
    G.sort()
    sns.boxplot(data=df, x="G", y="D", hue="Graph", showfliers=False)

    plt.xlabel("#Samples")
    plt.ylabel("ΔMatches")
    plt.tight_layout()
    plt.show()
    plt.savefig("x.png")


# def main2():
#     gaf1_fn = sys.argv[1]  # FULL
#     gaf2_fn = sys.argv[2]  # 1OUT
#     gaf3_fn = sys.argv[3]  # FULL (A)
#     gaf4_fn = sys.argv[4]  # 1OUT (A)
#     N = sys.argv[5]
#     l = float(sys.argv[6])

#     gaf1 = parse_gaf(gaf1_fn, l)
#     gaf2 = parse_gaf(gaf2_fn, l)
#     gaf3 = parse_gaf(gaf3_fn, l)
#     gaf4 = parse_gaf(gaf4_fn, l)

#     Ks = {1: "FULL/1OUT", 2: "FULL/FULL-AUG", 3: "FULL/1OUT-AUG", 4: "1OUT/1OUT-AUG"}
#     data = []
#     nal = 0
#     for k in set(gaf1.keys()) & set(gaf2.keys()) & set(gaf3.keys()) & set(gaf4.keys()):
#         nal += 1
#         m1 = gaf1[k][0]
#         m2 = gaf2[k][0]
#         m3 = gaf3[k][0]
#         m4 = gaf4[k][0]

#         # if min(m1[1], m2[1]) / max(m1[1], m2[1]) < 0.95:
#         #     continue
#         # if m1[0] > m2[0]:
#         #     if m1[0] - m2[0] > 1000:
#         #         print(k, m1, m2)
#         data.append([1, m1, m2])
#         data.append([2, m1, m3])
#         data.append([3, m1, m4])
#         data.append([4, m2, m4])
#         if m3 - m1 > 100:
#             print(k, m3 - m1, gaf1[k], gaf3[k], gaf3[k][1] - gaf1[k][1])
#     print(len(data))
#     df = pd.DataFrame(data, columns=["k", "GAF1", "GAF2"])

#     fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(9, 9))
#     axes = [ax1, ax2, ax3, ax4]
#     for i, ax in enumerate(axes, 1):
#         sns.scatterplot(data=df[df["k"] == i], x="GAF1", y="GAF2", alpha=0.23, ax=ax)

#         corr, _ = pearsonr(df[df["k"] == i]["GAF1"], df[df["k"] == i]["GAF2"])
#         corr = round(corr, 4)
#         ax.set_xlabel(Ks[i].split("/")[0])
#         ax.set_ylabel(Ks[i].split("/")[1])
#         ax.legend([f"PCC={corr}"], loc="lower right")
#         ax.set_xlim([2000, 20000])
#         ax.set_ylim([2000, 20000])
#         ax.set_xticks(np.arange(2000, 20001, 3000))
#         ax.set_yticks(np.arange(2000, 20001, 3000))
#         # ax1.set_title(len(df[df["k"] == 1]))

#     title = f"{N} sample(s) - {nal} alignments"
#     plt.suptitle(title)
#     plt.tight_layout()
#     plt.show()
#     plt.savefig("corr.png")


if __name__ == "__main__":
    main()
