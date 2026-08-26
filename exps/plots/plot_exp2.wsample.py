import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
from functools import reduce
import numpy as np

sns.set_theme()
# sns.set_style("whitegrid")

def get_pr(pdf, rdf, label):
    pdf.loc[pdf["qual"] == -1, "qual"] = 0
    pdf = (
        pdf.groupby(["tool", "n", "coverage"])[["supp", "qual", "l"]]
        .apply(
            lambda x: pd.Series(
                {
                    "P": (x["qual"] * x["l"]).sum() / x["l"].sum(),
                }
            )
        )
        .reset_index()
    )

    rdf = (
        rdf.groupby(["tool", "n", "coverage"])[["nm"]]
        .apply(
            lambda x: pd.Series(
                {
                    "R": (x["nm"] == 0).sum() / x["nm"].count(),
                }
            )
        )
        .reset_index()
    )

    df = pd.merge(pdf, rdf, on=["tool", "n", "coverage"], how="inner")
    df["type"] = label
    return df

def clean_df(df):
    tools = [
        "mgcactus",
        "mgcactus-real",
        "original-full",
        # "palss-d0.1-w2",
        # "palss-d0.25-w2",
        # "palss-d0.33-w2",
        "palss-d0.5-w2",
        "palss-d0.5",
    ]
    df = df[df["tool"].isin(tools)]

    # df.loc[df["tool"].eq("palss-d0.5"), "tool"] = "palss"
    df.loc[df["tool"].eq("palss-d0.5-w2"), "tool"] = "palss"
    df.loc[df["tool"].eq("original-full"), "tool"] = "Truth"
    df.loc[df["tool"].eq("mgcactus"), "tool"] = "MGC"
    df.loc[df["tool"].eq("mgcactus-real"), "tool"] = "MGCr"

    df = df.rename(columns={
        "tool": "Graph",
        "coverage": "Coverage"
    })

    return df

def load_df(pt, rt, sample, c, cov, l):
    pdf = pd.read_csv(pt)
    pdf = pdf[pdf["n"].isin([8, 32, 64])]
    pdf = pdf[pdf["kind"] == "vertex"]
    pdf = pdf[pdf["clen"] == l]
    pdf = pdf[pdf["coverage"] == cov]
    pdf["l"] = pdf["l"].astype(int)

    rdf = pd.read_csv(rt)
    rdf = rdf[rdf["cov"] >= c]
    rdf = rdf[rdf["coverage"] == cov]
    rdf = rdf[rdf["n"].isin([8, 32, 64])]

    df_all = get_pr(pdf, rdf, "All")
    df_all["Sample"] = sample
    df_smp = get_pr(pdf[pdf["type"] == "Simple"], rdf[rdf["type"] == "Simple"], "Simple")
    df_smp["Sample"] = sample
    df_cpx = get_pr(pdf[pdf["type"] == "Complex"], rdf[rdf["type"] == "Complex"], "Complex")
    df_cpx["Sample"] = sample

    df = pd.concat([df_all, df_smp, df_cpx])
    return clean_df(df)

# def load_p(fp, name, cov, l):
#     pdf = pd.read_csv(fp)
#     pdf = pdf[pdf["kind"] == "vertex"]
#     pdf = pdf[pdf["coverage"] == cov]
#     pdf = pdf[pdf["clen"] == l]
#     pdf["l"] = pdf["l"].astype(int)
#     pdf.loc[pdf["qual"] == -1, "qual"] = 0
#     # pdf = pdf[pdf["n"].isin([8, 32, 64])]
#     pdf = (
#         pdf.groupby(["tool", "n", "type", "coverage"])[["supp", "qual", "l"]]
#         .apply(
#             lambda x: pd.Series(
#                 {
#                     "supported_vertices": (x["supp"] > 0).sum(),
#                     "total_vertices": x["supp"].count(),
#                     "P_vertices": (x["supp"] > 0).sum() / x["supp"].count(),
#                     "supported_bases": (x["qual"] * x["l"]).sum(),
#                     "total_bases": x["l"].sum(),
#                     "P_bases": (x["qual"] * x["l"]).sum() / x["l"].sum(),
#                 }
#             )
#         )
#         .reset_index()
#     )
#     pdf["sample"] = name
#     print(pdf.head())
#     return pdf

# def load_r(fp, name, cov):
#     rdf = pd.read_csv(fp)
#     rdf = rdf[rdf["cov"] >= 0.995]
#     pdf = rdf[rdf["coverage"] == cov]
#     # rdf = rdf[rdf["n"].isin([8, 32, 64])]
#     rdf = (
#         rdf.groupby(["tool", "n", "coverage", "type"])[["nm"]]
#         .apply(
#             lambda x: pd.Series(
#                 {
#                     "perfect": (x["nm"] == 0).sum(),
#                     "alignments": x["nm"].count(),
#                     "R": (x["nm"] == 0).sum() / x["nm"].count(),
#                 }
#             )
#         )
#         .reset_index()
#     )
#     rdf["sample"] = name
#     print(rdf.head())
#     return rdf


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("PT1")
    parser.add_argument("RT1")
    parser.add_argument("PT2")
    parser.add_argument("RT2")

    parser.add_argument("-x", default="A", type=str)
    parser.add_argument("-y", default="B", type=str)
    #
    parser.add_argument("--cov", default=0.995, type=float)
    parser.add_argument("--coverage", default=10, type=int)
    # parser.add_argument("-b", "--bases", action="store_true")
    parser.add_argument("-l", default=50000, type=int)
    #
    parser.add_argument("-o", default="", type=str)
    parser.add_argument("--title", default="", type=str)

    args = parser.parse_args()

    df1 = load_df(args.PT1, args.RT1, args.x, args.cov, args.coverage, args.l)
    df2 = load_df(args.PT2, args.RT2, args.y, args.cov, args.coverage, args.l)
    
    df = pd.concat([df1, df2])
    print(df)

    ### === Plotting ===
    Ns = df["n"].unique()
    Ns.sort()
    nn = len(Ns) if len(Ns) > 1 else 2

    # XXX: order is hardcoded
    # x_order = ["Truth", "palss", "MGC", "MGCr"]
    x_order = ["Truth", "palss", "MGC"]
    hue_order = [5, 10, 15]

    # soft_palette = sns.color_palette("pastel", n_colors=len(hue_order))
    # palette = sns.color_palette("bright", n_colors=len(hue_order))

    fig, axes = plt.subplots(3, nn, sharex=True, sharey=True, figsize=(11, 9))
    for row, t in enumerate(["All", "Simple", "Complex"]):
        for col, n in enumerate(Ns):
            ax = axes[row][col]
            subdf = df[(df["n"] == n) & (df["type"] == t)]
            sns.scatterplot(
                data=subdf,
                x="P",
                y="R",
                hue="Graph",
                style="Sample",
                s=75,
                ax=ax,
                legend=True if row == 0 and col == 0 else False,
            )

            # F1 iso-curves
            # xmin, xmax = ax.get_xlim()  # 0.01, 1.05
            # ymin, ymax = ax.get_ylim()  # 0.01, 1.05
            xmin, xmax = 0.01, 1.05
            ymin, ymax = 0.01, 1.05
            recall_grid = np.linspace(xmin, xmax, 500)
            precision_grid = np.linspace(ymin, ymax, 500)
            f_scores = np.linspace(0.2, 0.9, 8)  # F1 levels to draw

            # recall_grid = np.linspace(0.01, 1.05, 500)
            # precision_grid = np.linspace(0.01, 1.05, 500)
            R, P = np.meshgrid(recall_grid, precision_grid)
            F1 = 2 * (P * R) / (P + R)
            cs = ax.contour(
                R,
                P,
                F1,
                levels=f_scores,
                colors="0.6",
                linestyles="dashed",
                linewidths=0.8,
            )
            plt.clabel(cs, inline=1, fmt="F1=%.2f", fontsize=11)

            if row == 0:
                ax.set_title(f"n = {n}")
            if col == 0:
                ax.set_ylabel(f"{t}")
            ax.set_xlabel("")
    fig.supxlabel("Correctness (% supported bases)")
    fig.supylabel("Completeness (% perfect alignments)")

    if args.title != "":
        plt.suptitle(args.title)

    plt.tight_layout()
    if args.o == "":
        plt.show()
    else:
        plt.savefig(args.o)


if __name__ == "__main__":
    main()
