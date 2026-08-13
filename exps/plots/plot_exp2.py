import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns

sns.set_theme()
# sns.set_style("whitegrid")


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

    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("PT")
    parser.add_argument("RT")
    parser.add_argument("-c", default=0.995, type=float)
    parser.add_argument("-t", default="simple", type=str)
    parser.add_argument("-l", default=50000, type=int)
    parser.add_argument("-s", action="store_true")
    #
    parser.add_argument("-o", default="", type=str)
    parser.add_argument("--title", default="", type=str)

    args = parser.parse_args()

    pdf = pd.read_csv(args.PT)
    pdf = pdf[pdf["n"].isin([8, 32, 64])]
    if args.t == "simple":
        pdf = pdf[pdf["type"] == "Simple"]
    elif args.t == "complex":
        pdf = pdf[pdf["type"] == "Complex"]
    else:
        pass

    rdf = pd.read_csv(args.RT)
    rdf = rdf[rdf["cov"] >= args.c]
    rdf = rdf[rdf["n"].isin([8, 32, 64])]
    if args.t == "simple":
        rdf = rdf[rdf["type"] == "Simple"]
    elif args.t == "complex":
        rdf = rdf[rdf["type"] == "Complex"]
    else:
        pass

    ### === PRECISION ===
    pdf_v = pdf[pdf["kind"] == "vertex"]
    pdf_v = pdf_v[pdf_v["clen"] == args.l]
    pdf_v["l"] = pdf_v["l"].astype(int)
    pdf_v.loc[pdf_v["qual"] == -1, "qual"] = 0
    pdf_v = (
        pdf_v.groupby(["tool", "n", "coverage"])[["supp", "qual", "l"]]
        .apply(
            lambda x: pd.Series(
                {
                    "supported_vertices": (x["supp"] > 0).sum(),
                    "total_vertices": x["supp"].count(),
                    "P_vertices": (x["supp"] > 0).sum() / x["supp"].count(),
                    "supported_bases": (x["qual"] * x["l"]).sum(),
                    "total_bases": x["l"].sum(),
                    "P_bases": (x["qual"] * x["l"]).sum() / x["l"].sum(),
                }
            )
        )
        .reset_index()
    )

    pdf_e = pdf[pdf["kind"] == "edge"]
    pdf_e = (
        pdf_e.groupby(["tool", "n", "coverage"])[["supp", "qual", "l"]]
        .apply(
            lambda x: pd.Series(
                {
                    "supported_edges": (x["supp"] > 0).sum(),
                    "total_edges": x["supp"].count(),
                    "P_edges": (x["supp"] > 0).sum() / x["supp"].count(),
                }
            )
        )
        .reset_index()
    )

    # print(pdf_v.head())
    # print(pdf_e.head())

    ### === RECALL ===
    rdf = (
        rdf.groupby(["tool", "n", "coverage"])[["nm"]]
        .apply(
            lambda x: pd.Series(
                {
                    "perfect": (x["nm"] == 0).sum(),
                    "alignments": x["nm"].count(),
                    "R": (x["nm"] == 0).sum() / x["nm"].count(),
                }
            )
        )
        .reset_index()
    )
    # print(rdf.head())

    ### === Plotting ===
    pdf_v = clean_df(pdf_v)
    pdf_e = clean_df(pdf_e)
    rdf = clean_df(rdf)

    Ns = pdf_v["n"].unique()
    Ns.sort()
    nn = len(Ns) if len(Ns) > 1 else 2

    # XXX: order is hardcoded
    x_order = ["Truth", "palss", "MGC", "MGCr"]
    hue_order = [5, 10, 15]

    soft_palette = sns.color_palette("pastel", n_colors=len(hue_order))
    palette = sns.color_palette("bright", n_colors=len(hue_order))

    fig, axes = plt.subplots(4, nn, sharex=True, sharey="row", figsize=(11, 9))
    if not args.s:
        for col, n in enumerate(Ns):
            ### Precision with number of vertices ###
            row = 0
            ax = sns.barplot(
                data=pdf_v[(pdf_v["n"] == n)],
                x="tool",
                order=x_order,
                hue="coverage",
                y="total_vertices",
                palette=soft_palette,
                legend=False,
                ax=axes[row][col],
            )

            # annotation
            lookup = (
                pdf_v[(pdf_v["n"] == n)]
                .set_index(["tool", "coverage"])["P_vertices"]
                .to_dict()
            )
            lookup = {k: round(v * 100, 1) for k, v in lookup.items()}

            x = len(x_order)
            h = len(hue_order)
            for i in range(x * h):
                # ii = i + x * h
                p = ax.patches[i]

                # XXX: here we use x since it's larger
                xcat = x_order[i % x]
                huecat = hue_order[i // x]
                ratio_val = lookup[(xcat, huecat)]
                # print(xcat, huecat, ratio_val)

                ax.annotate(
                    f"{ratio_val:.1f}",
                    (p.get_x() + p.get_width() / 2, p.get_height()),
                    ha="center",
                    va="bottom",
                    rotation=90,
                    xytext=(1, 3),
                    textcoords="offset points",
                )

            sns.barplot(
                data=pdf_v[(pdf_v["n"] == n)],
                x="tool",
                order=x_order,
                hue="coverage",
                hue_order=hue_order,
                y="supported_vertices",
                palette=palette,
                legend=False,
                ax=axes[row][col],
            )

            ### Precision with bases ###
            row = 1
            ax = sns.barplot(
                data=pdf_v[(pdf_v["n"] == n)],
                x="tool",
                order=x_order,
                hue="coverage",
                y="total_bases",
                palette=soft_palette,
                legend=False,
                ax=axes[row][col],
            )

            # annotation
            lookup = (
                pdf_v[(pdf_v["n"] == n)]
                .set_index(["tool", "coverage"])["P_bases"]
                .to_dict()
            )
            lookup = {k: round(v * 100, 1) for k, v in lookup.items()}

            x = len(x_order)
            h = len(hue_order)
            for i in range(x * h):
                # ii = i + x * h
                p = ax.patches[i]

                # XXX: here we use x since it's larger
                xcat = x_order[i % x]
                huecat = hue_order[i // x]
                ratio_val = lookup[(xcat, huecat)]

                ax.annotate(
                    f"{ratio_val:.1f}",
                    (p.get_x() + p.get_width() / 2, p.get_height()),
                    ha="center",
                    va="bottom",
                    rotation=90,
                    xytext=(1, 3),
                    textcoords="offset points",
                )

            sns.barplot(
                data=pdf_v[(pdf_v["n"] == n)],
                x="tool",
                order=x_order,
                hue="coverage",
                y="supported_bases",
                palette=palette,
                legend=False,
                ax=axes[row][col],
            )

            ### Precision with edges ###
            row = 2
            ax = sns.barplot(
                data=pdf_e[(pdf_e["n"] == n)],
                x="tool",
                order=x_order,
                hue="coverage",
                y="total_edges",
                palette=soft_palette,
                legend=False,
                ax=axes[row][col],
            )

            # annotation
            lookup = (
                pdf_e[(pdf_e["n"] == n)]
                .set_index(["tool", "coverage"])["P_edges"]
                .to_dict()
            )
            lookup = {k: round(v * 100, 1) for k, v in lookup.items()}

            x = len(x_order)
            h = len(hue_order)
            for i in range(x * h):
                # ii = i + x * h
                p = ax.patches[i]

                # XXX: here we use x since it's larger
                xcat = x_order[i % x]
                huecat = hue_order[i // x]
                ratio_val = lookup[(xcat, huecat)]

                ax.annotate(
                    f"{ratio_val:.1f}",
                    (p.get_x() + p.get_width() / 2, p.get_height()),
                    ha="center",
                    va="bottom",
                    rotation=90,
                    xytext=(1, 3),
                    textcoords="offset points",
                )

            sns.barplot(
                data=pdf_e[(pdf_e["n"] == n)],
                x="tool",
                order=x_order,
                hue="coverage",
                y="supported_edges",
                palette=palette,
                legend=False,
                ax=axes[row][col],
            )

            ### Recall ###
            row = 3
            ax = sns.barplot(
                data=rdf[(rdf["n"] == n)],
                x="tool",
                order=x_order,
                hue="coverage",
                y="alignments",
                palette=soft_palette,
                legend=False,
                ax=axes[row][col],
            )

            # annotation
            lookup = rdf[(rdf["n"] == n)].set_index(["tool", "coverage"])["R"].to_dict()
            lookup = {k: round(v * 100, 1) for k, v in lookup.items()}

            x = len(x_order)
            h = len(hue_order)
            for i in range(x * h):
                # ii = i + x * h
                p = ax.patches[i]

                # XXX: here we use x since it's larger
                xcat = x_order[i % x]
                huecat = hue_order[i // x]
                ratio_val = lookup[(xcat, huecat)]

                ax.annotate(
                    f"{ratio_val:.1f}",
                    (p.get_x() + p.get_width() / 2, p.get_height()),
                    ha="center",
                    va="bottom",
                    rotation=90,
                    xytext=(1, 3),
                    textcoords="offset points",
                )

            sns.barplot(
                data=rdf[(rdf["n"] == n)],
                x="tool",
                order=x_order,
                hue="coverage",
                y="perfect",
                palette=palette,
                legend=False,
                ax=axes[row][col],
            )
    else:
        for col, n in enumerate(Ns):
            # Precision with number of vertices
            row = 0
            ax = sns.barplot(
                data=pdf_v[(pdf_v["n"] == n)],
                x="tool",
                order=x_order,
                hue="coverage",
                y="P_vertices",
                palette=palette,
                legend=False,
                ax=axes[row][col],
            )
            for p in ax.patches:
                ax.annotate(
                    f"{round(p.get_height()*100, 1):.1f}",
                    (p.get_x() + p.get_width() / 2, p.get_height()),
                    ha="center",
                    va="bottom",
                    rotation=90,
                    xytext=(1, 3),
                    textcoords="offset points",
                )

            # Precision with bases
            row = 1
            ax = sns.barplot(
                data=pdf_v[(pdf_v["n"] == n)],
                x="tool",
                order=x_order,
                hue="coverage",
                y="P_bases",
                palette=palette,
                legend=False,
                ax=axes[row][col],
            )
            for p in ax.patches:
                ax.annotate(
                    f"{round(p.get_height()*100, 1):.1f}",
                    (p.get_x() + p.get_width() / 2, p.get_height()),
                    ha="center",
                    va="bottom",
                    rotation=90,
                    xytext=(1, 3),
                    textcoords="offset points",
                )

            # Precision with edges
            row = 2
            ax = sns.barplot(
                data=pdf_e[(pdf_e["n"] == n)],
                x="tool",
                order=x_order,
                hue="coverage",
                y="P_edges",
                palette=palette,
                legend=False,
                ax=axes[row][col],
            )
            for p in ax.patches:
                ax.annotate(
                    f"{round(p.get_height()*100, 1):.1f}",
                    (p.get_x() + p.get_width() / 2, p.get_height()),
                    ha="center",
                    va="bottom",
                    rotation=90,
                    xytext=(1, 3),
                    textcoords="offset points",
                )

            # Recall
            row = 3
            ax = sns.barplot(
                data=rdf[(rdf["n"] == n)],
                x="tool",
                order=x_order,
                hue="coverage",
                y="R",
                palette=palette,
                legend=False,
                ax=axes[row][col],
            )
            for p in ax.patches:
                ax.annotate(
                    f"{round(p.get_height()*100, 1):.1f}",
                    (p.get_x() + p.get_width() / 2, p.get_height()),
                    ha="center",
                    va="bottom",
                    rotation=90,
                    xytext=(1, 3),
                    textcoords="offset points",
                )

    if args.title != "":
        plt.suptitle(args.title)

    # XXX: a lot of hardcoded stuff ##################################

    axes[0][0].set_ylabel("P (vertices)")
    axes[1][0].set_ylabel("P (bases)")
    axes[2][0].set_ylabel("P (edges)")
    axes[3][0].set_ylabel("R")
    for i in range(3):
        axes[3][i].set_xlabel(f"Graph (n={Ns[i]})")

    labels = [f"{x}x" for x in hue_order]
    palette = sns.color_palette(palette, n_colors=len(labels))
    color_map = dict(zip(labels, palette))
    handles = [Patch(facecolor=color_map[l], edgecolor="none", label=l) for l in labels]
    fig.legend(
        handles=handles,
        title="Coverage",
        loc="lower center",
        bbox_to_anchor=(0.5375, -0.01),
        ncol=len(handles),
        frameon=False,
    )
    fig.tight_layout(rect=[0, 0.05, 1, 1])

    ##################################################################

    # plt.tight_layout()
    if args.o == "":
        plt.show()
    else:
        plt.savefig(args.o)


if __name__ == "__main__":
    main()
