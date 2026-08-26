import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import FuncFormatter

import venn

# from matplotlib_venn import venn2
# from matplotlib_venn.layout.venn2 import DefaultLayoutAlgorithm

sns.set_theme()
# sns.set_style("whitegrid")


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("CSV")
    parser.add_argument("-c", type=float, default=0.995)
    parser.add_argument("-m", type=int, default=5)
    parser.add_argument("--peak", action="store_true")
    parser.add_argument("-o", type=str, default="")

    args = parser.parse_args()

    df = pd.read_csv(args.CSV)

    print(df["Coverage"].gt(args.c).mean() * 100)
    print(df.describe(percentiles=[0.05, 0.1, 0.25, 0.5, 0.75]))

    df = df[df["Coverage"] >= args.c]

    graphs = df["Reference"].unique()
    for g in graphs:
        zeros = len(df[(df["Reference"] == g) & (df["NM"] == 0)])
        size = len(df[df["Reference"] == g])
        print(g, zeros, size, zeros / size)

    df.loc[df["NM"].values > args.m, "NM"] = args.m
    if not args.peak:
        df = df[df["NM"] != args.m]

    wide = df.pivot_table(
        index="Read",
        columns="Reference",
        values=["NM"],
        aggfunc="first",
    )
    # flatten columns: nm_<tool>, run_<tool>
    wide.columns = [f"{metric}_{tool}" for metric, tool in wide.columns]
    wide = wide.reset_index()

    wide = wide[wide["NM_Assembly"] < wide["NM_HPRCv2"]]

    wide_long = wide.melt(id_vars="Read", var_name="var", value_name="value")

    # split "nm_tool" or "run_tool" into metric and tool
    wide_long[["metric", "Reference"]] = wide_long["var"].str.split(
        "_", n=1, expand=True
    )

    # metric -> nm/run; value -> numeric/string; then pivot back to columns
    subdf = (
        wide_long.pivot_table(
            index=["Read", "Reference"],
            columns="metric",
            values="value",
            aggfunc="first",
        )
        .reset_index()
        .rename_axis(None, axis=1)
    )

    subdf = subdf.dropna(subset=["NM"])
    # if args.noruns:
    #     df = df[df["run"] == False]

    # for m in range(args.m):
    for g in graphs:
        zeros = len(subdf[(subdf["Reference"] == g) & (subdf["NM"] == 0)])
        size = len(subdf[subdf["Reference"] == g])
        print(0, g, zeros, size, zeros / size)

    # refs = ["palss-00", "original", "assembly", "reference"]
    # df = df[df["Reference"].isin(refs)]
    # subdf = subdf[subdf["Reference"].isin(refs)]

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(10, 4))

    sns.histplot(
        data=df,
        x="NM",
        hue="Reference",
        hue_order=graphs,
        element="poly",
        discrete=True,
        fill=False,
        cumulative=True,
        ax=ax1,
    )

    ax1.set_xlabel("NM ≤ x")
    ax1.set_ylabel("Number of alignments")
    ax1.yaxis.set_major_formatter(
        FuncFormatter(
            lambda value, position: f"{value / 1000:g}k" if value != 0 else "0"
        )
    )

    sns.histplot(
        data=subdf,
        x="NM",
        hue="Reference",
        hue_order=graphs,
        element="poly",
        discrete=True,
        fill=False,
        legend=False,
        cumulative=True,
        ax=ax2,
    )
    ax2.set_xlabel("NM ≤ x")
    ax2.set_ylabel("")
    ax2.yaxis.set_major_formatter(
        FuncFormatter(
            lambda value, position: f"{value / 1000:g}k" if value != 0 else "0"
        )
    )
    # plt.tight_layout()
    # plt.show()
    # plt.close()

    # fig, ax = plt.subplots(1, 1, figsize=(6, 6))

    zeros_assembly = set(
        subdf[(subdf["Reference"] == "Assembly") & (subdf["NM"] == 0)]["Read"].unique()
    )
    zeros_palss = set(
        subdf[(subdf["Reference"] == "palss") & (subdf["NM"] == 0)]["Read"].unique()
    )

    # total = len(zeros_assembly | zeros_palss)
    # venn2(
    #     (zeros_assembly, zeros_palss),
    #     set_labels=("Assembly", "palss"),
    #     ax=ax3,
    #     layout_algorithm=DefaultLayoutAlgorithm(fixed_subset_sizes=(2, 2, 2)),
    #     subset_label_formatter=lambda v: f"{v}\n{v/total:.1%}",
    # )

    dataset_dict = {"Assembly": zeros_assembly, "palss": zeros_palss}
    venn.venn(
        dataset_dict,
        fmt="{size}\n({percentage:.1f}%)",
        cmap=sns.color_palette()[2:],
        fontsize=10,
        legend_loc="upper left",
        ax=ax3,
        # legend_loc=None,
        figsize=(9, 7),
    )
    ax3.set_aspect("equal", adjustable="box")
    ax3.margins(0)
    ax3.axis("off")

    # fig = ax3.get_figure()
    # fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    # ax3.set_position([0, 0, 1, 1])
    # ax3.axis("off")

    ax1.set_title("(a)")
    ax2.set_title("(b)")
    ax3.set_title("(c)", pad=18)

    # ax1.set_anchor("N")
    # ax2.set_anchor("N")
    # ax3.set_anchor("N")

    # ax3.set_aspect("auto")

    # fig.align_titles()
    # plt.tight_layout(pad=0)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
