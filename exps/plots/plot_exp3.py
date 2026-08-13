import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# import venn
# from matplotlib_venn import venn3


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("CSV")

    parser.add_argument("-c", type=float, default=0.995)
    # parser.add_argument("-n", type=int, default=8)
    parser.add_argument("-m", type=int, default=5)
    parser.add_argument("--peak", action="store_true")
    # parser.add_argument("--better", action="store_true")
    # parser.add_argument("--noruns", action="store_true")
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

    wide = wide[wide["NM_assembly"] < wide["NM_original"]]

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

    for g in graphs:
        zeros = len(subdf[(subdf["Reference"] == g) & (subdf["NM"] == 0)])
        size = len(subdf[subdf["Reference"] == g])
        print(g, zeros, size, zeros / size)

    fig, (ax1, ax2) = plt.subplots(1, 2)

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
    plt.tight_layout()
    if args.o == "":
        plt.show()
    else:
        plt.savefig(args.o)

    # plt.close()

    # # df = subdf

    # zeros_assembly = set(
    #     df[(df["tool"] == "assembly") & (df["NM"] == 0)]["name"].unique()
    # )
    # zeros_palss0 = set(
    #     df[(df["tool"] == "palss-s0") & (df["NM"] == 0)]["name"].unique()
    # )
    # zeros_palss3 = set(
    #     df[(df["tool"] == "palss-s3") & (df["NM"] == 0)]["name"].unique()
    # )
    # zeros_original = set(
    #     df[(df["tool"] == "original") & (df["NM"] == 0)]["name"].unique()
    # )

    # labels = venn.get_labels(
    #     [zeros_assembly, zeros_palss0, zeros_palss3, zeros_original],
    #     fill=["number"],  # , "logic"],
    # )
    # fig, ax = venn.venn4(labels, names=("Assembly", "palss-s0", "palss-s3", "original"))
    # plt.show()

    # plt.close()

    # df = subdf

    # zeros_assembly = set(
    #     df[(df["tool"] == "assembly") & (df["NM"] == 0)]["name"].unique()
    # )
    # zeros_palss0 = set(
    #     df[(df["tool"] == "palss-s0") & (df["NM"] == 0)]["name"].unique()
    # )
    # zeros_palss3 = set(
    #     df[(df["tool"] == "palss-s3") & (df["NM"] == 0)]["name"].unique()
    # )

    # labels = venn.get_labels(
    #     [zeros_assembly, zeros_palss0, zeros_palss3],
    #     fill=["number"],  # , "logic"],
    # )
    # fig, ax = venn.venn3(labels, names=("Assembly", "palss-s0", "palss-s3"))

    # plt.show()


if __name__ == "__main__":
    main()
