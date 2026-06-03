import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns


def hist_plot(df, palss=False):
    Ns = df["n"].unique()
    Ns.sort()

    fig, axes = plt.subplots(1, len(Ns), sharey=True, figsize=(13, 7))
    for i, n in enumerate(Ns):
        sns.histplot(
            df[df["n"] == n],
            x="nm",
            hue="Graph",
            discrete=True,
            element="step",
            fill=False,
            legend=True if i == 0 else False,
            binrange=[0, 10],
            palette="bright",
            ax=axes[i],
        )
        axes[i].set_ylabel("")
        axes[i].set_title(n)

    sns.move_legend(
        axes[0],
        "upper right",
        # bbox_to_anchor=(0.5, 1),
        # ncol=3,
        title=None,
        # frameon=False,
    )
    plt.suptitle("NM (per best/longest alignment)")
    plt.tight_layout()
    plt.show()
    # plt.savefig("nm-hist.png")


def dump_wide(df):
    values = [0, 1, 2]
    df = df[["n", "Graph", "nm"]]
    df = df[df["nm"].isin(values)]
    df["t"] = df["Graph"]
    df = df.pivot_table(
        values="nm", index=["Graph", "n"], columns="nm", aggfunc="count", fill_value=-1
    ).reset_index()
    df.columns.name = None
    df = df.rename(columns={x: str(x) for x in values})
    print(df)


def bar_plot(df):
    Ns = df["n"].unique()
    Ns.sort()

    fig, axes = plt.subplots(1, len(Ns), sharey=True, figsize=(13, 7))
    for i, n in enumerate(Ns):
        sub_df = df[df["n"] == n]
        fdf = sub_df[sub_df["nm"].isin([0, 1])]
        sns.countplot(
            data=fdf,
            x="nm",
            hue="Graph",
            legend=True if i == 0 else False,
            palette="bright",
            ax=axes[i],
        )
        axes[i].set_ylabel("")
        axes[i].set_title(n)

    sns.move_legend(axes[0], "lower left")
    plt.tight_layout()
    plt.show()
    # plt.savefig("nm-bar.png")


def bar_plot_palss(df):
    df["full"] = df["Graph"].str.contains("full", case=False)

    Ns = df["n"].unique()
    Ns.sort()

    fig, axes = plt.subplots(2, len(Ns), sharey=True, figsize=(13, 7))
    for i, n in enumerate(Ns):
        sub_df = df[df["n"] == n]
        fdf = sub_df[sub_df["nm"].isin([0, 1])]
        for j, f in enumerate([True, False]):
            sns.countplot(
                data=fdf[fdf["full"] == j],
                x="nm",
                hue="d",
                legend=True if i == 0 else False,
                palette="bright",
                ax=axes[j][i],
            )
            # axes[0][i].set_ylabel("")
            # axes[0][i].set_title(n)

    # sns.move_legend(axes[0], "lower left")
    plt.tight_layout()
    plt.show()
    # plt.savefig("nm-bar-palss.png")


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("CSV")
    parser.add_argument("-p", "--palss", action="store_true")
    parser.add_argument("-c", "--cov", type=float, default=0.9)
    args = parser.parse_args()

    fn = args.CSV
    df = pd.read_csv(fn)

    df = df[df["cov"] >= args.cov]

    if args.palss:
        df = df[df["graph"] == "palss"]
    else:
        df = df[(df["d"] == 0.1) | (df["d"] == -1)]

    df.loc[df["fn"].str.contains("palss-full"), "Graph"] = "PALSS (full)"
    df.loc[df["fn"].str.contains("palss-oneout"), "Graph"] = "PALSS (1out)"
    df.loc[df["fn"].str.contains("mgcactus"), "Graph"] = "MGCACTUS"
    df.loc[df["fn"].str.contains("original-full"), "Graph"] = "Original"
    df.loc[df["fn"].str.contains("original-oneout"), "Graph"] = "1out"
    df.loc[df["fn"].str.contains("reads.bam"), "Graph"] = "Reference"
    df.loc[df["fn"].str.contains("reads.tohaps.bam"), "Graph"] = "Haplotypes"

    df = df.drop(columns="fn")

    if not args.palss:
        hist_plot(df, args.palss)
    dump_wide(df)
    if args.palss:
        bar_plot_palss(df)
    else:
        bar_plot(df)


if __name__ == "__main__":
    main()
