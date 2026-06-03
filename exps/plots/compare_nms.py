import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("CSV")
    parser.add_argument("-z", "--zeros", action="store_true")
    parser.add_argument("-g", "--greater", action="store_true")
    parser.add_argument("-l", "--lower", action="store_true")
    parser.add_argument("-a", required=True, type=str)
    parser.add_argument("-b", required=True, type=str)
    args = parser.parse_args()

    fn = args.CSV
    t1 = args.a
    t2 = args.b

    df = pd.read_csv(fn)

    df = df[(df["fn"] == t1) | (df["fn"] == t2)]
    df = df[df["cov"] >= 0.9]

    df = df[["fn", "n", "graph", "read", "nm"]]

    g1 = df[df["fn"] == t1]["graph"].unique()[0]
    g2 = df[df["fn"] == t2]["graph"].unique()[0]

    df.drop("fn", axis=1)

    Ns = df["n"].unique()

    wide_df = (
        df.pivot(index=["n", "read"], columns="graph", values="nm")
        .reset_index()
        .dropna()
    )

    for n in Ns:
        print(len(wide_df[wide_df["n"] == n]))
        print(len(wide_df[(wide_df["n"] == n) & (wide_df[g2] == 0)]))
        print(len(wide_df[(wide_df["n"] == n) & (wide_df[g1] == 0)]))
        print("")

    wide_df = wide_df[wide_df[g1] != wide_df[g2]]

    if args.zeros:
        wide_df = wide_df[(wide_df[g1] == 0) | (wide_df[g2] == 0)]
    elif args.greater:
        wide_df = wide_df[wide_df[g1] > wide_df[g2]]
    elif args.lower:
        wide_df = wide_df[wide_df[g1] < wide_df[g2]]
    print(wide_df.describe())

    """
    # Sorted scatterplot
    fig, axes = plt.subplots(
        len(Ns),
        2,
        sharey=True,
        figsize=(13, 7),
    )
    print(Ns)
    title = ""
    for i, n in enumerate(Ns):
        sorted_df = wide_df[wide_df["n"] == n].sort_values(g1)
        title += str(len(sorted_df)) + "/"
        sns.lineplot(
            data=sorted_df,
            x="read",
            y=g2,
            label=g2,
            color="red",
            # color="white",
            # markeredgecolor="red",
            marker=".",
            linestyle="None",
            ax=axes[i][0],
        )
        sns.lineplot(
            data=sorted_df,
            x="read",
            y=g1,
            label=g1,
            color="green",
            marker="None",
            # linestyle="None",
            ax=axes[i][0],
        )

        sorted_df = wide_df[wide_df["n"] == n].sort_values(g2)
        sns.lineplot(
            data=sorted_df,
            x="read",
            y=g1,
            label=g1,
            color="green",
            marker=".",
            linestyle="None",
            ax=axes[i][1],
        )
        sns.lineplot(
            data=sorted_df,
            x="read",
            y=g2,
            label=g2,
            color="red",
            marker=None,
            # linestyle="None",
            ax=axes[i][1],
        )

        axes[i][0].set_xticks([])
        axes[i][0].set_ylabel("NM")
        axes[i][1].set_xticks([])
        if i == len(Ns) - 1:
            axes[i][0].set_xlabel("Alignments")
            axes[i][1].set_xlabel("Alignments")
        else:
            axes[i][0].set_xlabel("")
            axes[i][1].set_xlabel("")
    plt.suptitle(title[:-1] + " alignments")
    plt.tight_layout()
    plt.show()
    """

    print(df)
    # df["nm"] = np.log(df["nm"])

    # reshaped_df = df.melt(
    #     id_vars="n", value_vars=["palss", "mgcactus"], var_name="Tool", value_name="nm"
    # )
    wide_df["delta"] = wide_df[g1] - wide_df[g2]
    # wide_df["neg"] = wide_df["delta"] < 0
    # wide_df["delta"] = wide_df["delta"].abs()

    T = 25
    fig, axes = plt.subplots(2, len(Ns), figsize=(13, 13))
    for i, n in enumerate(Ns):
        ndf = wide_df[wide_df["n"] == n]

        ### First row: scatterplot

        axes[0][i].axline((0, 0), slope=1, color="black")
        sns.scatterplot(
            ndf,
            x=g1,
            y=g2,
            facecolor="none",
            edgecolor="blue",
            # alpha=0.5,
            s=7,
            ax=axes[0][i],
        )

        axes[0][i].set_xscale("log")
        axes[0][i].set_yscale("log")

        ### Second row: barplot

        counts = ndf["delta"].value_counts()

        full_range = pd.Series(
            range(int(min(ndf["delta"])), int(max(ndf["delta"])) + 1)
        )
        counts = counts.reindex(full_range, fill_value=0)
        counts = counts.reset_index()

        counts.columns = ["delta", "count"]

        print(counts)
        sum_1 = counts.loc[counts["delta"] >= T, "count"].sum()
        sum_2 = counts.loc[counts["delta"] <= -T, "count"].sum()

        counts = counts.drop(
            counts[(counts["delta"] >= T) | (counts["delta"] <= -T)].index
        )

        new_rows = pd.DataFrame({"delta": [-T - 1, T + 1], "count": [sum_2, sum_1]})
        counts = pd.concat([counts, new_rows], ignore_index=True)
        counts = counts.drop(counts[counts["delta"] == 0].index)
        print(counts)
        counts["neg"] = counts["delta"] < 0

        counts.loc[counts["neg"] == True, f"Δ({g1}, {g2})"] = f"< 0"
        counts.loc[counts["neg"] == False, f"Δ({g1}, {g2})"] = f"> 0"

        print(counts)
        counts["delta"] = counts["delta"].abs()
        counts["delta"] = counts["delta"].astype(int)
        print(counts)

        sns.lineplot(
            data=counts,
            x="delta",
            y="count",
            hue=f"Δ({g1}, {g2})",
            # fill=False,
            # dodge=False,
            ax=axes[1][i],
        )
        for hue_value in counts[f"Δ({g1}, {g2})"].unique():
            print(hue_value)
            subset = counts[counts[f"Δ({g1}, {g2})"] == hue_value]
            subset = subset.sort_values(by="delta")
            axes[1][i].fill_between(
                subset["delta"],
                subset["count"],
                alpha=0.5,
                # color=palette[list(df["Type"].unique()).index(hue_value)],
            )

        xlabels = []
        for t in range(0, T):
            xlabels.append(str(t) if t % 5 == 0 else "")
        xlabels.append(f">={T}")

        print(list(range(T + 1)))
        print(xlabels)

        axes[1][i].set_xticks(
            ticks=range(T + 1),
            labels=xlabels,
        )

        axes[0][i].set_title(n)

    limits = [
        min(axes[0][0].get_xlim()[0], axes[0][0].get_ylim()[0]),
        min(axes[0][0].get_xlim()[1], axes[0][0].get_ylim()[1]),
    ]
    for i, n in enumerate(Ns[1:], 1):
        new_limits = [
            min(axes[0][i].get_xlim()[0], axes[0][i].get_ylim()[0]),
            min(axes[0][i].get_xlim()[1], axes[0][i].get_ylim()[1]),
        ]
        limits = [min(limits[0], new_limits[0]), max(limits[1], new_limits[1])]
    for i, n in enumerate(Ns):
        axes[0][i].set_xlim(limits)
        axes[0][i].set_ylim(limits)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
