import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("CSV")

    parser.add_argument("-c", type=float, default=0.995)
    # parser.add_argument("-n", type=int, default=8)
    parser.add_argument("-m", type=int, default=200)
    parser.add_argument("-e", action="store_true")
    parser.add_argument("-o", type=str, default="")

    args = parser.parse_args()

    df = pd.read_csv(args.CSV)
    # df = df[df["n"] == args.n]
    df = df[df["cov"] > 0]

    print(df["cov"].gt(args.c).mean() * 100)
    print(df.describe(percentiles=[0.05, 0.1, 0.25, 0.5, 0.75]))

    df = df[df["cov"] >= args.c]

    graphs = df["tool"].unique()
    for g in graphs:
        zeros = len(df[(df["tool"] == g) & (df["nm"] == 0)])
        size = len(df[df["tool"] == g])
        print(g, zeros, size, zeros / size)
    # name1, name2, name3, name4 = graphs

    # df = df.pivot_table(index=["read", "n"], columns="tool", values="nm", aggfunc="first")
    # df = df.reset_index()
    # # df = df.dropna()
    # df = df.fillna(args.m)

    # # if not args.e:
    # #     # df = df[df[name1] != df[name2]]
    # #     df = df[df["Original"] > df["Augmented"]]

    # df = df.melt(
    #     id_vars=["read", "n"], value_vars=graphs, var_name="tool", value_name="nm"
    # )
    df["nm"].values[df["nm"].values > args.m] = args.m
    df = df[df["nm"] != args.m]

    Ns = df["n"].unique()
    Ns.sort()
    fig, axes = plt.subplots(3, len(Ns), sharex=True, sharey="row", figsize=(13, 6))
    for i, n in enumerate(Ns):
        subdf = df[df["n"] == n]
        sns.histplot(
            data=subdf, x="nm", hue="tool", element="step", discrete=True, fill=False, legend=i==0, ax=axes[0][i],
        )
        sns.histplot(
            data=subdf[subdf["type"] == "Simple"], x="nm", hue="tool", element="step", discrete=True, fill=False, legend=False, ax=axes[1][i],
        )
        sns.histplot(
            data=subdf[subdf["type"] == "Complex"], x="nm", hue="tool", element="step", discrete=True, fill=False, legend=False, ax=axes[2][i],
        )
        axes[0][i].set_title(f"n={n}")
        # if i == 0:
        #     axes[i].legend(loc="lower right")

    # T = 50
    # fig, axes = plt.subplots(1, 2)

    # ### First row: scatterplot
    # ax = axes[0]
    # ax.axline((0, 0), slope=1, color="black")

    # sns.jointplot(data=df, x=name1, y=name2, hue="class")

    # sns.scatterplot(
    #     data=df,
    #     x=name1,
    #     y=name2,
    #     hue="class",
    #     # facecolor="none",
    #     # edgecolor="blue",
    #     # alpha=0.5,
    #     s=7,
    #     # ax=ax,
    # )

    # ### Second row: barplot
    # ax = axes[1]

    # counts = df["delta"].value_counts()

    # full_range = pd.Series(range(int(min(df["delta"])), int(max(df["delta"])) + 1))
    # counts = counts.reindex(full_range, fill_value=0)
    # counts = counts.reset_index()

    # counts.columns = ["delta", "count"]

    # sum_1 = counts.loc[counts["delta"] >= T, "count"].sum()
    # sum_2 = counts.loc[counts["delta"] <= -T, "count"].sum()

    # counts = counts.drop(counts[(counts["delta"] >= T) | (counts["delta"] <= -T)].index)

    # new_rows = pd.DataFrame({"delta": [-T - 1, T + 1], "count": [sum_2, sum_1]})
    # counts = pd.concat([counts, new_rows], ignore_index=True)
    # counts = counts.drop(counts[counts["delta"] == 0].index)
    # counts["neg"] = counts["delta"] < 0

    # counts.loc[counts["neg"] == True, f"Δ({name1}, {name2})"] = f"< 0"
    # counts.loc[counts["neg"] == False, f"Δ({name1}, {name2})"] = f"> 0"

    # counts["delta"] = counts["delta"].abs()
    # counts["delta"] = counts["delta"].astype(int)

    # sns.lineplot(
    #     data=counts,
    #     x="delta",
    #     y="count",
    #     hue=f"Δ({name1}, {name2})",
    #     # fill=False,
    #     # dodge=False,
    #     ax=ax,
    # )

    # #     for hue_value in counts[f"Δ({g1}, {g2})"].unique():
    # #         subset = counts[counts[f"Δ({g1}, {g2})"] == hue_value]
    # #         subset = subset.sort_values(by="delta")
    # #         axes[1][i].fill_between(
    # #             subset["delta"],
    # #             subset["count"],
    # #             alpha=0.5,
    # #             # color=palette[list(df["Type"].unique()).index(hue_value)],
    # #         )

    # #     xlabels = []
    # #     for t in range(0, T):
    # #         xlabels.append(str(t) if t % 5 == 0 else "")
    # #     xlabels.append(f">={T}")

    # #     axes[1][i].set_xticks(
    # #         ticks=range(T + 1),
    # #         labels=xlabels,
    # #     )

    # #     axes[0][i].set_title(n)

    # # limits = [
    # #     min(axes[0][0].get_xlim()[0], axes[0][0].get_ylim()[0]),
    # #     min(axes[0][0].get_xlim()[1], axes[0][0].get_ylim()[1]),
    # # ]
    # # for i, n in enumerate(Ns[1:], 1):
    # #     new_limits = [
    # #         min(axes[0][i].get_xlim()[0], axes[0][i].get_ylim()[0]),
    # #         min(axes[0][i].get_xlim()[1], axes[0][i].get_ylim()[1]),
    # #     ]
    # #     limits = [min(limits[0], new_limits[0]), max(limits[1], new_limits[1])]
    # # for i, n in enumerate(Ns):
    # #     axes[0][i].set_xlim(limits)
    # #     axes[0][i].set_ylim(limits)

    # # plt.show()

    # # fig, ax = plt.subplots(1, 1)

    # # # x0, x1 = ax.get_xlim()
    # # # y0, y1 = ax.get_ylim()
    # # # print(x0, x1)
    # # # print(y0, y1)

    # # xs = np.linspace(m, M, 2)
    # # ax.plot(xs, xs, color="k", linestyle="--")  # y=x

    # # sns.scatterplot(data=df, x=name1, y=name2, alpha=0.3, ax=ax)
    # # # ax.set_aspect("equal", adjustable="box")  # equal data unit scaling
    # # # ax.set_box_aspect(1)

    plt.tight_layout()
    if args.o == "":
        plt.show()
    else:
        plt.savefig(args.o)


if __name__ == "__main__":
    main()
