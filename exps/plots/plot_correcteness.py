import sys

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt


def main():
    csv_fn = sys.argv[1]
    df = pd.read_csv(csv_fn)

    # df = df[df["n"] == 32]
    # df = df[df["Graph"] == "full"]

    df["Perfect"] = df["NM"] == 0

    # df = df[df["Perfect"] == False]

    df["Errors"] = df["NM"]
    df["Added"] = df["NB"]  # df["NM"] + df["NB"]

    print(df.head())
    print(df.size)

    # nerrors = df["NM"].sum()
    # added = sum([len(x) for x in newv.values()])
    # # print(min(errors), max(errors))
    # print(nerrors, added, round(nerrors / added * 100, 3))

    # print(df["Perfect"].value_counts())

    # sns.boxplot(data=df, x="Length", hue="Perfect", log_scale=True)
    # plt.tight_layout()
    # # plt.show()
    # plt.savefig(f"{xxx}-1.png")
    # plt.close()

    # subdf.sort_values("Added")

    ndf = pd.pivot_table(
        df, values=["Errors", "Added"], index=["n", "Graph", "d"], aggfunc="sum"
    ).reset_index()
    ndf["Ratio"] = (ndf["Errors"] / ndf["Added"]).round(3)
    # df.columns.name = None
    # df = df.rename(columns={0: "0", 1: "1"})
    print(ndf)

    Ns = ndf["n"].unique()
    fig, axes = plt.subplots(2, len(Ns), sharex=True, sharey="row", figsize=(13, 7))
    for i, n in enumerate(Ns):
        for j, graph in enumerate(ndf["Graph"].unique()):
            print(n, graph)
            sdf = ndf[(ndf["n"] == n) & (ndf["Graph"] == graph)]
            sns.barplot(
                data=sdf, x="d", y="Added", color="black", ax=axes[j][i], label="Added"
            )
            for x, container in enumerate(axes[j][i].containers):
                axes[j][i].bar_label(container, labels=sdf["Ratio"])
            sns.barplot(
                data=sdf, x="d", y="Errors", color="red", ax=axes[j][i], label="Errors"
            )

            # for p in x.patches:
            #     print(sdf["Ratio"].iloc[int(p.get_x() // p.get_width())])
            #     axes[j][i].annotate(
            #         p.get_height(),
            #         (p.get_x() + p.get_width() / 2.0, p.get_height()),
            #         ha="center",
            #         va="bottom",
            #         fontsize=12,
            #         color="black",
            #         text=sdf["Ratio"].iloc[int(p.get_x() // p.get_width())],
            #     )

            if j == 0:
                axes[j][i].set_title(n)

    # sorted_x = df.sort_values("Added")["Consensus"]

    # sns.barplot(
    #     data=df,
    #     x="Consensus",
    #     y="Added",
    #     hue="d",
    #     order=sorted_x,
    #     palette="dark",
    # )  # , label="Added")

    # sns.barplot(
    #     data=df,
    #     x="Consensus",
    #     y="Errors",
    #     hue="d",
    #     order=sorted_x,
    #     palette="pastel",
    # )

    # sns.barplot(data=df, x="Consensus", y="Errors", label="Errors")
    # plt.xticks(plt.xticks(), [])
    # plt.ylabel("Count")
    plt.tight_layout()
    plt.show()
    # plt.savefig(f"{xxx}-2.png")


if __name__ == "__main__":
    main()
