import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns

sns.set_theme()
# sns.set_style("whitegrid")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("RT")
    parser.add_argument("-l", default=50000, type=int)
    parser.add_argument("--nm", default=100, type=int)
    parser.add_argument("-c", default=0.995, type=float)
    parser.add_argument("--steps", default=5, type=int)
    parser.add_argument("--title", default="", type=str)
    parser.add_argument("-o", default="", type=str)
    args = parser.parse_args()

    df = pd.read_csv(args.RT)
    df = df[df["clen"] == args.l]
    df = df[df["tool"] == "original-full"]
    df = df[df["coverage"] == 15]
    df = df[df["cov"] >= args.c]

    grouped_df = df.groupby(["n", "type"])[["nm"]]

    data = []
    for name, grp in grouped_df:
        n, t = name
        for nm in range(0, args.nm + 1, args.steps):
            R = grp[grp["nm"] <= nm].count().iloc[0] / grp.count().iloc[0]
            data.append([n, t, nm, round(R * 100, 2)])
    df = pd.DataFrame(data, columns=["n", "type", "nm", "R"])
    print(df)

    Ns = df["n"].unique()
    Ns.sort()
    nn = len(Ns) if len(Ns) > 1 else 2
    fig, axes = plt.subplots(1, nn, sharex=True, sharey=True, figsize=(11, 4))
    for col, n in enumerate(Ns):
        sns.lineplot(
            data=df[df["n"] == n],
            x="nm",
            y="R",
            style="type",
            legend=col == nn - 1,
            ax=axes[col],
        )

        axes[col].set_title(f"n={n}")
        axes[col].set_xlabel(f"NM ≤ x")

        if col == nn - 1:
            axes[col].legend(loc="lower right")
        if col == 0:
            axes[col].set_ylabel(f"% alignments")

    if args.title != "":
        plt.suptitle(args.title)
    plt.tight_layout()
    if args.o == "":
        plt.show()
    else:
        plt.savefig(args.o)


if __name__ == "__main__":
    main()
