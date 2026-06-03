import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns


def bar_plot(df, t, palss=False):
    Ns = df["n"].unique()
    Ns.sort()
    print(Ns)
    fig, axes = plt.subplots(2, len(Ns), sharex=True, sharey="row", figsize=(13, 7))
    for i, n in enumerate(Ns):
        for j, g in enumerate(["oneout", "full"]):
            sub_df = df[(df["graph"] == g) & (df["n"] == n)]
            subsub_df = sub_df[sub_df["supp"].isin([0, 1, 2])]
            sns.countplot(
                data=subsub_df,
                x="supp",
                hue="d" if palss else "augmentation",
                ax=axes[j][i],
                palette="bright",
            )
            # sns.move_legend(axes[j][i], "upper right")
            if j == 0:
                axes[j][i].set_title(f"{n}")
            # plt.legend(loc=3)
            # plt.tight_layout()
            # plt.show()
    plt.tight_layout()
    plt.show()
    # plt.savefig(f"support-{t}.{palss}.png")


def hist_plot(df, t, palss=False):
    Ns = df["n"].unique()
    Ns.sort()
    fig, axes = plt.subplots(2, len(Ns), sharex=True, sharey="row", figsize=(13, 7))
    for i, n in enumerate(Ns):
        for j, g in enumerate(["oneout", "full"]):
            sub_df = df[(df["graph"] == g) & (df["n"] == n)]
            sns.histplot(
                x="supp",
                hue="d" if palss else "augmentation",
                data=sub_df,
                discrete=True,
                element="step",
                # legend=None,
                fill=False,
                binrange=[0, 20],
                palette="bright",
                ax=axes[j][i],
            )
            if j == 0:
                axes[j][i].set_title(f"#samples: {n}")
                axes[j][i].set_xlabel("")
            if i == 0:
                axes[0][i].set_ylabel(f"oneout")
                axes[1][i].set_ylabel(f"full")
    plt.tight_layout()
    plt.show()
    # plt.savefig("support.png")


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("CSV")
    parser.add_argument("-p", "--palss", action="store_true")
    args = parser.parse_args()

    fn = args.CSV
    df = pd.read_csv(fn)

    if args.palss:
        df = df[df["augmentation"] == "palss"]
    else:
        df = df[(df["d"] == 0.1) | (df["d"] == -1)]

    print(df["augmentation"].unique())
    print(df.head())
    print(df.describe())

    bar_plot(df[df["kind"] == "vertex"], "vertices", args.palss)
    hist_plot(df[df["kind"] == "vertex"], "vertices", args.palss)
    bar_plot(df[df["kind"] == "edge"], "edges", args.palss)
    hist_plot(df[df["kind"] == "edge"], "edges", args.palss)

    return

    fig, axes = plt.subplots(2, len(Ns), figsize=(11, 7))
    for i, n in enumerate(Ns):
        for j, g in enumerate(["oneout", "full"]):
            sub_df = df[(df["graph"] == g) & (df["n"] == n)]
            sns.boxplot(
                data=sub_df,
                x="augmentation",
                y="l",
                log_scale=True,
                ax=axes[j][i],
            )
            # sns.move_legend(axes[j][i], "upper right")
            if j == 0:
                axes[j][i].set_title(f"{n}")
            # plt.legend(loc=3)
            # plt.tight_layout()
            # plt.show()
    plt.tight_layout()
    plt.show()
    # plt.savefig("support.png")

    fig, axes = plt.subplots(2, len(Ns), figsize=(11, 7))
    for i, n in enumerate(Ns):
        for j, g in enumerate(["oneout", "full"]):
            sub_df = df[(df["graph"] == g) & (df["n"] == n)]
            sub_df = sub_df[sub_df["supp"].isin([0, 1, 2])]
            sns.boxplot(
                data=sub_df,
                x="augmentation",
                y="l",
                hue="supp",
                log_scale=True,
                ax=axes[j][i],
            )
            if j == 0:
                axes[j][i].set_title(f"{n}")
    plt.tight_layout()
    plt.show()
    # plt.savefig("support.png")


if __name__ == "__main__":
    main()
