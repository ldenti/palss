import sys

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch
from matplotlib.cbook import boxplot_stats


def main():
    gaf_fn = sys.argv[1]
    data = {}
    for line in open(gaf_fn):
        line = line.strip("\n").split("\t")
        vidx, support = line[0].split(".")
        support = int(support)
        nm = int(line[12].split(":")[-1])
        if vidx not in data:
            data[vidx] = [support, nm]
        else:
            data[vidx] = [support, min(data[vidx][1], nm)]

    df = pd.DataFrame(data.values(), columns=["Support", "NM"])

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1)

    sns.boxplot(data=df, x="Support", y="NM", showfliers=True, ax=ax2)
    if ax2.get_ylim()[1] > 26:
        ax2.set_ylim(-1, 26)

    # get bars manually
    d = {}
    nm0 = 0
    for s, nm in data.values():
        d[s] = d[s] + 1 if s in d else 1
        nm0 += nm == 0
    print(nm0, len(data), nm0 / len(data))

    outliers = {}
    for s in d:
        outliers[s] = len(boxplot_stats(df[df["Support"] == s]["NM"])[0]["fliers"])

    df2 = []
    for k, v in d.items():
        df2.append([k, v, outliers[k]])
    df2 = pd.DataFrame(df2, columns=["Support", "Count", "Outliers"])

    sns.barplot(data=df2, x="Support", y="Count", ax=ax1)
    sns.barplot(data=df2, x="Support", y="Outliers", ax=ax1)

    cmap = sns.color_palette()
    legend = [
        Patch(facecolor=cmap[0], edgecolor="w", label="Non-Outliers"),
        Patch(facecolor=cmap[1], edgecolor="w", label="Outliers"),
    ]
    ax1.legend(handles=legend)

    sns.histplot(data=df, x="NM", ax=ax3, discrete=True, binrange=[0, 25])
    plt.tight_layout()
    plt.show()

    # plt.close()
    # sns.boxplot(data=df, x="Support", y="NM", showfliers=True)
    # plt.show()


if __name__ == "__main__":
    main()
