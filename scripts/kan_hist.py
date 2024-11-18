import sys

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt


def main():
    txt_fn = sys.argv[1]

    regions = {}
    i = 0
    for line in open(txt_fn):
        region = line.split(" ")[0]
        chrom = region.split(":")[0]
        s, e = region.split(":")[1].split("-")
        if chrom not in regions:
            regions[chrom] = []
        regions[chrom].append((int(s), int(e)))

    data = []
    nones = 0
    for chrom, kmers in regions.items():
        kmers.sort(key=lambda x: x[0])
        for (s1, e1), (s2, e2) in zip(kmers[:-1], kmers[1:]):
            if s2 < e1:
                nones += 1
            else:
                data.append([chrom, max(s2 - e1, 1)])
    print(nones, "/", nones + len(data), nones / (nones + len(data)))

    df = pd.DataFrame(data, columns=["Path", "d"])
    sns.histplot(df, x="d", hue="Path")
    plt.show()


if __name__ == "__main__":
    main()
