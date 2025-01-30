import sys
import pandas as pd
import seaborn as sns
import matplotlib
from matplotlib import pyplot as plt
import pickle
import glob
import os
import numpy as np

sns.set(style="whitegrid")


def main():
    wd = sys.argv[1]
    data = []
    runs = set()
    for f in glob.glob(os.path.join(wd, "*", "*.pkl")):
        n, fname = f.split("/")[-2:]
        n = int(n)
        t = os.path.splitext(fname)[0]
        if t.endswith("w3") or t.endswith("w5"):
            continue

        if "augmented" not in t:
            if "1out" in t:
                t = "1OUT"
            else:
                t = "FULL"
        else:
            if "1out" in t:
                t = "1OUT-AUG"
            else:
                t = "FULL-AUG"
        runs.add(t)

        fb = open(f, "rb")
        _ = pickle.load(fb)
        NM = pickle.load(fb)
        fb.close()

        for nm in NM:
            # if nm > 500:
            #     continue
            data.append([n, t, nm])

    df = pd.DataFrame(data, columns=["G", "Graph", "NM"])

    plt.figure(figsize=(8, 6), dpi=80)
    x = sns.boxplot(
        data=df,
        x="G",
        y="NM",
        hue="Graph",
        hue_order=["1OUT-AUG", "1OUT", "FULL", "FULL-AUG"],
        showfliers=False,
    )
    plt.xlabel("#Samples")
    yl = x.get_ylim()
    plt.yticks(np.arange(-5, int(yl[1]) + (int(yl[1]) % 5), 5))
    plt.tight_layout()

    plt.show()
    plt.savefig("nm.png")


if __name__ == "__main__":
    main()
