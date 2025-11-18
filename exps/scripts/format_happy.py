import sys
import os
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def parse_csv(csv_fn):
    data = []
    for line in open(csv_fn):
        if line[0] == "T":
            continue
        line = line.split(",")
        ty = line[0]
        flt = line[1]
        tp = int(line[3])
        fn = int(line[4])
        fp = int(line[6])
        r = round(float(line[10]), 3)
        p = round(float(line[11]), 3)
        f1 = round(float(line[13]) if line[13] != "" else 0.0, 3)
        data.append([ty, flt, tp, fp, fn, p, r, f1])
    return data


def main():
    wd = sys.argv[1]

    data = []
    for csv_fn in glob.glob(os.path.join(wd, "*", "*.summary.csv")):
        wbed = True
        if "nobed" in csv_fn:
            wbed = False
        rows = parse_csv(csv_fn)
        t1, _, t2 = csv_fn.split("/")[-2].split("-")
        for row in rows:
            data.append([t2, t1, wbed] + row)

    df = pd.DataFrame(
        data,
        columns=[
            "Truth",
            "Caller",
            "Confident",
            "Type",
            "Filter",
            "TP",
            "FP",
            "FN",
            "P",
            "R",
            "F1",
        ],
    )
    df.to_csv(sys.stdout, index=False)

    df = df[
        (df["Confident"] == True) & (df["Type"] == "SNP") & (df["Filter"] == "PASS")
    ]

    sns.scatterplot(df, x="P", y="R", hue="Caller", style="Truth")
    plt.show()


if __name__ == "__main__":
    main()
