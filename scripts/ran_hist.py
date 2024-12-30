import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

data = []
for line in open(sys.argv[1]):
    _ridx, anchors, total = line.strip("\n").split("\t")
    data.append([int(anchors), int(total)])
df = pd.DataFrame(data, columns=["Anchors", "Tot.Kmers"])

anchors_0 = len(df[df["Anchors"] == 0])
anchors_1 = len(df[df["Anchors"] == 1])
print("0 1 01 :", anchors_0, anchors_1, anchors_0 + anchors_1)
print("<10    :", len(df[df["Anchors"] < 10]))
print("Reads  :", len(df))

print(
    df.describe(
        percentiles=[
            0.1,
            0.25,
            0.4,
            0.5,
            0.6,
            0.75,
            0.90,
            0.95,
            0.96,
            0.97,
            0.98,
            0.99,
        ]
    ),
)

sns.histplot(
    df,
    #     # binrange=(1, df["d"].quantile(q=0.95)),
    #     # bins=max(1, int(df["d"].quantile(q=0.95) / 50)),
    #     legend=None,
)

plt.tight_layout()
# plt.show()
plt.savefig(sys.argv[2])
