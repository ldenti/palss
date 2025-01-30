import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

txt = sys.argv[1]
fai = sys.argv[2]
# l = int(sys.argv[2])

data = {}
for line in open(fai):
    ridx, l = line.split("\t")[:2]
    data[ridx] = [0, int(l)]

for line in open(txt):
    ridx, anchors, total = line.strip("\n").split("\t")
    anchors = int(anchors)
    total = int(total)
    # if total > l:
    #     data.append([anchors, total])
    #     if anchors < 2:
    #         print(_ridx)
    data[ridx][0] = anchors
    # data[ridx][1] = total

df = pd.DataFrame(list(data.values()), columns=["Anchors", "ReadLen"])

anchors_0 = len(df[df["Anchors"] == 0])
anchors_1 = len(df[df["Anchors"] == 1])
print("0 1 01 :", anchors_0, anchors_1, anchors_0 + anchors_1)
print("<10    :", len(df[df["Anchors"] < 10]))
print("Reads  :", len(df))

for x in np.arange(0, 1.01, 0.01):
    print(x, df["Anchors"].quantile(x))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 6), dpi=80)
sns.histplot(
    df,
    # binrange=df["Tot.Kmers"].quantile(q=1)),
    # bins=bins=max(1, int(df["Tot.Kmers"].quantile(q=1) / 50)),
    # legend=None,
    ax=ax1,
)

percs = np.linspace(0, 100, 1000)
qn_a = np.percentile(df["Anchors"], percs)
qn_b = np.percentile(df["ReadLen"], percs)

ax2.plot(qn_a, qn_b, ls="", marker="o", color="g", alpha=0.71)
x = np.linspace(np.min((qn_a.min(), qn_b.min())), np.max((qn_a.max(), qn_b.max())))
ax2.plot(x, x, color="k", ls="--")

# ax1.set_xlabel("#Anchors")
ax2.set_xlabel("#Anchors")
ax2.set_ylabel("ReadLen")

plt.tight_layout()
plt.show()
# plt.savefig(sys.argv[3])
