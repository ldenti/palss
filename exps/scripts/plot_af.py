import sys
from pysam import VariantFile
import seaborn as sns
import matplotlib.pyplot as plt

vcf_fp = sys.argv[1]
sample = sys.argv[2]

data = []
nsamples = len(VariantFile(vcf_fp).header.samples) * 2

for record in VariantFile(vcf_fp):
    a1, a2 = record.samples[sample]["GT"]
    for a in set([a1, a2]):
        if a == None or a == 0:
            continue
        data.append(record.info["AC"][a - 1])  #  / record.info["AN"])

specific = len([x for x in data if x == 1])
print(specific, len(data), specific / len(data))

sns.histplot(data, bins=nsamples, discrete=True)
plt.tight_layout()
plt.show()
plt.savefig("af.png")
