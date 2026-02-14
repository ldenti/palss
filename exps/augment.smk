# from snakemake.utils import min_version
# min_version("6.4.1")
from os.path import join as pjoin
import random


##### config file #####
configfile: "config/config.yaml"


seed = 23
random.seed(seed)

FA = config["fa"]
GBZ = config["gbz"]
# vg gbwt --gbz-input --samples --list-names $gbz | grep -Pv "CHM13|GRCh38" > $wd/all_samples.list
SAMPLES = config["samples"]
REALFQ = config["realfq"]
# TRF = config["trf"]
WD = config["wd"]
REF = "CHM13"  # "GRCh38"

coverage = 7.5  # coverage per haplotype

#

Ns = config["ns"]  # [1, 8, 32, 64]
Ds = config["ds"]  # [0.1, 0.25, 0.5, 1.0]
Ws = [2]  # , 3]


wildcard_constraints:
    t=r"full|oneout",
    w=r"\d+",
    d=r"\d\.\d+",


#

all_samples = []
for line in open(SAMPLES):
    all_samples.append(line.strip("\n"))

random.shuffle(all_samples)

sample = all_samples[0]
kept_samples = {}
for n in Ns:
    kept_samples[str(n)] = [REF] + all_samples[1 : n + 1]
print(sample, kept_samples)


include: "./rules/prepare_data.smk"
include: "./rules/palss.smk"
include: "./rules/analyze.smk"
include: "./rules/mgc.smk"


rule run:
    input:
        # pjoin(WD, sample + "-hap1.bam"),
        # pjoin(WD, sample + "-hap2.bam"),
        # pjoin(WD, sample + "-reads.fq.gz"),
        # pjoin(WD, sample + "-reads.bam"),
        # pjoin(WD, sample + "-reads.ec.fa"),
        # pjoin(WD, sample + "-reads.ec.bam"),
        # #
        # expand(pjoin(WD, "n{n}", "pangenome-oneout.gbz"), n=Ns),
        # expand(pjoin(WD, "n{n}", "pangenome-full.gbz"), n=Ns),
        #
        expand(
            pjoin(WD, "n{n}", "palss-{t}", "specific_strings.d{d}.bam"),
            n=Ns,
            d=Ds,
            t=["full", "oneout"],
        ),
        # expand(pjoin(WD, "n{n}", "palss-{t}", "pangenome-augmented.d{d}.w{w}.gfa"), n=Ns, d=Ds, w=Ws, t=["full", "oneout"]),
        # expand(pjoin(WD, "n{n}", "pangenome-mgcactus.gfa"), n=Ns),
        # expand(
        #     pjoin(WD, "n{n}", "graphaligner", "original-{t}.gaf"),
        #     t=["full", "oneout"],
        #     n=Ns,
        # ),
        # expand(
        #     pjoin(WD, "n{n}", "graphaligner", "palss-{t}.d{d}.w{w}.gaf"),
        #     t=["full", "oneout"],
        #     n=Ns,
        #     d=Ds,
        #     w=Ws,
        # ),
        pjoin(WD, "support.csv"),
        pjoin(WD, "nm.csv"),
        pjoin(WD, "correctness.csv"),
        pjoin(WD, "anchoring.csv"),
