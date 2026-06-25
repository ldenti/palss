from os.path import join as pjoin
import random


seed = config["seed"]
random.seed(seed)

FA = config["fa"]
GBZ = config["gbz"]
# vg gbwt --gbz-input --samples --list-names $gbz | grep -Pv "CHM13|GRCh38" > $wd/all_samples.list
SAMPLES = config["samples"]
REALFQ = config["realfq"]
# TRF = config["trf"]
WD = config["wd"]
REF = config["ref"]
BED = config["bed"]

coverages = config["coverage"]  # coverage per haplotype

# tools
cactus_activate = config["tools"]["cactus"]
HIFIASM_BIN = config["tools"]["hifiasm-custom"]

Ns = config["ns"]

#

Ds = config["ds"]
Ws = config["ws"]

#

graphs = ["oneout", "full"]
# SIZES = [10000, 50000]
SIZES = [50000]


wildcard_constraints:
    h=r"1|2",
    graph=r"full|oneout",
    cov=r"\d+",
    n=r"\d+",
    w=r"\d+",
    d=r"\d\.\d+",
    size=r"\d+",


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
include: "./rules/mgc.smk"
include: "./rules/palss.smk"
include: "./rules/analyses.smk"
#
include: "./rules/support.smk"
include: "./rules/nm.smk"
#
include: "./rules/downstream.smk"


rule run:
    input:
        # prepare_data.smk
        expand(pjoin(WD, "n{n}", "pangenome-{graph}.gbz"), n=Ns, graph=graphs),
        expand(pjoin(WD, sample + "-cov{cov}.fq.gz"), cov=coverages),
        expand(pjoin(WD, sample + "-cov{cov}.ec.fa"), cov=coverages),
        #
        # minigraph-cactus
        expand(
            pjoin(WD, "mgcactus", "n{n}", "cov{cov}", "pangenome-mgcactus.gfa"),
            n=Ns,
            cov=coverages,
        ),
        # 
        # contig alignment to reference
        pjoin(WD, sample + "-hap1.bam"),
        pjoin(WD, sample + "-hap2.bam"),
        # reads alignment to reference
        expand(pjoin(WD, sample + "-cov{cov}.bam"), cov=coverages),
        expand(pjoin(WD, sample + "-cov{cov}.ec.bam"), cov=coverages),
        # reads alignment to real contigs
        # pjoin(WD, sample + "-reads.tohaps.bam"),
        #
        # # hifiasm contigs aligned to reference
        # pjoin(WD, sample + ".asm.bp.hap1.p_ctg.bam"),
        # pjoin(WD, sample + ".asm.bp.hap2.p_ctg.bam"),
        #
        # PALSS
        expand(
            pjoin(WD, "palss", "n{n}", "cov{cov}", "augmented-{graph}.d{d}.w{w}.gfa"),
            n=Ns,
            cov=coverages,
            graph=["oneout"],
            d=Ds,
            w=Ws,
        ),
        #
        # SFS aligned to reference genome
        # expand(
        #     pjoin(WD, "n{n}", "palss-{graph}", "specific_strings.d{d}.bam"),
        #     n=Ns,
        #     d=Ds,
        #     graph=["oneout"],
        # ),
        #
        # PALSS consensus to real contigs
        expand(
            pjoin(
                WD,
                "palss",
                "n{n}",
                "cov{cov}",
                "anchoredconsensus-{graph}.d{d}.w{w}.to-contigs.bam",
            ),
            n=Ns,
            cov=coverages,
            graph=["oneout"],
            d=Ds,
            w=Ws,
        ),
        #
        pjoin(WD, "nm.csv"),
        pjoin(WD, "support.csv"),
        pjoin(WD, "downstream.csv"),
