# from snakemake.utils import min_version
# min_version("6.4.1")
from os.path import join as pjoin
import random


##### config file #####
configfile: "config/config.yaml"


seed = config["seed"]
random.seed(seed)

FA = config["fa"]
GBZ = config["gbz"]
# vg gbwt --gbz-input --samples --list-names $gbz | grep -Pv "CHM13|GRCh38" > $wd/all_samples.list
SAMPLES = config["samples"]
REALFQ = config["realfq"]
# TRF = config["trf"]
WD = config["wd"]
REF = "CHM13"  # "GRCh38"

coverage = config["coverage"] / 2  # coverage per haplotype


#

Ns = config["ns"]
Ds = config["ds"]
Ws = [2]  # , 3]


wildcard_constraints:
    t=r"full|oneout",
    w=r"\d+",
    d=r"\d\.\d+",
    pv=r"0|1",


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
include: "./rules/analyses.smk"
include: "./rules/mgc.smk"


rule run:
    input:
        # contig alignment to reference
        pjoin(WD, sample + "-hap1.bam"),
        pjoin(WD, sample + "-hap2.bam"),
        # reads alignment to reference
        pjoin(WD, sample + "-reads.bam"),
        pjoin(WD, sample + "-reads.ec.bam"),
        # reads alignment to real contigs
        pjoin(WD, sample + "-reads.tohaps.bam"),
        #
        # hifiasm contigs aligned to reference
        pjoin(WD, sample + ".asm.bp.hap1.p_ctg.bam"),
        pjoin(WD, sample + ".asm.bp.hap2.p_ctg.bam"),
        # minigraph-cactus
        # expand(pjoin(WD, "n{n}", "pangenome-mgcactus.gfa"), n=Ns),
        #
        # PALSS
        expand(
            pjoin(WD, "n{n}", "palss{pv}-{t}", "pangenome-augmented.d{d}.w{w}.gfa"),
            n=Ns,
            t=["full", "oneout"],
            d=Ds,
            w=Ws,
            pv=[0, 1],
        ),
        #
        expand(
            pjoin(WD, "n{n}", "palss{pv}-{t}", "specific_strings.d{d}.bam"),
            n=Ns,
            d=Ds,
            t=["full", "oneout"],
            pv=[0, 1],
        ),
        #
        # PALSS (oneout) unanchored contigs to FULL/ONEOUT graphs
        expand(
            pjoin(
                WD,
                "n{n}",
                "palss1-oneout",
                "specific_strings.d{d}.txt.reads_with_unanchored.bp.p_ctg.to-{t}.gaf",
            ),
            n=Ns,
            t=["full", "oneout"],
            d=Ds,
            w=Ws,
        ),
        # PALSS consensus to real contigs
        expand(
            pjoin(
                WD,
                "n{n}",
                "palss{pv}-{t}",
                "resulting-consensus.d{d}.w{w}.to-contigs.bam",
            ),
            n=Ns,
            t=["full", "oneout"],
            d=Ds,
            w=Ws,
            pv=[0, 1],
        ),
        #
        # real contigs alignments
        expand(
            pjoin(WD, "n{n}", "truecontigs-aln", "original-{t}.gaf"),
            n=Ns,
            t=["full", "oneout"],
        ),
        expand(
            pjoin(WD, "n{n}", "truecontigs-aln", "palss{pv}-{t}.d{d}.w{w}.gaf"),
            n=Ns,
            t=["full", "oneout"],
            d=Ds,
            w=Ws,
            pv=[0, 1],
        ),
        # expand(pjoin(WD, "n{n}", "truecontigs-aln", "mgcactus.gaf"), n=Ns),
        pjoin(WD, sample + "-haps.50k-overlapping.bam"),
