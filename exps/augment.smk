# from snakemake.utils import min_version
# min_version("6.4.1")
from os.path import join as pjoin
import random


##### config file #####
# configfile: "config/config.yaml"


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

coverage = config["coverage"] / 2  # coverage per haplotype

# tools
cactus_activate = config["tools"]["cactus"]
HIFIASM_BIN = config["tools"]["hifiasm-custom"]

#

Ns = config["ns"]
Ds = config["ds"]
Ws = [2]  # , 3]
Ms = config["ms"]
Cs = config["cs"]

#

graphs = ["oneout", "full"]
SIZES = [10000, 50000]


wildcard_constraints:
    h=r"1|2",
    graph=r"full|oneout",
    w=r"\d+",
    d=r"\d\.\d+",
    m=r"\d+",
    c=r"\d\.\d+",
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


rule run:
    input:
        # prepare_data.smk
        expand(pjoin(WD, "n{n}", "pangenome-{graph}.gbz"), n=Ns, graph=graphs),
        pjoin(WD, sample + "-reads.fq.gz"),
        pjoin(WD, sample + "-reads.ec.fa"),
        #
        # minigraph-cactus
        expand(pjoin(WD, "n{n}", "pangenome-mgcactus.gfa"), n=Ns),
        #
        #
        #
        # # contig alignment to reference
        # pjoin(WD, sample + "-hap1.bam"),
        # pjoin(WD, sample + "-hap2.bam"),
        # reads alignment to reference
        pjoin(WD, sample + "-reads.bam"),
        pjoin(WD, sample + "-reads.ec.bam"),
        # # reads alignment to real contigs
        # pjoin(WD, sample + "-reads.tohaps.bam"),
        # #
        # # hifiasm contigs aligned to reference
        # pjoin(WD, sample + ".asm.bp.hap1.p_ctg.bam"),
        # pjoin(WD, sample + ".asm.bp.hap2.p_ctg.bam"),
        # #
        # PALSS
        expand(
            pjoin(
                WD,
                "n{n}",
                "palss-{graph}",
                "pangenome-augmented.d{d}.w{w}.c{c}.m{m}.gfa",
            ),
            n=Ns,
            graph=["oneout"],
            d=Ds,
            w=Ws,
            c=Cs,
            m=Ms,
        ),
        expand(
            pjoin(WD, "n{n}", "palss-{graph}", "pangenome-augmented.d{d}.w{w}.gfa"),
            n=Ns,
            graph=["oneout"],
            d=Ds,
            w=Ws,
        ),
        #
        # SFS aligned to reference genome
        expand(
            pjoin(WD, "n{n}", "palss-{graph}", "specific_strings.d{d}.bam"),
            n=Ns,
            d=Ds,
            graph=["oneout"],
        ),
        # PALSS consensus to real contigs
        expand(
            pjoin(
                WD,
                "n{n}",
                "palss-{graph}",
                "anchored-consensus.d{d}.w{w}.to-contigs.bam",
            ),
            n=Ns,
            graph=["oneout"],
            d=Ds,
            w=Ws,
        ),
        expand(
            pjoin(
                WD,
                "n{n}",
                "palss-{graph}",
                "unanchored-consensus.d{d}.w{w}.c{c}.m{m}.to-contigs.bam",
            ),
            n=Ns,
            graph=["oneout"],
            d=Ds,
            w=Ws,
            c=Cs,
            m=Ms,
        ),


### NM (AKA RECALL)
rule get_nm:
    input:
        expand(
            pjoin(
                WD,
                "n{n}",
                "truecontigs-aln",
                "palss-{graph}.d{d}.w{w}.{size}.gaf",
            ),
            n=Ns,
            graph=["oneout"],
            d=Ds,
            w=Ws,
            size=SIZES,
        ),
        expand(
            pjoin(
                WD,
                "n{n}",
                "truecontigs-aln",
                "palss-{graph}.d{d}.w{w}.c{c}.m{m}.{size}.gaf",
            ),
            n=Ns,
            graph=["oneout"],
            d=Ds,
            w=Ws,
            c=Cs,
            m=Ms,
            size=SIZES,
        ),
        expand(
            pjoin(WD, "n{n}", "truecontigs-aln", "original-{graph}.{size}.gaf"),
            n=Ns,
            graph=["full", "oneout"],
            size=SIZES,
        ),
        expand(
            pjoin(WD, "n{n}", "truecontigs-aln", "mgcactus.{size}.gaf"),
            n=Ns,
            size=SIZES,
        ),
        expand(pjoin(WD, sample + "-haps.{size}-overlapping.bam"), size=SIZES),
    output:
        pjoin(WD, "nm.csv"),
    conda:
        "./envs/pysam.yaml"
    shell:
        """
        python3 ./utils/get_nm.py {WD} {sample} > {output}
        """


### SUPPORT (AKA PRECISION)
rule get_support:
    input:
        expand(
            pjoin(WD, "n{n}", "tables", "{graph}.{size}.support.csv"),
            n=Ns,
            graph=["full", "oneout"],
            size=SIZES,
        ),
        expand(
            pjoin(WD, "n{n}", "tables", "mgc.{size}.support.csv"),
            n=Ns,
            size=SIZES,
        ),
        #
        expand(
            pjoin(
                WD,
                "n{n}",
                "tables",
                "palss-{graph}.d{d}.w{w}.{size}.support.csv",
            ),
            n=Ns,
            graph=["oneout"],
            d=Ds,
            w=Ws,
            size=SIZES,
        ),
        expand(
            pjoin(
                WD,
                "n{n}",
                "tables",
                "palss-{graph}.d{d}.w{w}.c{c}.m{m}.{size}.support.csv",
            ),
            n=Ns,
            graph=["oneout"],
            d=Ds,
            w=Ws,
            c=Cs,
            m=Ms,
            size=SIZES,
        ),
    output:
        pjoin(WD, "support.csv"),
    shell:
        """
        head -1 {input[0]} > {output}
        for i in {input} ; do sed '1d' $i ; done >> {output}
        """


rule get_support_original:
    input:
        gfa=pjoin(WD, "n{n}", "pangenome-{graph}.gfa"),
        gaf=pjoin(WD, "n{n}", "truecontigs-aln", "original-{graph}.{size}.gaf"),
        bed=BED,
    output:
        csv=pjoin(WD, "n{n}", "tables", "{graph}.{size}.support.csv"),
    conda:
        "./envs/intervaltree.yaml"
    threads: workflow.cores / 4
    shell:
        """
        python3 ./utils/get_support.py -s {sample} -n {wildcards.n} {input.gfa} {input.gaf} > {output.csv}.unflagged
        python3 ./utils/flag_vertices.py {input.gfa} {output.csv}.unflagged {input.bed} > {output.csv}
        """


rule get_support_mgc:
    input:
        gfa=pjoin(WD, "n{n}", "pangenome-mgcactus.gfa"),
        gaf=pjoin(WD, "n{n}", "truecontigs-aln", "mgcactus.{size}.gaf"),
        bed=BED,
    output:
        csv=pjoin(WD, "n{n}", "tables", "mgc.{size}.support.csv"),
    conda:
        "./envs/intervaltree.yaml"
    threads: workflow.cores / 4
    shell:
        """
        python3 ./utils/get_support.py -s {sample} -n {wildcards.n} {input.gfa} {input.gaf} > {output.csv}.unflagged
        python3 ./utils/flag_vertices.py {input.gfa} {output.csv}.unflagged {input.bed} > {output.csv}
        """


rule get_support_palss:
    input:
        gfa=pjoin(
            WD,
            "n{n}",
            "palss-{graph}",
            "pangenome-augmented.d{d}.w{w}.gfa",
        ),
        gaf=pjoin(
            WD,
            "n{n}",
            "truecontigs-aln",
            "palss-{graph}.d{d}.w{w}.{size}.gaf",
        ),
        bed=BED,
    output:
        csv=pjoin(
            WD,
            "n{n}",
            "tables",
            "palss-{graph}.d{d}.w{w}.{size}.support.csv",
        ),
    conda:
        "./envs/intervaltree.yaml"
    threads: workflow.cores / 4
    shell:
        """
        python3 ./utils/get_support.py -n {wildcards.n} {input.gfa} {input.gaf} > {output.csv}.unflagged
        python3 ./utils/flag_vertices.py {input.gfa} {output.csv}.unflagged {input.bed} > {output.csv}
        """


rule get_support_palss_refine:
    input:
        gfa=pjoin(
            WD,
            "n{n}",
            "palss-{graph}",
            "pangenome-augmented.d{d}.w{w}.c{c}.m{m}.gfa",
        ),
        gaf=pjoin(
            WD,
            "n{n}",
            "truecontigs-aln",
            "palss-{graph}.d{d}.w{w}.c{c}.m{m}.{size}.gaf",
        ),
        bed=BED,
    output:
        csv=pjoin(
            WD,
            "n{n}",
            "tables",
            "palss-{graph}.d{d}.w{w}.c{c}.m{m}.{size}.support.csv",
        ),
    conda:
        "./envs/intervaltree.yaml"
    threads: workflow.cores / 4
    shell:
        """
        python3 ./utils/get_support.py -n {wildcards.n} {input.gfa} {input.gaf} > {output.csv}.unflagged
        python3 ./utils/flag_vertices.py {input.gfa} {output.csv}.unflagged {input.bed} > {output.csv}
        """
