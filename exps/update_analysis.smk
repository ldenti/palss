from os.path import join as pjoin
import gzip
import random

seed=23
random.seed(seed)
K = 27 # kmer size
L = 512 # max vertex size
coverage=7.5 # coverage per haplotype
n = 1 # number of runs

FA = config["fa"]
VCF = config["vcf"]
WD = config["wd"]
REALFQ = config["fq"] # for sampling based simulation

nh = 0
SAMPLES = []
for line in gzip.open(VCF, mode="rt"):
    if line.startswith("#CHROM"):
        line = line.strip("\n").split("\t")[10:] # remove hg38 from the list
        nh = 1 + 1 + (len(line)-1)*2
        random.shuffle(line)
        SAMPLES = line[:n]
        break

print(SAMPLES, nh)

rule run:
    input:
        # pjoin(WD, sample, "pangenome.gfa"),
        expand(pjoin(WD, "{sample}", "pangenome-augmented.k27.gfa"), sample=SAMPLES),

# ======================= #
# === READ SIMULATION === #
# ======================= #
rule get_sample_vcf:
    input:
        vcf=VCF,
    output:
        vcf=pjoin(WD, "{sample}", "sample.vcf.gz"),
    conda:
        "./envs/bcftools.yml"
    shell:
        """
        bcftools view -s {wildcards.sample} {input.vcf} | bcftools view -c 1 -Oz > {output.vcf}
        tabix -p vcf {output.vcf}
        """

rule get_haplotype:
    input:
        fa=FA,
        vcf=rules.get_sample_vcf.output.vcf,
    output:
        fa=pjoin(WD, "{sample}", "hap{h}.fa"),
    conda:
        "./envs/bcftools.yml"
    shell:
        """
        bcftools consensus -s {wildcards.sample} -H {wildcards.h} --fasta-ref {input.fa} {input.vcf} > {output.fa}
        samtools faidx {output.fa}
        """

rule pbsim3:
    input:
        fa=pjoin(WD, "{sample}", "hap{h}.fa"),
        fq=REALFQ,
    output:
        bam=pjoin(WD, "{sample}", "pbsim3", "hap{h}_0001.fastq"),
    params:
        oprefix=pjoin(WD, "{sample}", "pbsim3", "hap{h}"),
        cov=coverage,
    threads: workflow.cores / 2
    conda:
        "./envs/pbsim3.yml"
    shell:
        """
        pbsim --prefix {params.oprefix} --strategy wgs --method sample --sample {input.fq} --depth {params.cov} --genome {input.fa}
        """

rule combine:
    input:
        expand(pjoin(WD, "{{sample}}", "pbsim3", "hap{h}_0001.fastq"), h=[1, 2]),
    output:
        pjoin(WD, "{sample}", "hifi.fq"),
    params:
        oprefix=pjoin(WD, "{sample}", "pbsim3"),
    shell:
        """
        cat {params.oprefix}/*.fastq > {output}
        """

# ============================== #
# === PANGENOME CONSTRUCTION === #
# ============================== #
rule remove_sample:
    input:
        vcf=VCF,
        vcf2=rules.get_sample_vcf.output.vcf,
    output:
        vcf=pjoin(WD, "{sample}", "all-nosample.vcf.gz"),
    conda:
        "./envs/bcftools.yml"
    shell:
        """
        python3 scripts/remove_sample.py {input.vcf} {input.vcf2} {wildcards.sample} | bgzip -c > {output.vcf}
        sleep 1
        tabix -p vcf {output.vcf}
        """

rule vg_construct:
    input:
        fa=FA,
        vcf=rules.remove_sample.output.vcf,
    output:
        vg=pjoin(WD, "{sample}", "vg", "pangenome.walts.vg"),
    threads: workflow.cores
    conda:
        "./envs/vg.yml"
    shell:
        """
        vg construct --threads {threads} --reference {input.fa} --alt-paths --node-max {L} --vcf {input.vcf} > {output.vg}
        """

rule vg_droppaths:
    input:
        vg=rules.vg_construct.output.vg,
    output:
        vg=pjoin(WD, "{sample}", "vg", "pangenome.ref.vg"),
    threads: workflow.cores / 4
    conda:
        "./envs/vg.yml"
    shell:
        """
        vg paths --drop-paths --variant-paths -x {input.vg} > {output.vg}
        """

rule vg_gbwt:
    input:
        vcf=rules.remove_sample.output.vcf,
        vg=rules.vg_construct.output.vg,
    output:
        gbwt=pjoin(WD, "{sample}", "vg", "haplotypes.gbwt"),
    threads: workflow.cores
    conda:
        "./envs/vg.yml"
    shell:
        """
        vg gbwt --num-jobs {threads} --discard-overlaps --vcf-input {input.vcf} --xg-name {input.vg} --output {output.gbwt}
        """

rule vg_extractgam:
    input:
        gbwt=rules.vg_gbwt.output.gbwt,
        vg=rules.vg_droppaths.output.vg,
    output:
        gam=pjoin(WD, "{sample}", "vg", "haplotypes.gam"),
    threads: workflow.cores / 2
    conda:
        "./envs/vg.yml"
    shell:
        """
        vg paths --extract-gam --gbwt {input.gbwt} -x {input.vg} > {output.gam}
        """

rule vg_augment:
    input:
        vg=rules.vg_droppaths.output.vg,
        gam=rules.vg_extractgam.output.gam,
    output:
        vg=pjoin(WD, "{sample}", "pangenome.vg"),
    threads: workflow.cores / 2
    conda:
        "./envs/vg.yml"
    shell:
        """
        vg augment --label-paths {input.vg} {input.gam} > {output.vg}
        """


rule vg_view:
    input:
        vg=rules.vg_augment.output.vg,
    output:
        gfa=pjoin(WD, "{sample}", "pangenome.gfa"),
    threads: 1
    conda:
        "./envs/vg.yml"
    shell:
        """
        vg view {input.vg} > {output.gfa}
        """

# ============== #
# === METHOD === #
# ============== #

rule remove_ns:
    input:
        fq=pjoin(WD, "{sample}", "hifi.fq"),
    output:
        fq=pjoin(WD, "{sample}", "hifi.clean.fq"),
    threads: workflow.cores / 2
    conda:
        "./envs/biopython.yml"
    shell:
        """
        python3 scripts/remove_n.py {input.fq} > {output.fq}
        """

rule hifiasm:
    input:
        fq=rules.remove_ns.output.fq,
    output:
        fa=pjoin(WD, "{sample}", "hifi.ec.fa"),
    params:
        prefix=pjoin(WD, "{sample}", "hifi"),
    threads: workflow.cores
    conda:
        "./envs/hifiasm.yml"
    log:
        time=pjoin(WD, "{sample}", "times", "hifiasm.time"),
    shell:
        """
        /usr/bin/time -vo {log.time} hifiasm -t{threads} --write-ec --bin-only {input.fq} -o {params.prefix}
        """

rule get_paths:
    input:
        vg=rules.vg_augment.output.vg,
    output:
        fa=pjoin(WD, "{sample}", "pangenome.paths.fa"),
    conda:
        "./envs/vg.yml"
    shell:
        """
        vg paths --extract-fasta --xg {input.vg} > {output.fa}
        """


rule index_paths:
    input:
        fa=rules.get_paths.output.fa,
    output:
        fmd=pjoin(WD, "{sample}", "pangenome.paths.fa.fmd"),
    log:
        time=pjoin(WD, "{sample}", "times", "ropebwt3.time"),
    threads: workflow.cores
    shell:
        """
        /usr/bin/time -vo {log.time} ../build/rb3-prefix/src/rb3/ropebwt3 build -m 2G -t {threads} -d {input.fa} > {output.fmd}
        """

rule sketch:
    input:
        gfa=rules.vg_view.output.gfa,
        fmd=rules.index_paths.output.fmd,
    output:
        skt=pjoin(WD, "{sample}", "pangenome.k{k}.skt"),
    log:
        time=pjoin(WD, "{sample}", "times", "sketch-k{k}.time"),
    threads: workflow.cores
    shell:
        """
        /usr/bin/time -vo {log.time} ../pansv sketch -g{nh} -k{wildcards.k} {input.gfa} {input.fmd} > {output.skt}
        """

rule search:
    input:
        gfa=rules.vg_view.output.gfa,
        skt=rules.sketch.output.skt,
        fmd=rules.index_paths.output.fmd,
        fa=rules.hifiasm.output.fa,
    output:
        txt=pjoin(WD, "{sample}", "specific-k{k}.txt"),
    log:
        time=pjoin(WD, "{sample}", "times", "search-k{k}.time"),
    threads: workflow.cores / 2
    shell:
        """
        /usr/bin/time -vo {log.time} ../pansv search -k{wildcards.k} {input.gfa} {input.skt} {input.fmd} {input.fa} > {output.txt}
        """

rule call:
    input:
        gfa=rules.vg_view.output.gfa,
        skt=rules.sketch.output.skt,
        txt=rules.search.output.txt,
        fa=rules.hifiasm.output.fa,
    output:
        gaf=pjoin(WD, "{sample}", "specific-k{k}.gaf"),
    log:
        time=pjoin(WD, "{sample}", "times", "call-k{k}.time"),
    threads: workflow.cores / 2
    shell:
        """
        /usr/bin/time -vo {log.time} ../pansv call -k{wildcards.k} {input.gfa} {input.skt} {input.txt} {input.fa} > {output.gaf}
        """

rule augment:
    input:
        gfa=rules.vg_view.output.gfa,
        gaf=rules.call.output.gaf,
    output:
        gfa=pjoin(WD, "{sample}", "pangenome-augmented.k{k}.gfa"),
    threads: workflow.cores / 2
    conda:
        "./envs/vg.yml"
    shell:
        """
        vg augment --min-coverage 1 --gaf {input.gfa} {input.gaf} > {output.gfa}
        """
