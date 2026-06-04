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
WD = config["wd"]

Cs = [2.5, 5, 7.5, 10, 15, 30]  # coverage per haplotype
N = 1

#

samples = []
for line in open(SAMPLES):
    samples.append(line.strip("\n"))
samples = samples[:N]

#


rule run:
    input:
        pjoin(WD, "results.csv"),


rule extract_haplotype:
    input:
        gbz=GBZ,
    output:
        fa=pjoin(WD, "{sample}", "hap-{h}.fa"),
    shell:
        """
        vg paths --paths-by "{wildcards.sample}#{wildcards.h}" --extract-fasta --xg {input.gbz} > {output.fa}
        """


rule cat_haplotypes:
    input:
        fa1=pjoin(WD, "{sample}", "hap-1.fa"),
        fa2=pjoin(WD, "{sample}", "hap-2.fa"),
    output:
        fa=pjoin(WD, "{sample}", "haps.fa"),
    shell:
        """
        cat {input.fa1} {input.fa2} > {output.fa}
        """


rule pbsim3:
    input:
        fa=rules.extract_haplotype.output.fa,
        fq=REALFQ,
    output:
        fq=pjoin(WD, "{sample}", "c{c}", "pbsim3", "hap{h}_0001.fq.gz"),
    params:
        oprefix=pjoin(WD, "{sample}", "c{c}", "pbsim3", "hap{h}"),
    threads: workflow.cores / 2
    conda:
        "./envs/pbsim3.yaml"
    shell:
        """
        pbsim --id-prefix "S{wildcards.h}_" --prefix {params.oprefix} --strategy wgs --method sample --sample {input.fq} --depth {wildcards.c} --genome {input.fa} 2> {params.oprefix}.log
        """


rule combine:
    input:
        pjoin(WD, "{sample}", "c{c}", "pbsim3", "hap1_0001.fq.gz"),
        pjoin(WD, "{sample}", "c{c}", "pbsim3", "hap2_0001.fq.gz"),
    output:
        fq=pjoin(WD, "{sample}", "c{c}", "reads.fq.gz"),
    params:
        oprefix=pjoin(WD, "{sample}", "c{c}", "pbsim3"),
    threads: workflow.cores
    shell:
        """
        i=0
        for fq in {params.oprefix}/*.fq.gz
        do
            python3 ./utils/remove_n.py $fq | gzip -c > $fq.clean.gz &
            i=$((i+1))
            if (( i % {threads} == 0 )); then
                wait
            fi
        done
        wait
        #
        cat {params.oprefix}/*.clean.gz > {output.fq}
        rm {params.oprefix}/*.clean.gz
        """


rule hifiasm:
    input:
        fq=pjoin(WD, "{sample}", "c{c}", "reads.fq.gz"),
    output:
        fa1=pjoin(WD, "{sample}", "c{c}", "hifiasm.asm.bp.hap1.p_ctg.fa"),
        fa2=pjoin(WD, "{sample}", "c{c}", "hifiasm.asm.bp.hap2.p_ctg.fa"),
        fa=pjoin(WD, "{sample}", "c{c}", "hifiasm.p_ctg.fa"),
    params:
        prefix=pjoin(WD, "{sample}", "c{c}", "hifiasm.asm"),
    threads: workflow.cores
    conda:
        "./envs/hifiasm.yaml"
    shell:
        """
        hifiasm --hom-cov {wildcards.c} -o {params.prefix} -t{threads} {input.fq}
        awk '/^S/{{print ">"$2;print $3}}' {params.prefix}.bp.hap1.p_ctg.gfa > {output.fa1}
        awk '/^S/{{print ">"$2;print $3}}' {params.prefix}.bp.hap2.p_ctg.gfa > {output.fa2}
        awk '/^S/{{print ">"$2;print $3}}' {params.prefix}.bp.p_ctg.gfa > {output.fa}
        # cat {output.fa1} {output.fa2} > {output.fa}
        """


rule index_reference:
    input:
        fa=pjoin(WD, "{sample}", "haps.fa"),
    output:
        pjoin(WD, "{sample}", "haps.fa.sa"),
    conda:
        "./envs/bwa.yaml"
    threads: workflow.cores
    shell:
        """
        bwa index {input.fa}
        """


rule align_contigs_to_real_bwa:
    input:
        rfa=pjoin(WD, "{sample}", "haps.fa"),
        sa=pjoin(WD, "{sample}", "haps.fa.sa"),
        cfa=pjoin(WD, "{sample}", "c{c}", "hifiasm.p_ctg.fa"),
    output:
        bam=pjoin(WD, "{sample}", "c{c}", "hifiasm.p_ctg.bwa.bam"),
    threads: workflow.cores
    conda:
        "./envs/bwa.yaml"
    shell:
        """
        bwa mem -t{threads} {input.rfa} {input.cfa} | samtools view -bS | samtools sort > {output.bam}
        samtools index {output.bam}
        """


rule align_contigs_to_real_minimap2:
    input:
        rfa=pjoin(WD, "{sample}", "haps.fa"),
        cfa=pjoin(WD, "{sample}", "c{c}", "hifiasm.p_ctg.fa"),
    output:
        bam=pjoin(WD, "{sample}", "c{c}", "hifiasm.p_ctg.mm2.bam"),
    threads: workflow.cores
    conda:
        "./envs/minimap2.yaml"
    shell:
        """
        minimap2 -t{threads} -ax asm10 --MD --eqx -Y {input.rfa} {input.cfa} | samtools view -bS | samtools sort > {output.bam}
        samtools index {output.bam}
        """


rule summarize:
    input:
        expand(
            pjoin(WD, "{sample}", "c{c}", "hifiasm.p_ctg.{aln}.bam"),
            sample=samples,
            c=Cs,
            aln=["bwa", "mm2"],
        ),
    output:
        pjoin(WD, "results.csv"),
    conda:
        "./envs/pysam.yaml"
    shell:
        """
        python3 utils/evaluate_bams.py {WD} > {output}
        """
