rule get_haplotype:
    input:
        fa=FA,
        vcf=VCF,
    output:
        fa=pjoin(WD, sample, "hap{h}.fa"),
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        bcftools consensus -s {sample} -H {wildcards.h} --fasta-ref {input.fa} {input.vcf} > {output.fa}
        """


rule pbsim3:
    input:
        fa=pjoin(WD, sample, "hap{h}.fa"),
        fq=REALFQ,
    output:
        bam=pjoin(WD, sample, "pbsim3", "hap{h}_0001.fastq"),
    params:
        oprefix=pjoin(WD, sample, "pbsim3", "hap{h}"),
        cov=coverage,
    conda:
        "../envs/pbsim3.yml"
    shell:
        """
        pbsim --prefix {params.oprefix} --strategy wgs --method sample --sample {input.fq} --depth {params.cov} --genome {input.fa}
        """


# rule ccs:
#     input:
#         bam=pjoin(WD, sample, "reads", "hap{h}.bam"),
#     output:
#         fq=pjoin(WD, sample, "reads", "hap{h}.fq.gz"),
#     threads: workflow.cores
#     conda:
#         "envs/css.yml"
#     shell:
#         """
#         ccs --num-threads {threads} {input.bam} {output.fq}
#         """


rule combine:
    input:
        expand(pjoin(WD, sample, "pbsim3", "hap{h}_0001.fastq"), h=[1, 2]),
    output:
        pjoin(WD, sample, "hifi.fq"),
    shell:
        """
        cat {input} > {output}
        """
