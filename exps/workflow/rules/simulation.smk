rule get_sample_vcf:
    input:
        vcf=VCF,
    output:
        vcf=pjoin(WD, SAMPLE, "variations.vcf.gz"),
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        $CONDA_PREFIX/bin/bcftools view -s {SAMPLE} {input.vcf} | \
            $CONDA_PREFIX/bin/bcftools view -c 1 | \
            $CONDA_PREFIX/bin/bcftools +$CONDA_PREFIX/libexec/bcftools/missing2ref.so | \
            $CONDA_PREFIX/bin/bcftools +$CONDA_PREFIX/libexec/bcftools/remove-overlaps.so -Oz> {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule get_haplotype:
    input:
        fa=FA,
        vcf=pjoin(WD, SAMPLE, "variations.vcf.gz"),
    output:
        fa=pjoin(WD, SAMPLE, "hap{h}.fa"),
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        bcftools consensus -s {SAMPLE} -H {wildcards.h} --fasta-ref {input.fa} {input.vcf} > {output.fa}
        """


rule pbsim3:
    input:
        fa=pjoin(WD, SAMPLE, "hap{h}.fa"),
        fq=REALFQ,
    output:
        bam=pjoin(WD, SAMPLE, "pbsim3", "hap{h}_0001.fastq"),
    params:
        oprefix=pjoin(WD, SAMPLE, "pbsim3", "hap{h}"),
        cov=coverage,
    conda:
        "../envs/pbsim3.yml"
    shell:
        """
        pbsim --prefix {params.oprefix} --strategy wgs --method sample --sample {input.fq} --depth {params.cov} --genome {input.fa}
        """


# rule ccs:
#     input:
#         bam=pjoin(WD, SAMPLE, "reads", "hap{h}.bam"),
#     output:
#         fq=pjoin(WD, SAMPLE, "reads", "hap{h}.fq.gz"),
#     threads: workflow.cores
#     conda:
#         "envs/css.yml"
#     shell:
#         """
#         ccs --num-threads {threads} {input.bam} {output.fq}
#         """


rule combine:
    input:
        expand(pjoin(WD, SAMPLE, "pbsim3", "hap{h}_0001.fastq"), h=[1, 2]),
    output:
        pjoin(WD, SAMPLE, "hifi.fq"),
    shell:
        """
        cat {input} > {output}
        """
