rule dipcall:
    input:
        fa=FA,
        # bed=PARBED,
        hap1=HAP1,
        hap2=HAP2,
    output:
        vcf=pjoin(WD, SAMPLE, "dipcall", "prefix.dip.vcf.gz"),
        bed=pjoin(WD, SAMPLE, "dipcall", "prefix.dip.bed"),
        bam1=pjoin(WD, SAMPLE, "dipcall", "prefix.hap1.bam"),
        bam2=pjoin(WD, SAMPLE, "dipcall", "prefix.hap2.bam"),
    params:
        wdir=pjoin(WD, SAMPLE, "dipcall"),
        prefix=pjoin(WD, SAMPLE, "dipcall", "prefix"),
        mak=pjoin(WD, SAMPLE, "dipcall.mak"),
    threads: workflow.cores
    conda:
        "../envs/dipcall.yml"
    shell:
        """
        run-dipcall -t {threads} {params.prefix} {input.fa} {input.hap1} {input.hap2} > {params.mak}
        mkdir -p {params.wdir}
        make -j 2 -f {params.mak}
        tabix -p vcf {output.vcf}
        samtools index {output.bam1}
        samtools index {output.bam2}
        """


rule clean_dipcall:
    input:
        vcf=pjoin(WD, SAMPLE, "dipcall", "prefix.dip.vcf.gz"),
    output:
        vcf=pjoin(WD, SAMPLE, "dipcall.vcf.gz"),
    conda:
        "../envs/dipcall.yml"
    shell:
        """
        bcftools view -v indels -i '(ILEN <= -50 || ILEN >= 50)' {input.vcf} | bcftools norm -Oz --multiallelics -both > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule merge_bam:
    input:
        bam1=pjoin(WD, SAMPLE, "dipcall", "prefix.hap1.bam"),
        bam2=pjoin(WD, SAMPLE, "dipcall", "prefix.hap2.bam"),
    output:
        bam=pjoin(WD, SAMPLE, "dipcall", "prefix.dip.bam"),
    conda:
        "../envs/samtools.yml"
    shell:
        """
        samtools merge {output.bam} {input.bam1} {input.bam2}
        samtools index {output.bam}
        """


rule cutesv_asm:
    input:
        fa=FA,
        bam=pjoin(WD, SAMPLE, "dipcall", "prefix.dip.bam"),
    output:
        vcf=pjoin(WD, SAMPLE, "cutesv-asm.vcf"),
    params:
        tmp=pjoin(WD, SAMPLE, "cutesv-asm.tmp"),
    conda:
        "../envs/cutesv.yml"
    shell:
        """
        mkdir -p {params.tmp}
        cuteSV {input.bam} {input.fa} {output.vcf}.pre {params.tmp} -s 1 --genotype --report_readid -p -1 -mi 500 -md 500 --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5
        python3 $CONDA_PREFIX/lib/python3.6/site-packages/cuteSV/diploid_calling.py {output.vcf}.pre {output.vcf}
        """


rule svim_asm:
    input:
        fa=FA,
        bam1=pjoin(WD, SAMPLE, "dipcall", "prefix.hap1.bam"),
        bam2=pjoin(WD, SAMPLE, "dipcall", "prefix.hap2.bam"),
    output:
        vcf=pjoin(WD, SAMPLE, "svim-asm.vcf"),
    params:
        wd=pjoin(WD, SAMPLE, "svim-asm"),
    conda:
        "../envs/svimasm.yml"
    shell:
        """
        mkdir -p {params.wd}
        svim-asm diploid {params.wd} {input.bam1} {input.bam2} {input.fa}
        mv {params.wd}/variants.vcf {output.vcf}
        """
