### MINIMAP2 ###
rule minimap:
    input:
        fa=FA,
        fq=pjoin(WD, sample, "hifi.fq"),
    output:
        bam=pjoin(WD, sample, "hifi.bam"),
    threads: workflow.cores
    conda:
        "../envs/minimap2.yml"
    log:
        time=pjoin(WD, sample, "TIMES", "minimap2.time"),
    shell:
        """
        minimap2 -ax map-hifi --MD --eqx -Y -R '@RG\\tID:{sample}\\tSM:{sample}' -t {threads} {input.fa} {input.fq} | samtools view -bS | samtools sort -@ $(({threads}-1)) -T {output.bam} > {output.bam}
        samtools index {output.bam}
        """


### SNIFFLES2 ###
rule sniffles:
    input:
        fa=FA,
        bam=pjoin(WD, sample, "hifi.bam"),
        bed=TRF,
    output:
        vcf=pjoin(WD, sample, "sniffles", "variations.vcf"),
    params:
        support="auto",  # lambda wildcards: "2" if wildcards.support == "2" else "auto",
    conda:
        "../envs/sniffles.yml"
    threads: workflow.cores
    log:
        time=pjoin(WD, sample, "TIMES", "sniffles.time"),
    shell:
        """
        sniffles --phase --threads {threads} --reference {input.fa} --input {input.bam} --vcf {output.vcf} --tandem-repeats {input.bed} --mapq 20 --minsvlen 50 --minsupport {params.support}
        """


### CUTESV ###
rule cutesv:
    input:
        fa=FA,
        bam=pjoin(WD, sample, "hifi.bam"),
        bed=TRF,
    output:
        vcf=pjoin(WD, sample, "cutesv", "variations-w{w}.vcf"),
    params:
        tmp=pjoin(WD, sample, "cutesv", "tmp-w{w}"),
    conda:
        "../envs/cutesv.yml"
    threads: 1
    log:
        time=pjoin(WD, sample, "TIMES", "cutesv-w{w}.time"),
    shell:
        """
        mkdir -p {params.tmp}
        cuteSV --genotype --min_size 50 --min_mapq 20 --min_support {wildcards.w} --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 {input.bam} {input.fa} {output.vcf} {params.tmp}
        """
