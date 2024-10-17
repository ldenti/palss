rule get_truth:
    input:
        vcf=VCF,
    output:
        vcf=pjoin(WD, SAMPLE, "truth.vcf.gz"),
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        bcftools view -Oz -s {SAMPLE} -v indels -i '(ILEN <= -50 || ILEN >= 50)' {input.vcf} > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule gzip:
    input:
        "{x}.vcf",
    output:
        "{x}.vcf.gz",
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        $CONDA_PREFIX/bin/bgzip -k {input}
        tabix -p vcf {output}
        """


rule truvari_pansv:
    input:
        fa=FA,
        vcf=pjoin(WD, SAMPLE, "pansv-l{l}", "{graph}-calls.k{k}.vcf"),
        truth=pjoin(WD, SAMPLE, "truth.vcf.gz"),
        bed=lambda wildcards: TIERS[wildcards.mode],
    output:
        directory(pjoin(WD, SAMPLE, "truvari", "{mode}", "pansv-l{l}.k{k}.{graph}")),
    params:
        bed_flag=lambda wildcards: (
            "" if wildcards.mode == "full" else "--includebed " + modes[wildcards.mode]
        ),
    conda:
        "../envs/truvari.yml"
    shell:
        """
        truvari bench -b {input.truth} -c {input.vcf} -o {output} -f {input.fa} --includebed {input.bed} --passonly
        """


rule truvari_others:
    input:
        fa=FA,
        vcf=pjoin(WD, SAMPLE, "{caller}", "variations.vcf"),
        truth=pjoin(WD, SAMPLE, "truth.vcf.gz"),
        bed=lambda wildcards: TIERS[wildcards.mode],
    output:
        directory(pjoin(WD, SAMPLE, "truvari", "{mode}", "{caller}")),
    params:
        bed_flag=lambda wildcards: (
            "" if wildcards.mode == "full" else "--includebed " + modes[wildcards.mode]
        ),
    conda:
        "../envs/truvari.yml"
    shell:
        """
        truvari bench -b {input.truth} -c {input.vcf} -o {output} -f {input.fa} --includebed {input.bed} --passonly
        """


# truvari-wgt > --gtcomp
