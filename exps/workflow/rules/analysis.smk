rule get_truth:
    input:
        vcf=pjoin(WD, SAMPLE, "variations.vcf.gz"),
    output:
        vcf=pjoin(WD, SAMPLE, "truth.vcf.gz"),
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        bcftools view -Oz -v indels -i '(ILEN <= -50 || ILEN >= 50)' {input.vcf} > {output.vcf}
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
        vcf=pjoin(WD, SAMPLE, "pansv-l{l}", "{graph}-calls.k{k}.vcf.gz"),
        truth=pjoin(WD, SAMPLE, "truth.vcf.gz"),
        bed=lambda wildcards: TIERS[wildcards.mode],
    output:
        directory(pjoin(WD, SAMPLE, "truvari", "{mode}", "pansv-l{l}.k{k}.{graph}")),
    params:
        bed_flag=lambda wildcards: (
            "" if wildcards.mode == "full" else "--includebed " + modes[wildcards.mode]
        ),
        tmp=pjoin(WD, SAMPLE, "truvari", "{mode}", "pansv-l{l}.k{k}.{graph}.tmp"),
    conda:
        "../envs/truvari.yml"
    shell:
        """
        mkdir -p {params.tmp}
        export TMPDIR="{params.tmp}"
        truvari bench -b {input.truth} -c {input.vcf} -o {output} --passonly {params.bed_flag}
        rm -rf {params.tmp}
        """


rule truvari_others:
    input:
        fa=FA,
        vcf=pjoin(WD, SAMPLE, "{caller}.vcf.gz"),
        truth=pjoin(WD, SAMPLE, "truth.vcf.gz"),
        bed=lambda wildcards: TIERS[wildcards.mode],
    output:
        directory(pjoin(WD, SAMPLE, "truvari", "{mode}", "{caller}")),
    params:
        bed_flag=lambda wildcards: (
            "" if wildcards.mode == "full" else "--includebed " + modes[wildcards.mode]
        ),
        tmp=pjoin(WD, SAMPLE, "truvari", "{mode}", "{caller}.tmp")
    conda:
        "../envs/truvari.yml"
    shell:
        """
        mkdir -p {params.tmp}
        export TMPDIR="{params.tmp}"
        truvari bench -b {input.truth} -c {input.vcf} -o {output} --passonly {params.bed_flag}
        rm -rf {params.tmp}
        """


# truvari-wgt > --gtcomp
