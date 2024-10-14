rule remove_sample:
    input:
        vcf=VCF,
    output:
        vcf=pjoin(WD, sample, "variations.vcf.gz"),
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        bcftools view -Oz --samples ^{sample} {input.vcf} > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule build_reference_graph:
    input:
        fa=FA,
    output:
        gfa=pjoin(WD, sample, "pansv-l{l}", "reference.gfa"),
    threads: workflow.cores
    conda:
        "../envs/vg.yml"
    shell:
        """
        vg construct -t {threads} -r {input.fa} --node-max {wildcards.l} | vg view - > {output.gfa}
        """


rule build_variation_graph:
    input:
        fa=FA,
        vcf=pjoin(WD, sample, "variations.vcf.gz"),
    output:
        vg=pjoin(WD, sample, "pansv-l{l}", "variations.walts.vg"),
    threads: workflow.cores
    conda:
        "../envs/vg.yml"
    shell:
        """
        vg construct -t {threads} -r {input.fa} -v {input.vcf} --alt-paths --node-max {wildcards.l} > {output.vg}
        """


rule drop_alts:
    input:
        vg=pjoin(WD, sample, "pansv-l{l}", "variations.walts.vg"),
    output:
        vg=pjoin(WD, sample, "pansv-l{l}", "variations.vg"),
    conda:
        "../envs/vg.yml"
    shell:
        """
        vg paths --drop-paths --variant-paths -x {input.vg} | vg convert --packed-out - > {output.vg}
        """


rule build_gbwt_referencepath:
    input:
        vg=pjoin(WD, sample, "pansv-l{l}", "variations.vg"),
    output:
        gbwt=pjoin(WD, sample, "pansv-l{l}", "referencepath.gbwt"),
    threads: workflow.cores
    conda:
        "../envs/vg.yml"
    shell:
        """
        vg gbwt --index-paths --num-jobs {threads} --xg-name {input.vg} --output {output.gbwt}
        """


rule build_gbwt:
    input:
        vcf=pjoin(WD, sample, "variations.vcf.gz"),
        vg=pjoin(WD, sample, "pansv-l{l}", "variations.walts.vg"),
    output:
        gbwt=pjoin(WD, sample, "pansv-l{l}", "variations.gbwt"),
        gbwtgraph=pjoin(WD, sample, "pansv-l{l}", "variations.gbwtgraph"),
    conda:
        "../envs/vg.yml"
    shell:
        """
        vg gbwt --discard-overlaps --vcf-input {input.vcf} --xg-name {input.vg} --output {output.gbwt} --graph-name {output.gbwtgraph}
        """


rule merge_gbwt:
    input:
        gbwt1=pjoin(WD, sample, "pansv-l{l}", "variations.gbwt"),
        gbwt2=pjoin(WD, sample, "pansv-l{l}", "referencepath.gbwt"),
    output:
        gbwt=pjoin(WD, sample, "pansv-l{l}", "all.gbwt"),
    conda:
        "../envs/vg.yml"
    shell:
        """
        vg gbwt --merge -o {output.gbwt} {input.gbwt1} {input.gbwt2}
        """


rule build_variation_graph_final:
    input:
        gbwt=pjoin(WD, sample, "pansv-l{l}", "all.gbwt"),
        gbwtgraph=pjoin(WD, sample, "pansv-l{l}", "variations.gbwtgraph"),
    output:
        gfa=pjoin(WD, sample, "pansv-l{l}", "variations.gfa"),
    conda:
        "../envs/vg.yml"
    shell:
        """
        vg convert --gbwt-in {input.gbwt} {input.gbwtgraph} | vg ids -s - | vg view - > {output.gfa}
        """


rule get_paths:
    input:
        gfa=pjoin(WD, sample, "pansv-l{l}", "{graph}.gfa"),
    output:
        fa=pjoin(WD, sample, "pansv-l{l}", "{graph}.paths.fa"),
    conda:
        "../envs/vg.yml"
    shell:
        """
        vg paths --extract-fasta --xg {input.gfa} > {output.fa}
        """


rule index_paths:
    input:
        fa=pjoin(WD, sample, "pansv-l{l}", "{graph}.paths.fa"),
    output:
        fmd=pjoin(WD, sample, "pansv-l{l}", "{graph}.paths.fa.fmd"),
    log:
        time=pjoin(WD, sample, "TIMES", "pansv-l{l}", "{graph}-index.time"),
    shell:
        """
        ../build/rb3-prefix/src/rb3/ropebwt3 build -m 2G -d {input.fa} > {output.fmd}
        """


rule pansv_ref:
    input:
        fq=pjoin(WD, sample, "hifi.fq"),
    output:
        fq=pjoin(WD, sample, "hifi.hifiasm.fq"),
    params:
        prefix=pjoin(WD, sample, "hifi.hifiasm"),
    threads: workflow.cores
    conda:
        "../envs/hifiasm.yml"
    log:
        time=pjoin(WD, sample, "TIMES", "hifiasm.time"),
    shell:
        """
        /usr/bin/time -v hifiasm -t{threads} --write-ec --bin-only {input.fq} -o {params.prefix}
        """


rule pansv:
    input:
        gfa=pjoin(WD, sample, "pansv-l{l}", "{graph}.gfa"),
        fmd=pjoin(WD, sample, "pansv-l{l}", "{graph}.paths.fa.fmd"),
        fq=pjoin(WD, sample, "hifi.hifiasm.fq"),
    output:
        txt=pjoin(WD, sample, "pansv-l{l}", "{graph}-calls.k{k}.txt"),
    params:
        k="{k}",
    log:
        time=pjoin(WD, sample, "TIMES", "pansv-l{l}", "{graph}-calling.k{k}.time"),
    shell:
        """
        /usr/bin/time -vo {log.time} ../pansv {input.gfa} {input.fmd} {input.fq} {params.k} > {output.txt}
        """
