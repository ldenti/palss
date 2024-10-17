# TODO: FIX THIS
rule remove_sample:
    input:
        vcf=VCF,
    output:
        vcf=pjoin(WD, SAMPLE, "variations-1out.vcf.gz"),
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        $CONDA_PREFIX/bin/bcftools view -Oz --samples ^{SAMPLE} {input.vcf} > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule build_reference_graph:
    input:
        fa=FA,
    output:
        gfa=pjoin(WD, SAMPLE, "pansv-l{l}", "reference.gfa"),
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
        vcf=pjoin(WD, SAMPLE, "variations-1out.vcf.gz"),
    output:
        vg=pjoin(WD, SAMPLE, "pansv-l{l}", "variations.walts.vg"),
    threads: workflow.cores
    conda:
        "../envs/vg.yml"
    shell:
        """
        vg construct -t {threads} -r {input.fa} -v {input.vcf} --alt-paths --node-max {wildcards.l} > {output.vg}
        """


rule drop_alts:
    input:
        vg=pjoin(WD, SAMPLE, "pansv-l{l}", "variations.walts.vg"),
    output:
        vg=pjoin(WD, SAMPLE, "pansv-l{l}", "variations.vg"),
    conda:
        "../envs/vg.yml"
    shell:
        """
        vg paths --drop-paths --variant-paths -x {input.vg} | vg convert --packed-out - > {output.vg}
        """


rule build_gbwt_referencepath:
    input:
        vg=pjoin(WD, SAMPLE, "pansv-l{l}", "variations.vg"),
    output:
        gbwt=pjoin(WD, SAMPLE, "pansv-l{l}", "referencepath.gbwt"),
    threads: workflow.cores
    conda:
        "../envs/vg.yml"
    shell:
        """
        vg gbwt --index-paths --num-jobs {threads} --xg-name {input.vg} --output {output.gbwt}
        """


rule build_gbwt:
    input:
        vcf=pjoin(WD, SAMPLE, "variations-1out.vcf.gz"),
        vg=pjoin(WD, SAMPLE, "pansv-l{l}", "variations.walts.vg"),
    output:
        gbwt=pjoin(WD, SAMPLE, "pansv-l{l}", "variations.gbwt"),
        gbwtgraph=pjoin(WD, SAMPLE, "pansv-l{l}", "variations.gbwtgraph"),
    conda:
        "../envs/vg.yml"
    shell:
        """
        vg gbwt --discard-overlaps --vcf-input {input.vcf} --xg-name {input.vg} --output {output.gbwt} --graph-name {output.gbwtgraph}
        """


rule merge_gbwt:
    input:
        gbwt1=pjoin(WD, SAMPLE, "pansv-l{l}", "variations.gbwt"),
        gbwt2=pjoin(WD, SAMPLE, "pansv-l{l}", "referencepath.gbwt"),
    output:
        gbwt=pjoin(WD, SAMPLE, "pansv-l{l}", "all.gbwt"),
    conda:
        "../envs/vg.yml"
    shell:
        """
        vg gbwt --merge -o {output.gbwt} {input.gbwt1} {input.gbwt2}
        """


rule build_variation_graph_final:
    input:
        gbwt=pjoin(WD, SAMPLE, "pansv-l{l}", "all.gbwt"),
        gbwtgraph=pjoin(WD, SAMPLE, "pansv-l{l}", "variations.gbwtgraph"),
    output:
        gfa=pjoin(WD, SAMPLE, "pansv-l{l}", "variations.gfa"),
    conda:
        "../envs/vg.yml"
    shell:
        """
        vg convert --gbwt-in {input.gbwt} {input.gbwtgraph} | vg ids -s - | vg view - > {output.gfa}
        """


rule get_paths:
    input:
        gfa=pjoin(WD, SAMPLE, "pansv-l{l}", "{graph}.gfa"),
    output:
        fa=pjoin(WD, SAMPLE, "pansv-l{l}", "{graph}.paths.fa"),
    conda:
        "../envs/vg.yml"
    shell:
        """
        vg paths --extract-fasta --xg {input.gfa} > {output.fa}
        """


rule index_paths:
    input:
        fa=pjoin(WD, SAMPLE, "pansv-l{l}", "{graph}.paths.fa"),
    output:
        fmd=pjoin(WD, SAMPLE, "pansv-l{l}", "{graph}.paths.fa.fmd"),
    log:
        time=pjoin(WD, SAMPLE, "TIMES", "pansv-l{l}", "{graph}-index.time"),
    threads: workflow.cores
    shell:
        """
        ../build/rb3-prefix/src/rb3/ropebwt3 build -t {threads} -d {input.fa} > {output.fmd}
        """


rule pansv_ref:
    input:
        fq=pjoin(WD, SAMPLE, "hifi.fq"),
    output:
        fq=pjoin(WD, SAMPLE, "hifi.hifiasm.ec.fa"),
    params:
        prefix=pjoin(WD, SAMPLE, "hifi.hifiasm"),
    threads: workflow.cores
    conda:
        "../envs/hifiasm.yml"
    log:
        time=pjoin(WD, SAMPLE, "TIMES", "hifiasm.time"),
    shell:
        """
        /usr/bin/time -v hifiasm -t{threads} --write-ec --bin-only {input.fq} -o {params.prefix}
        """


rule pansv:
    input:
        gfa=pjoin(WD, SAMPLE, "pansv-l{l}", "{graph}.gfa"),
        fmd=pjoin(WD, SAMPLE, "pansv-l{l}", "{graph}.paths.fa.fmd"),
        fq=pjoin(WD, SAMPLE, "hifi.hifiasm.ec.fa"),
    output:
        txt=pjoin(WD, SAMPLE, "pansv-l{l}", "{graph}-calls.k{k}.txt"),
    params:
        k="{k}",
    log:
        time=pjoin(WD, SAMPLE, "TIMES", "pansv-l{l}", "{graph}-calling.k{k}.time"),
        log=pjoin(WD, SAMPLE, "pansv-l{l}", "{graph}-calling.k{k}.log"),
    shell:
        """
        /usr/bin/time -vo {log.time} ../pansv -k {params.k} {input.gfa} {input.fmd} {input.fq} > {output.txt} 2> {log.log}
        """


rule pansv_convert:
    input:
        gfa=pjoin(WD, SAMPLE, "pansv-l{l}", "{graph}.gfa"),
        txt=pjoin(WD, SAMPLE, "pansv-l{l}", "{graph}-calls.k{k}.txt"),
    output:
        vcf=pjoin(WD, SAMPLE, "pansv-l{l}", "{graph}-calls.k{k}.vcf"),
    log:
        time=pjoin(WD, SAMPLE, "TIMES", "pansv-l{l}", "{graph}-converting.k{k}.time"),
    shell:
        """
        /usr/bin/time -vo {log.time} python3 ../format_vcf.py {input.gfa} {input.txt} | bcftools sort > {output.vcf}
        """
