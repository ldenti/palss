from os.path import join as pjoin

FA = config["fa"]
VCF = config["vcf"]
WD = config["wd"]

nths = workflow.cores

Ns = [4] # [1,2,4,8,16]
Ks = [27]

rule run:
    input:
        expand(pjoin(WD, "{n}", "paths-anchors.k{k}.txt"), n=Ns, k=Ks),
        expand(pjoin(WD, "{n}", "reference-anchors.k{k}.txt"), n=Ns, k=Ks),

rule select_samples:
    input:
        vcf=VCF,
    output:
        txt=pjoin(WD, "{n}", "samples.list"),
    shell:
        """
        bcftools view -h {input.vcf} | tail -1 | cut -f10- | tr "\\t" "\\n" | sort -R | head -{wildcards.n} > {output}
        """

rule extract_vcf:
    input:
        vcf=VCF,
        txt=rules.select_samples.output.txt,
    output:
        vcf=pjoin(WD, "{n}", "variations.vcf.gz"),
    shell:
        """
        bcftools view -Oz -S {input.txt} {input.vcf} > {output.vcf}
        tabix -p vcf {output.vcf}
        """

rule build_pangenome:
    input:
        fa=FA,
        vcf=rules.extract_vcf.output.vcf,
    output:
        gfa=pjoin(WD, "{n}", "pangenome.gfa"),
    log:
        time=pjoin(WD, "TIMES", "{n}", "vg-construct.time"),
    threads:
        workflow.cores
    shell:
        """
        /usr/bin/time -vo {log.time} vg construct --threads {threads} --reference {input.fa} --node-max 512 --vcf {input.vcf} | vg view - > {output.gfa}
        """

rule get_paths:
    input:
        gfa=rules.build_pangenome.output.gfa,
    output:
        fa=pjoin(WD, "{n}", "pangenome.paths.fa"),
    log:
        time=pjoin(WD, "TIMES", "{n}", "vg-paths.time"),
    shell:
        """
        /usr/bin/time -vo {log.time} vg paths -F -x {input.gfa} > {output.fa}
        """

rule rb3_index:
    input:
        fa=rules.get_paths.output.fa,
    output:
        fmd=pjoin(WD, "{n}", "pangenome.paths.fa.fmd"),
    log:
        time=pjoin(WD, "TIMES", "{n}", "rb3-build.time"),
    threads:
        workflow.cores
    shell:
        """
        /usr/bin/time -vo {log.time} ../build/rb3-prefix/src/rb3/ropebwt3 build -t {threads} -d {input.fa} > {output.fmd}
        """

rule sketch:
    input:
        gfa=rules.build_pangenome.output.gfa,
        fmd=rules.rb3_index.output.fmd,
    output:
        skt=pjoin(WD, "{n}", "pangenome-k{k}.skt"),
    log:
        time = pjoin(WD, "TIMES", "{n}", "sketch-k{k}.time"),
    threads:
        workflow.cores
    shell:
        """
        /usr/bin/time -vo {log.time} ../pansv sketch -g {wildcards.n} -k {wildcards.k} -v50000 -@{threads} {input.gfa} {input.fmd} > {output.skt}
        """

rule kan_paths:
    input:
        skt=rules.sketch.output.skt,
        fa=rules.get_paths.output.fa,
    output:
        txt=pjoin(WD, "{n}", "paths-anchors.k{k}.bed"),
    log:
        time = pjoin(WD, "TIMES", "{n}", "kan-paths-k{k}.time"),
    shell:
        """
        /usr/bin/time -vo {log.time} ../pansv kan {input.skt} {input.fa} > {output.txt}
        """

rule kan_ref:
    input:
        skt=rules.sketch.output.skt,
        fa=FA,
    output:
        bed=pjoin(WD, "{n}", "reference-anchors.k{k}.bed"),
    log:
        time = pjoin(WD, "TIMES", "{n}", "kan-reference-k{k}.time"),
    shell:
        """
        /usr/bin/time -vo {log.time} ../pansv kan {input.skt} {input.fa} > {output.bed}
        """

rule run_py:
    input:
        bed="{x}.bed",
    output:
        bed="{x}.txt",
        png="{x}.png",
    conda:
        "workflow/envs/seaborn.yml"
    shell:
        """
        python3 ../scripts/kan_hist.py {input.bed} -o {wildcards.x}
        """

