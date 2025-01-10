from os.path import join as pjoin

FA = config["fa"]
GFA = config["gfa"]
FQ = config["fq"]
WD = config["wd"]
n = config["nh"]

nths = workflow.cores

Ks = [23, 27, 31]

rule run:
    input:
        expand(pjoin(WD, "k{k}", "reference.kan.png"), k=Ks),
        expand(pjoin(WD, "k{k}", "sample.nuk.png"), k=Ks),

rule faidx:
    input:
        FA,
    output:
        FA + ".fai",
    shell:
        """
        samtools faidx {input}
        """

rule get_paths:
    input:
        gfa=GFA,
    output:
        fa=pjoin(WD, "paths.fa"),
    shell:
        """
        vg paths --extract-fasta --xg {input.gfa} > {output.fa}
        """

rule rb3_index:
    input:
        fa=rules.get_paths.output.fa,
    output:
        fmd=pjoin(WD, "pangenome.paths.fa.fmd"),
    log:
        time=pjoin(WD, "TIMES", "rb3-build.time"),
    threads:
        workflow.cores
    shell:
        """
        /usr/bin/time -vo {log.time} ../build/rb3-prefix/src/rb3/ropebwt3 build -t {threads} -d {input.fa} > {output.fmd}
        """

rule sketch:
    input:
        gfa=GFA,
        fmd=rules.rb3_index.output.fmd,
    output:
        skt=pjoin(WD, "k{k}", "pangenome.skt"),
    #params:
    #    n = n;
    log:
        time = pjoin(WD, "TIMES", "k{k}", "sketch.time"),
    threads:
        workflow.cores
    shell:
        """
        /usr/bin/time -vo {log.time} ../pansv sketch -g {n} -k {wildcards.k} -v50000 -@{threads} {input.gfa} {input.fmd} > {output.skt}
        """

rule kan_ref:
    input:
        skt=rules.sketch.output.skt,
        fa=FA,
    output:
        bed=pjoin(WD, "k{k}", "reference-anchors.bed"),
    log:
        time = pjoin(WD, "TIMES", "k{k}", "kan-reference.time"),
    shell:
        """
        /usr/bin/time -vo {log.time} ../pansv kan -k{wildcards.k} {input.skt} {input.fa} > {output.bed}
        """

rule kan_py:
    input:
        bed=rules.kan_ref.output.bed,
        fai=FA + ".fai",
    output:
        pjoin(WD, "k{k}", "reference.kan.png"),
    params:
        prefix=pjoin(WD, "k{k}", "reference.kan"),
    conda: "workflow/envs/seaborn.yml"
    shell:
        """
        python3 ../scripts/kan_hist.py {input.bed} {input.fai} -o {params.prefix}
        """


rule kan_reads:
    input:
        skt=rules.sketch.output.skt,
        fq=FQ,
    output:
        nuk=pjoin(WD, "k{k}", "sample.nuk"),
    log:
        time = pjoin(WD, "TIMES", "k{k}", "kan-reads.time"),
    shell:
        """
        /usr/bin/time -vo {log.time} ../pansv chreads -k{wildcards.k} {input.skt} {input.fq} > {output.nuk}
        """

rule kan_reads_py:
    input:
        nuk=rules.kan_reads.output.nuk,
    output:
        png=pjoin(WD, "k{k}", "sample.nuk.png"),
        txt=pjoin(WD, "k{k}", "sample.nuk.txt"),
    conda: "workflow/envs/seaborn.yml"
    shell:
        """
        python3 ../scripts/ran_hist.py {input.nuk} {output.png} > {output.txt}
        """
