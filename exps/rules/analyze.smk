"""
  * pangenome: pjoin(WD, "n{n}", "palss-{t}", "pangenome-augmented.d{d}.w{w}.gfa")
  * reads: pjoin(WD, sample + "-reads.ec.fa"),
"""


rule graphaligner_original:
    input:
        gfa=pjoin(WD, "n{n}", "pangenome-{t}.gfa"),
        reads=pjoin(WD, sample + "-reads.fq.gz"),
    output:
        gaf=pjoin(WD, "n{n}", "graphaligner", "original-{t}.gaf"),
    conda:
        "../envs/graphaligner.yaml"
    threads: workflow.cores / 2
    shell:
        """
        GraphAligner --graph {input.gfa} --reads {input.reads} --alignments-out {output.gaf} --preset vg --threads {threads}
        """


rule graphaligner_postpalss:
    input:
        gfa=pjoin(WD, "n{n}", "palss-{t}", "pangenome-augmented.d{d}.w{w}.gfa"),
        reads=pjoin(WD, sample + "-reads.fq.gz"),
    output:
        gaf=pjoin(WD, "n{n}", "graphaligner", "palss-{t}.d{d}.w{w}.gaf"),
    conda:
        "../envs/graphaligner.yaml"
    threads: workflow.cores / 2
    shell:
        """
        GraphAligner --graph {input.gfa} --reads {input.reads} --alignments-out {output.gaf} --preset vg --threads {threads}
        """


rule graphaligner_postmgc:
    input:
        gfa=pjoin(WD, "n{n}", "pangenome-mgcactus.gfa"),
        reads=pjoin(WD, sample + "-reads.fq.gz"),
    output:
        gaf=pjoin(WD, "n{n}", "graphaligner", "mgcactus.gaf"),
    conda:
        "../envs/graphaligner.yaml"
    threads: workflow.cores / 2
    shell:
        """
        GraphAligner --graph {input.gfa} --reads {input.reads} --alignments-out {output.gaf} --preset vg --threads {threads}
        """


rule get_support:
    input:
        expand(
            pjoin(WD, "n{n}", "graphaligner", "palss-{t}.d{d}.w{w}.gaf"),
            t=["full", "oneout"],
            n=Ns,
            d=Ds,
            w=Ws,
        ),
    output:
        csv=pjoin(WD, "support.csv"),
    shell:
        """
        python3 ./utils/get_support.py {WD} > {output.csv}
        """


rule get_nm:
    input:
        expand(
            pjoin(WD, "n{n}", "graphaligner", "original-{t}.gaf"),
            t=["full", "oneout"],
            n=Ns,
        ),
        expand(
            pjoin(WD, "n{n}", "graphaligner", "palss-{t}.d{d}.w{w}.gaf"),
            t=["full", "oneout"],
            n=Ns,
            d=Ds,
            w=Ws,
        ),
        expand(pjoin(WD, "n{n}", "graphaligner", "mgcactus.gaf"), n=Ns),
    output:
        csv=pjoin(WD, "nm.csv"),
    shell:
        """
        python3 ./utils/get_nm.py {WD} > {output.csv}
        """
