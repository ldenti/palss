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
    threads: workflow.cores
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
    threads: workflow.cores
    shell:
        """
        GraphAligner --graph {input.gfa} --reads {input.reads} --alignments-out {output.gaf} --preset vg --threads {threads}
        """


rule graphaligner_postmgc:
    input:
        gfa=pjoin(WD, "n{n}", "pangenome-mgcactus.gfa"),
        reads=pjoin(WD, sample + "-reads.fq.gz"),
    output:
        gaf=pjoin(WD, "n{n}", "graphaligner", "mgc.gaf"),
    conda:
        "../envs/graphaligner.yaml"
    threads: workflow.cores
    shell:
        """
        GraphAligner --graph {input.gfa} --reads {input.reads} --alignments-out {output.gaf} --preset vg --threads {threads}
        """
