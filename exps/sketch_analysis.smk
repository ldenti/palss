from os.path import join as pjoin


FA = config["fa"]
FQ = config["fq"]
GBZ = config["gbz"]
WD = config["wd"]

Ds = [0.1, 0.25, 0.33, 0.5]


rule run:
    input:
        expand(pjoin(WD, "d{d}", "missed_regions.bed"), d=Ds),
        expand(pjoin(WD, "d{d}", "reads.txt"), d=Ds),


rule palss_sketch:
    input:
        gbz=GBZ,
    output:
        skt=pjoin(WD, "d{d}", "pangenome.skt"),
    log:
        time=pjoin(WD, "times", "palss-sketch.d{d}.time"),
    threads: workflow.cores
    shell:
        """
        /usr/bin/time -vo {log.time} ../palss sketch -@{threads} -d{wildcards.d} {input.gbz} > {output.skt}
        """


rule palss_kan_reference:
    input:
        skt=rules.palss_sketch.output.skt,
        fa=FA,
    output:
        bed=pjoin(WD, "d{d}", "missed_regions.bed"),
    threads: workflow.cores / 2
    shell:
        """
        ../palss kan {input.skt} {input.fa} > {output.bed}
        """


rule palss_kan_reads:
    input:
        skt=rules.palss_sketch.output.skt,
        fq=FQ,
    output:
        txt=pjoin(WD, "d{d}", "reads.txt"),
    threads: workflow.cores / 2
    shell:
        """
        ../palss kan -q {input.skt} {input.fq} > {output.txt}
        """
