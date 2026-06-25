"""
  * pangenome: pjoin(WD, "n{n}", "palss-{graph}", "pangenome-augmented.d{d}.w{w}.gfa")
  * reads: pjoin(WD, sample + "-cov{cov}.fq.gz"),
  * corrected reads: pjoin(WD, sample + "-cov{cov}.ec.fa"),
"""


### DOWNSTREAM (AKA READS ALIGNMENT)
rule get_complex_reads:
    input:
        bam=pjoin(WD, sample + "-cov{cov}.bam"),
        bed=BED,
    output:
        txt=pjoin(WD, sample + "-cov{cov}.overlapping-complex.list"),
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        bedtools intersect -f 0.33 -u -a {input.bam} -b {input.bed} | samtools view | cut -f1 > {output.txt}
        """


rule get_downstream:
    input:
        expand(
            pjoin(WD, "n{n}", "downstream", "palss-{graph}.d{d}.w{w}.cov{cov}.csv"),
            n=Ns,
            cov=coverages,
            graph=["oneout"],
            d=Ds,
            w=Ws,
        ),
        expand(
            pjoin(WD, "n{n}", "downstream", "mgcactus.cov{cov}.csv"),
            n=Ns,
            cov=coverages,
        ),
    output:
        pjoin(WD, "downstream.csv"),
    shell:
        """
        head -1 {input[0]} > {output}
        for i in {input} ; do sed '1d' $i ; done >> {output}
        """


rule graphaligner_reads_palss:
    input:
        gfa=pjoin(WD, "palss", "n{n}", "cov{cov}", "augmented-{graph}.d{d}.w{w}.gfa"),
        reads=pjoin(WD, sample + "-cov{cov}.fq.gz"),
    output:
        gaf=pjoin(WD, "n{n}", "downstream", "palss-{graph}.d{d}.w{w}.cov{cov}.gaf"),
    conda:
        "../envs/graphaligner.yaml"
    threads: workflow.cores / 2
    shell:
        """
        GraphAligner --graph {input.gfa} --reads {input.reads} --alignments-out {output.gaf} --preset vg --threads {threads}
        """


rule graphaligner_reads_mgc:
    input:
        gfa=pjoin(WD, "mgcactus", "n{n}", "cov{cov}", "pangenome-mgcactus.gfa"),
        reads=pjoin(WD, sample + "-cov{cov}.fq.gz"),
    output:
        gaf=pjoin(WD, "n{n}", "downstream", "mgcactus.cov{cov}.gaf"),
    conda:
        "../envs/graphaligner.yaml"
    threads: workflow.cores / 2
    shell:
        """
        GraphAligner --graph {input.gfa} --reads {input.reads} --alignments-out {output.gaf} --preset vg --threads {threads}
        """


rule get_downstream_palss:
    input:
        gaf=rules.graphaligner_reads_palss.output.gaf,
        txt=rules.get_complex_reads.output.txt,
    output:
        csv=pjoin(WD, "n{n}", "downstream", "palss-{graph}.d{d}.w{w}.cov{cov}.csv"),
    conda:
        "../envs/pysam.yaml"
    shell:
        """
        python3 ./utils/get_nm.py -t palss-d{wildcards.d}-w{wildcards.w} -l 0 -c {wildcards.cov} -n {wildcards.n} {input.gaf} {input.txt} > {output.csv}
        """


rule get_downstream_mgc:
    input:
        gaf=rules.graphaligner_reads_mgc.output.gaf,
        txt=rules.get_complex_reads.output.txt,
    output:
        csv=pjoin(WD, "n{n}", "downstream", "mgcactus.cov{cov}.csv"),
    conda:
        "../envs/pysam.yaml"
    shell:
        """
        python3 ./utils/get_nm.py -t mgcactus -l 0 -c {wildcards.cov} -n {wildcards.n} {input.gaf} {input.txt} > {output.csv}
        """
