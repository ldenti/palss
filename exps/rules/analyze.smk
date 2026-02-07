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


# ====================================================================


rule index_haplotypes:
    input:
        fa=pjoin(WD, sample + "-haps.fa"),
    output:
        pjoin(WD, sample + "-haps.fa.sa"),
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        bwa index {input.fa}
        """


rule gaf2fa:
    input:
        gaf=pjoin(WD, "n{n}", "palss-{t}", "resulting-consensus.d{d}.w{w}.gaf"),
    output:
        fa=pjoin(WD, "n{n}", "palss-{t}", "resulting-consensus.d{d}.w{w}.fa"),
    shell:
        """
        cut -f1,17 {input.gaf} | sed "s/^/>/" | sed "s/\\tqs:Z:/\\n/g" > {output.fa}
        """


rule align_consensuses_to_haplotypes:
    input:
        hfa=pjoin(WD, sample + "-haps.fa"),
        hsa=pjoin(WD, sample + "-haps.fa.sa"),
        qfa=rules.gaf2fa.output.fa,
    output:
        bam=pjoin(WD, "n{n}", "palss-{t}", "resulting-consensus.d{d}.w{w}.bam"),
    conda:
        "../envs/bwa.yaml"
    threads: workflow.cores
    shell:
        """
        bwa mem -t{threads} {input.hfa} {input.qfa} | samtools view -bS | samtools sort > {output.bam}
        samtools index {output.bam}
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
        expand(
            pjoin(WD, "n{n}", "graphaligner", "original-{t}.gaf"),
            t=["full", "oneout"],
            n=Ns,
        ),
        expand(pjoin(WD, "n{n}", "graphaligner", "mgcactus.gaf"), n=Ns),
    output:
        csv=pjoin(WD, "support.csv"),
    shell:
        """
        python3 ./utils/get_support.py {WD} {sample} > {output.csv}
        """


rule get_nm:
    input:
        pjoin(WD, sample + "-reads.bam"),
        pjoin(WD, sample + "-reads.tohaps.bam"),
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
    conda:
        "../envs/pysam.yaml"
    shell:
        """
        python3 ./utils/get_nm.py {WD} > {output.csv}
        """
