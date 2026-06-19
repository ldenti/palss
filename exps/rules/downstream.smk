"""
  * pangenome: pjoin(WD, "n{n}", "palss-{graph}", "pangenome-augmented.d{d}.w{w}.gfa")
  * reads: pjoin(WD, sample + "-cov{cov}.fq.gz"),
  * corrected reads: pjoin(WD, sample + "-cov{cov}.ec.fa"),
"""

# rule align_reads:
#     input:
#         fa=FA,
#         fq=pjoin(WD, sample + "-reads.fq.gz"),
#     output:
#         bam=pjoin(WD, sample + "-reads.bam"),
#     conda:
#         "../envs/minimap2.yaml"
#     threads: workflow.cores
#     shell:
#         """
#         minimap2 -t{threads} --MD -ax map-hifi --eqx {input.fa} {input.fq} | samtools view -bS | samtools sort > {output.bam}
#         samtools index {output.bam}
#         """

# rule align_corrected_reads:
#     input:
#         fa=FA,
#         reads=pjoin(WD, sample + "-reads.ec.fa"),
#     output:
#         bam=pjoin(WD, sample + "-reads.ec.bam"),
#     conda:
#         "../envs/minimap2.yaml"
#     threads: workflow.cores
#     shell:
#         """
#         minimap2 -t{threads} --MD -ax map-hifi --eqx {input.fa} {input.reads} | samtools view -bS | samtools sort > {output.bam}
#         samtools index {output.bam}
#         """

# rule align_reads_to_haplotypes:
#     input:
#         fa=pjoin(WD, sample + "-haps.fa"),
#         fq=pjoin(WD, sample + "-reads.fq.gz"),
#     output:
#         bam=pjoin(WD, sample + "-reads.tohaps.bam"),
#     conda:
#         "../envs/minimap2.yaml"
#     threads: workflow.cores
#     shell:
#         """
#         minimap2 -t{threads} --MD -ax map-hifi --eqx {input.fa} {input.fq} | samtools view -bS | samtools sort > {output.bam}
#         samtools index {output.bam}
#         """


# rule graphaligner_original:
#     input:
#         gfa=pjoin(WD, "n{n}", "pangenome-{t}.gfa"),
#         reads=pjoin(WD, sample + "-reads.fq.gz"),
#     output:
#         gaf=pjoin(WD, "n{n}", "graphaligner", "original-{t}.gaf"),
#     conda:
#         "../envs/graphaligner.yaml"
#     threads: workflow.cores / 2
#     shell:
#         """
#         GraphAligner --graph {input.gfa} --reads {input.reads} --alignments-out {output.gaf} --preset vg --threads {threads}
#         """


# rule graphaligner_original_ec:
#     input:
#         gfa=pjoin(WD, "n{n}", "pangenome-{t}.gfa"),
#         reads=pjoin(WD, sample + "-reads.ec.fa"),
#     output:
#         gaf=pjoin(WD, "n{n}", "graphaligner.ec", "original-{t}.gaf"),
#     conda:
#         "../envs/graphaligner.yaml"
#     threads: workflow.cores / 2
#     shell:
#         """
#         GraphAligner --graph {input.gfa} --reads {input.reads} --alignments-out {output.gaf} --preset vg --threads {threads}
#         """


rule graphaligner_postpalss:
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


rule graphaligner_postmgc:
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
