# align PALSS consensuses to original haplotypes
#######################################################################
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


rule align_consensuses_to_haplotypes:
    input:
        hfa=pjoin(WD, sample + "-haps.fa"),
        hsa=pjoin(WD, sample + "-haps.fa.sa"),
        qfa=pjoin(WD, "n{n}", "palss-{t}", "resulting-consensus.d{d}.w{w}.fa"),
    output:
        bam=pjoin(
            WD, "n{n}", "palss-{t}", "resulting-consensus.d{d}.w{w}.to-contigs.bam"
        ),
    conda:
        "../envs/bwa.yaml"
    threads: workflow.cores
    shell:
        """
        bwa mem -t{threads} {input.hfa} {input.qfa} | samtools view -bS | samtools sort > {output.bam}
        samtools index {output.bam}
        """


# # align unanchored contigs to both graphs
# #######################################################################
# rule align_unanchored_contigs:
#     input:
#         gfa=pjoin(WD, "n{n}", "pangenome-{t}.gfa"),
#         fa=pjoin(
#             WD,
#             "n{n}",
#             "palss-oneout",
#             "specific_strings.d{d}.txt.reads_with_unanchored.bp.p_ctg.fa",
#         ),
#     output:
#         gaf=pjoin(
#             WD,
#             "n{n}",
#             "palss-oneout",
#             "specific_strings.d{d}.txt.reads_with_unanchored.bp.p_ctg.to-{t}.gaf",
#         ),
#     conda:
#         "../envs/graphaligner.yaml"
#     threads: workflow.cores / 2
#     shell:
#         """
#         GraphAligner --graph {input.gfa} --reads {input.fa} --alignments-out {output.gaf} --preset vg --threads {threads}
#         """


# Split true contigs and align to all graphs and reference
#######################################################################


rule split_haplotypes:
    input:
        fa=pjoin(WD, sample + "-haps.fa"),
    output:
        fa=pjoin(WD, sample + "-haps.50k-overlapping.fa"),
    shell:
        """
        python3 ./utils/split_haplotypes.py {input.fa} > {output.fa}
        """


rule hapsegs_to_original:
    input:
        gfa=pjoin(WD, "n{n}", "pangenome-{t}.gfa"),
        fa=rules.split_haplotypes.output.fa,
    output:
        gaf=pjoin(WD, "n{n}", "truecontigs-aln", "original-{t}.gaf"),
    conda:
        "../envs/graphaligner.yaml"
    threads: workflow.cores / 2
    shell:
        """
        GraphAligner --graph {input.gfa} --reads {input.fa} --alignments-out {output.gaf} --preset vg --threads {threads}
        """


rule hapsegs_to_palss:
    input:
        gfa=pjoin(WD, "n{n}", "palss-{t}", "pangenome-augmented.d{d}.w{w}.id{iden}.gfa"),
        fa=rules.split_haplotypes.output.fa,
    output:
        gaf=pjoin(WD, "n{n}", "truecontigs-aln", "palss-{t}.d{d}.w{w}.id{iden}.gaf"),
    conda:
        "../envs/graphaligner.yaml"
    threads: workflow.cores / 2
    shell:
        """
        GraphAligner --graph {input.gfa} --reads {input.fa} --alignments-out {output.gaf} --preset vg --threads {threads}
        """


rule hapsegs_to_mgc:
    input:
        gfa=pjoin(WD, "n{n}", "pangenome-mgcactus.gfa"),
        fa=rules.split_haplotypes.output.fa,
    output:
        gaf=pjoin(WD, "n{n}", "truecontigs-aln", "mgcactus.gaf"),
    conda:
        "../envs/graphaligner.yaml"
    threads: workflow.cores / 2
    shell:
        """
        GraphAligner --graph {input.gfa} --reads {input.fa} --alignments-out {output.gaf} --preset vg --threads {threads}
        """


rule hapsegs_to_ref:
    input:
        fa=FA,
        qfa=rules.split_haplotypes.output.fa,
    output:
        bam=pjoin(WD, sample + "-haps.50k-overlapping.bam"),
    conda:
        "../envs/minimap2.yaml"
    threads: workflow.cores / 2
    shell:
        """
        minimap2 -t{threads} --MD -ax asm5 --eqx {input.fa} {input.qfa} | samtools view -bS | samtools sort > {output.bam}
        samtools index {output.bam}
        """
