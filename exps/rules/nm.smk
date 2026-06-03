### NM (AKA RECALL)
rule get_complex_contigs:
    input:
        bam=pjoin(WD, sample + "-haps.{size}-overlapping.bam"),
        bed=BED,
    output:
        txt=pjoin(WD, sample + "-haps.{size}-overlapping.complex.list"),
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        bedtools intersect -f 0.33 -u -a {input.bam} -b {input.bed} | samtools view | cut -f1 > {output.txt}
        """


rule get_nm:
    input:
        expand(
            pjoin(WD, "n{n}", "tables-nm", "{graph}.{size}.csv"),
            n=Ns,
            graph=["full", "oneout"],
            size=SIZES,
        ),
        expand(
            pjoin(WD, "n{n}", "tables-nm", "mgc.cov{cov}.{size}.csv"),
            n=Ns,
            cov=coverages,
            size=SIZES,
        ),
        #
        expand(
            pjoin(
                WD,
                "n{n}",
                "tables-nm",
                "palss-{graph}.d{d}.w{w}.cov{cov}.{size}.csv",
            ),
            n=Ns,
            cov=coverages,
            graph=["oneout"],
            d=Ds,
            w=Ws,
            size=SIZES,
        ),
        # expand(
        #     pjoin(
        #         WD,
        #         "n{n}",
        #         "tables-nm",
        #         "palss-{graph}.d{d}.w{w}.c{c}.m{m}.{size}.csv",
        #     ),
        #     n=Ns,
        #     graph=["oneout"],
        #     d=Ds,
        #     w=Ws,
        #     c=Cs,
        #     m=Ms,
        #     size=SIZES,
        # ),
        # XXX: this is useless, we could do this just once
        expand(
            pjoin(WD, "n{n}", "tables-nm", "reference.{size}.csv"),
            n=Ns,
            cov=coverages,
            size=SIZES,
        ),
    output:
        pjoin(WD, "nm.csv"),
    shell:
        """
        head -1 {input[0]} > {output}
        for i in {input} ; do sed '1d' $i ; done >> {output}
        """


rule get_nm_reference:
    input:
        bam=pjoin(WD, sample + "-haps.{size}-overlapping.bam"),
        txt=pjoin(WD, sample + "-haps.{size}-overlapping.complex.list"),
    output:
        csv=pjoin(WD, "n{n}", "tables-nm", "reference.{size}.csv"),
    params:
        cov=",".join([str(x) for x in coverages]),
    conda:
        "../envs/pysam.yaml"
    threads: workflow.cores / 4
    shell:
        """
        python3 ./utils/get_nm.py -t reference -l {wildcards.size} -c {params.cov} -n {wildcards.n} {input.bam} {input.txt} > {output.csv}
        """


rule get_nm_original:
    input:
        gaf=pjoin(WD, "n{n}", "truecontigs-aln", "original-{graph}.{size}.gaf"),
        txt=pjoin(WD, sample + "-haps.{size}-overlapping.complex.list"),
    output:
        csv=pjoin(WD, "n{n}", "tables-nm", "{graph}.{size}.csv"),
    params:
        cov=",".join([str(x) for x in coverages]),
    conda:
        "../envs/pysam.yaml"
    threads: workflow.cores / 4
    shell:
        """
        python3 ./utils/get_nm.py -t original-{wildcards.graph} -l {wildcards.size} -c {params.cov} -n {wildcards.n} {input.gaf} {input.txt}  > {output.csv}
        """


rule get_nm_mgc:
    input:
        gaf=pjoin(WD, "n{n}", "truecontigs-aln", "mgcactus.cov{cov}.{size}.gaf"),
        txt=pjoin(WD, sample + "-haps.{size}-overlapping.complex.list"),
    output:
        csv=pjoin(WD, "n{n}", "tables-nm", "mgc.cov{cov}.{size}.csv"),
    conda:
        "../envs/pysam.yaml"
    threads: workflow.cores / 4
    shell:
        """
        python3 ./utils/get_nm.py -t mgcactus -l {wildcards.size} -c {wildcards.cov} -n {wildcards.n} {input.gaf} {input.txt} > {output.csv}
        """


rule get_nm_palss:
    input:
        gaf=pjoin(
            WD,
            "n{n}",
            "truecontigs-aln",
            "palss-{graph}.d{d}.w{w}.cov{cov}.{size}.gaf",
        ),
        txt=pjoin(WD, sample + "-haps.{size}-overlapping.complex.list"),
    output:
        csv=pjoin(
            WD,
            "n{n}",
            "tables-nm",
            "palss-{graph}.d{d}.w{w}.cov{cov}.{size}.csv",
        ),
    conda:
        "../envs/pysam.yaml"
    threads: workflow.cores / 4
    shell:
        """
        python3 ./utils/get_nm.py -t palss-d{wildcards.d} -l {wildcards.size} -c {wildcards.cov} -n {wildcards.n} {input.gaf} {input.txt} > {output.csv}
        """


# rule get_nm_palss_refine:
#     input:
#         gaf=pjoin(
#             WD,
#             "n{n}",
#             "truecontigs-aln",
#             "palss-{graph}.d{d}.w{w}.c{c}.m{m}.{size}.gaf",
#         ),
#         txt=pjoin(WD, sample + "-haps.{size}-overlapping.complex.list"),
#     output:
#         csv=pjoin(
#             WD,
#             "n{n}",
#             "tables-nm",
#             "palss-{graph}.d{d}.w{w}.c{c}.m{m}.{size}.csv",
#         ),
#     conda:
#         "../envs/pysam.yaml"
#     threads: workflow.cores / 4
#     shell:
#         """
#         python3 ./utils/get_nm.py {input.gaf} {input.txt} {wildcards.n} > {output.csv}
#         """
