### SUPPORT (AKA PRECISION)
rule get_support:
    input:
        expand(
            pjoin(WD, "n{n}", "tables-supp", "{graph}.{size}.csv"),
            n=Ns,
            graph=["full", "oneout"],
            size=SIZES,
        ),
        expand(
            pjoin(WD, "n{n}", "tables-supp", "mgc.cov{cov}.{size}.csv"),
            n=Ns,
            cov=coverages,
            size=SIZES,
        ),
        #
        expand(
            pjoin(
                WD,
                "n{n}",
                "tables-supp",
                "palss-{graph}.d{d}.w{w}.cov{cov}.{size}.csv",
            ),
            n=Ns,
            cov=coverages,
            graph=["oneout"],
            d=Ds,
            w=Ws,
            size=SIZES,
        ),
    output:
        pjoin(WD, "support.csv"),
    shell:
        """
        head -1 {input[0]} > {output}
        for i in {input} ; do sed '1d' $i ; done >> {output}
        """


rule get_support_original:
    input:
        gfa=pjoin(WD, "n{n}", "pangenome-{graph}.gfa"),
        gaf=pjoin(WD, "n{n}", "truecontigs-aln", "original-{graph}.{size}.gaf"),
        bed=BED,
    output:
        csv=pjoin(WD, "n{n}", "tables-supp", "{graph}.{size}.csv"),
    params:
        cov=",".join([str(x) for x in coverages]),
    conda:
        "../envs/intervaltree.yaml"
    threads: workflow.cores / 4
    shell:
        """
        python3 ./utils/get_support.py -t original-{wildcards.graph} -l {wildcards.size} -c {params.cov} -s {sample} -n {wildcards.n} {input.gfa} {input.gaf} > {output.csv}.unflagged
        python3 ./utils/flag_vertices.py {input.gfa} {output.csv}.unflagged {input.bed} > {output.csv}
        """


rule get_support_mgc:
    input:
        gfa=pjoin(WD, "mgcactus", "n{n}", "cov{cov}", "pangenome-mgcactus.gfa"),
        gaf=pjoin(WD, "n{n}", "truecontigs-aln", "mgcactus.cov{cov}.{size}.gaf"),
        bed=BED,
    output:
        csv=pjoin(WD, "n{n}", "tables-supp", "mgc.cov{cov}.{size}.csv"),
    conda:
        "../envs/intervaltree.yaml"
    threads: workflow.cores / 4
    shell:
        """
        python3 ./utils/get_support.py -t mgcactus -l {wildcards.size} -c {wildcards.cov} -s {sample} -n {wildcards.n} {input.gfa} {input.gaf} > {output.csv}.unflagged
        python3 ./utils/flag_vertices.py {input.gfa} {output.csv}.unflagged {input.bed} > {output.csv}
        """


rule get_support_palss:
    input:
        gfa=pjoin(WD, "palss", "n{n}", "cov{cov}", "augmented-{graph}.d{d}.w{w}.gfa"),
        gaf=pjoin(
            WD,
            "n{n}",
            "truecontigs-aln",
            "palss-{graph}.d{d}.w{w}.cov{cov}.{size}.gaf",
        ),
        bed=BED,
    output:
        csv=pjoin(
            WD,
            "n{n}",
            "tables-supp",
            "palss-{graph}.d{d}.w{w}.cov{cov}.{size}.csv",
        ),
    conda:
        "../envs/intervaltree.yaml"
    threads: workflow.cores / 4
    shell:
        """
        python3 ./utils/get_support.py -t palss-d{wildcards.d}-w{wildcards.w} -l {wildcards.size} -c {wildcards.cov} -n {wildcards.n} {input.gfa} {input.gaf} > {output.csv}.unflagged
        python3 ./utils/flag_vertices.py {input.gfa} {output.csv}.unflagged {input.bed} > {output.csv}
        """
