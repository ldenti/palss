"""
  * pangenome: pjoin(WD, "n{n}", "pangenome-{t}.gbz")
  * reads: pjoin(WD, sample + "-reads.ec.fa"),
"""


rule palss_sketch:
    input:
        gbz=pjoin(WD, "n{n}", "pangenome-{t}.gbz"),
    output:
        skt=pjoin(WD, "n{n}", "pangenome-{t}.d{d}.skt"),
    log:
        time=pjoin(WD, "times", "n{n}", "palss-{t}", "sketch.d{d}.time"),
    threads: workflow.cores / 2
    shell:
        """
        /usr/bin/time -vo {log.time} ../palss sketch -d{wildcards.d} {input.gbz} > {output.skt}
        """


rule palss_fmd:
    input:
        gbz=pjoin(WD, "n{n}", "pangenome-{t}.gbz"),
    output:
        fmd=pjoin(WD, "n{n}", "pangenome-{t}.fmd"),
    threads: workflow.cores
    log:
        time=pjoin(WD, "times", "n{n}", "palss-{t}", "fmd.time"),
    shell:
        """
        # -m 2G
        /usr/bin/time -vo {log.time} sh -c 'LD_LIBRARY_PATH="$PWD/../lib" ../build/gbwtgraph-prefix/src/gbwtgraph/bin/gbz_extract {input.gbz} | ../build/rb3-prefix/src/rb3/ropebwt3 build -t {threads} -Ld - > {output.fmd}'
        """


rule palss_search:
    input:
        gbz=pjoin(WD, "n{n}", "pangenome-{t}.gbz"),
        skt=rules.palss_sketch.output.skt,
        fmd=rules.palss_fmd.output.fmd,
        fa=pjoin(WD, sample + "-reads.ec.fa"),
    output:
        sfs=pjoin(WD, "n{n}", "palss-{t}", "specific_strings.d{d}.txt"),
    threads: workflow.cores
    log:
        time=pjoin(WD, "times", "n{n}", "palss-{t}", "search.d{d}.time"),
    shell:
        """
        /usr/bin/time -vo {log.time} ../palss sfs -@{threads} {input.gbz} {input.skt} {input.fmd} {input.fa} > {output.sfs}
        """


rule palss_sam:
    input:
        gbz=pjoin(WD, "n{n}", "pangenome-{t}.gbz"),
        sfs=rules.palss_search.output.sfs,
    output:
        bam=pjoin(WD, "n{n}", "palss-{t}", "specific_strings.d{d}.bam"),
    shell:
        """
        ../palss sam {input.gbz} {input.sfs} | samtools view -bS | samtools sort > {output.bam}
        samtools index {output.bam}
        """


rule palss_align:
    input:
        gbz=pjoin(WD, "n{n}", "pangenome-{t}.gbz"),
        sfs=rules.palss_search.output.sfs,
    output:
        gaf=pjoin(WD, "n{n}", "palss-{t}", "consensus.d{d}.gaf"),
    log:
        time=pjoin(WD, "times", "n{n}", "palss-{t}", "align.d{d}.time"),
    threads: workflow.cores / 4
    shell:
        """
        /usr/bin/time -vo {log.time} ../palss align {input.gbz} {input.sfs} > {output.gaf}
        """


rule palss_augment:
    input:
        gfa=pjoin(WD, "n{n}", "pangenome-{t}.gfa"),
        gaf=rules.palss_align.output.gaf,
    output:
        gfa=pjoin(WD, "n{n}", "palss-{t}", "pangenome-augmented.d{d}.w{w}.gfa"),
    params:
        wd=pjoin(WD, "n{n}", "palss-{t}", "augment.d{d}.w{w}.wd"),
    conda:
        "../envs/graphaligner.yaml"
    log:
        time=pjoin(WD, "times", "n{n}", "palss-{t}", "augment.d{d}.w{w}.time"),
    threads: workflow.cores
    shell:
        """
        /usr/bin/time -vo {log.time} bash ../scripts/augment.sh {input.gfa} {input.gaf} {wildcards.w} {params.wd} {threads} > {output.gfa}
        """
