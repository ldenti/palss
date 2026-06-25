"""
  * pangenome: pjoin(WD, "n{n}", "pangenome-{graph}.gbz")
  * reads: pjoin(WD, sample + "-cov{cov}.fq.gz"),
  * corrected reads: pjoin(WD, sample + "-cov{cov}.ec.fa"),
"""


rule palss_sketch:
    input:
        gbz=pjoin(WD, "n{n}", "pangenome-{graph}.gbz"),
    output:
        skt=pjoin(WD, "palss", "pangenome-{graph}.n{n}.d{d}.skt"),
    log:
        time=pjoin(WD, "times", "palss", "sketch-{graph}.n{n}.d{d}.time"),
    threads: workflow.cores
    shell:
        """
        /usr/bin/time -vo {log.time} ../palss sketch -@{threads} -d{wildcards.d} {input.gbz} > {output.skt}
        """


rule palss_fmd:
    input:
        gbz=pjoin(WD, "n{n}", "pangenome-{graph}.gbz"),
    output:
        fmd=pjoin(WD, "palss", "pangenome-{graph}.n{n}.fmd"),
    threads: workflow.cores
    log:
        time=pjoin(WD, "times", "palss", "fmd-{graph}.n{n}.time"),
    shell:
        """
        /usr/bin/time -vo {log.time} sh -c 'LD_LIBRARY_PATH="$PWD/../lib" ../build/gbwtgraph-prefix/src/gbwtgraph/bin/gbz_extract {input.gbz} | ../build/rb3-prefix/src/rb3/ropebwt3 build -t {threads} -Ld - > {output.fmd}'
        """


rule palss_search:
    input:
        gbz=pjoin(WD, "n{n}", "pangenome-{graph}.gbz"),
        skt=rules.palss_sketch.output.skt,
        fmd=rules.palss_fmd.output.fmd,
        reads=pjoin(WD, sample + "-cov{cov}.ec.fa"),
    output:
        sfs=pjoin(WD, "palss", "n{n}", "cov{cov}", "specificstrings-{graph}.d{d}.txt"),
    threads: workflow.cores
    log:
        time=pjoin(WD, "times", "palss", "search-{graph}.n{n}.cov{cov}.d{d}.time"),
    shell:
        """
        /usr/bin/time -vo {log.time} ../palss sfs -@{threads} {input.gbz} {input.skt} {input.fmd} {input.reads} > {output.sfs}
        """


# rule palss_sam:
#     input:
#         gbz=pjoin(WD, "n{n}", "pangenome-{graph}.gbz"),
#         sfs=pjoin(WD, "n{n}", "palss-{graph}", "cov{cov}", "specific_strings.d{d}.txt"),
#     output:
#         bam=pjoin(WD, "n{n}", "palss-{graph}", "cov{cov}", "specific_strings.d{d}.bam"),
#     shell:
#         """
#         ../palss sam {input.gbz} {input.sfs} | samtools view -bS | samtools sort > {output.bam}
#         samtools index {output.bam}
#         """


rule palss_align:
    input:
        gbz=pjoin(WD, "n{n}", "pangenome-{graph}.gbz"),
        sfs=rules.palss_search.output.sfs,
    output:
        gaf=pjoin(WD, "palss", "n{n}", "cov{cov}", "consensus-{graph}.d{d}.gaf"),
    log:
        time=pjoin(WD, "times", "palss", "align-{graph}.n{n}.cov{cov}.d{d}.time"),
    threads: workflow.cores
    shell:
        """
        /usr/bin/time -vo {log.time} ../palss align -@{threads} {input.gbz} {input.sfs} > {output.gaf}
        """


rule palss_augment:
    input:
        pg=pjoin(WD, "n{n}", "pangenome-{graph}.pg"),
        gaf=rules.palss_align.output.gaf,
    output:
        gfa=pjoin(WD, "palss", "n{n}", "cov{cov}", "augmented-{graph}.d{d}.w{w}.gfa"),
        gaf=pjoin(
            WD, "palss", "n{n}", "cov{cov}", "anchoredconsensus-{graph}.d{d}.w{w}.gaf"
        ),
    params:
        wd=pjoin(WD, "palss", "n{n}", "cov{cov}", "augment-wd.{graph}.d{d}.w{w}"),
    log:
        log=pjoin(WD, "palss", "n{n}", "cov{cov}", "augment-wd.{graph}.d{d}.w{w}.log"),
        time=pjoin(WD, "times", "palss", "augment-{graph}.n{n}.cov{cov}.d{d}.w{w}.time"),
    threads: workflow.cores
    shell:
        """
        /usr/bin/time -vo {log.time} ../palss augment -s {wildcards.w} -w {params.wd} -g {output.gaf} {input.pg} {input.gaf} > {output.gfa}.unchop 2> {log.log}
        vg mod --unchop {output.gfa}.unchop > {output.gfa}
        """


rule gaf2fa:
    input:
        gaf=rules.palss_augment.output.gaf,
    output:
        fa=pjoin(
            WD, "palss", "n{n}", "cov{cov}", "anchoredconsensus-{graph}.d{d}.w{w}.fa"
        ),
    shell:
        """
        cut -f1,17 {input.gaf} | sed "s/^/>/" | sed "s/\\tqs:Z:/\\n/g" > {output.fa}
        """
