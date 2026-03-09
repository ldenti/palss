"""
  * pangenome: pjoin(WD, "n{n}", "pangenome-{t}.gbz")
  * reads: pjoin(WD, sample + "-reads.fq"), pjoin(WD, sample + "-reads.ec.fa"),
"""


rule palss_sketch:
    input:
        gbz=pjoin(WD, "n{n}", "pangenome-{t}.gbz"),
    output:
        skt=pjoin(WD, "n{n}", "pangenome-{t}.d{d}.skt"),
    log:
        time=pjoin(WD, "times", "n{n}", "palss-{t}", "sketch.d{d}.time"),
    threads: workflow.cores
    shell:
        """
        /usr/bin/time -vo {log.time} ../palss sketch -@{threads} -d{wildcards.d} {input.gbz} > {output.skt}
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
        /usr/bin/time -vo {log.time} sh -c 'LD_LIBRARY_PATH="$PWD/../lib" ../build/gbwtgraph-prefix/src/gbwtgraph/bin/gbz_extract {input.gbz} | ../build/rb3-prefix/src/rb3/ropebwt3 build -t {threads} -m 2G -Ld - > {output.fmd}'
        """


rule palss_search:
    input:
        gbz=pjoin(WD, "n{n}", "pangenome-{t}.gbz"),
        skt=rules.palss_sketch.output.skt,
        fmd=rules.palss_fmd.output.fmd,
        reads=pjoin(WD, sample + "-reads.ec.fa"),
    output:
        sfs=pjoin(WD, "n{n}", "palss-{t}", "specific_strings.d{d}.txt"),
    threads: workflow.cores
    log:
        time=pjoin(WD, "times", "n{n}", "palss-{t}", "search.d{d}.time"),
    shell:
        """
        /usr/bin/time -vo {log.time} ../palss sfs -@{threads} {input.gbz} {input.skt} {input.fmd} {input.reads} > {output.sfs}
        """


rule palss_sam:
    input:
        gbz=pjoin(WD, "n{n}", "pangenome-{t}.gbz"),
        sfs=pjoin(WD, "n{n}", "palss-{t}", "specific_strings.d{d}.txt"),
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
        sfs=pjoin(WD, "n{n}", "palss-{t}", "specific_strings.d{d}.txt"),
    output:
        gaf=pjoin(WD, "n{n}", "palss-{t}", "consensus.d{d}.gaf"),
    log:
        time=pjoin(WD, "times", "n{n}", "palss-{t}", "align.d{d}.time"),
    threads: workflow.cores
    shell:
        """
        /usr/bin/time -vo {log.time} ../palss align -@{threads} {input.gbz} {input.sfs} > {output.gaf}
        """


rule palss_augment:
    input:
        gfa=pjoin(WD, "n{n}", "pangenome-{t}.gfa"),
        gaf=rules.palss_align.output.gaf,
    output:
        gfa=pjoin(WD, "n{n}", "palss-{t}", "pangenome-augmented.d{d}.w{w}.gfa"),
        gaf=pjoin(WD, "n{n}", "palss-{t}", "resulting-consensus.d{d}.w{w}.gaf"),
    params:
        wd=pjoin(WD, "n{n}", "palss-{t}", "augment.d{d}.w{w}.wd"),
        log=pjoin(WD, "n{n}", "palss-{t}", "augment.d{d}.w{w}.log"),
    conda:
        "../envs/graphaligner.yaml"
    log:
        time=pjoin(WD, "times", "n{n}", "palss-{t}", "augment.d{d}.w{w}.time"),
    threads: max(1, workflow.cores / 4)
    shell:
        """
        /usr/bin/time -vo {log.time} bash ../scripts/augment.sh {input.gfa} {input.gaf} {wildcards.w} {params.wd} {threads} > {output.gfa} 2> {params.log}
        mv {params.wd}/resulting_consensus.gaf {output.gaf}
        """


rule palss_refine:
    input:
        gfa=rules.palss_augment.output.gfa,
        fa=pjoin(WD, sample + "-reads.ec.fa"),
        sfs=rules.palss_search.output.sfs,
    output:
        gfa=pjoin(WD, "n{n}", "palss-{t}", "pangenome-augmented.d{d}.w{w}.id{iden}.gfa"),
    params:
        wd=pjoin(WD, "n{n}", "palss-{t}", "augment-refine.d{d}.w{w}.id{iden}.wd"),
        log=pjoin(WD, "n{n}", "palss-{t}", "augment-refine.d{d}.w{w}.id{iden}.log"),
    conda:
        "../envs/graphaligner.yaml"
    log:
        time=pjoin(
            WD,
            "times",
            "n{n}",
            "palss-{t}",
            "augment-refine.d{d}.w{w}.id{iden}.time",
        ),
    threads: workflow.cores
    shell:
        """
        /usr/bin/time -vo {log.time} bash ../scripts/refine_augmentation.sh {input.gfa} {input.fa} {input.sfs} {params.wd} {threads} {wildcards.iden} > {output.gfa} 2> {params.log}
        """


rule gaf2fa:
    input:
        gaf=rules.palss_augment.output.gaf,
    output:
        fa=pjoin(WD, "n{n}", "palss-{t}", "resulting-consensus.d{d}.w{w}.fa"),
    shell:
        """
        cut -f1,17 {input.gaf} | sed "s/^/>/" | sed "s/\\tqs:Z:/\\n/g" > {output.fa}
        """
