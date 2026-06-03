from os.path import join as pjoin


FA = config["fa"]
FQ = config["fq"]
#
GFA = config["gfa"]
GBZ = config["gbz"]
#
WD = config["wd"]


rule run:
    input:
        pjoin(WD, "minimap2.bam"),
        pjoin(WD, "graphaligner-full.gaf"),
        pjoin(WD, "graphaligner-palss.gaf"),


### PALSS ###########################
rule hifiasm_ec:
    input:
        fq=FQ,
    output:
        fa=pjoin(WD, "reads.ec.fa"),
    params:
        prefix=pjoin(WD, "reads"),
    threads: workflow.cores
    conda:
        "./envs/hifiasm.yaml"
    log:
        time=pjoin(WD, "times", "hifiasm-ec.time"),
        log=pjoin(WD, "hifiasm-ec.log"),
    shell:
        """
        /usr/bin/time -vo {log.time} hifiasm -t{threads} --write-ec --bin-only {input.fq} -o {params.prefix} 2> {log.log}
        """


rule palss_sketch:
    input:
        gbz=GBZ,
    output:
        skt=pjoin(WD, "pangenome.skt"),
    log:
        time=pjoin(WD, "times", "palss-sketch.time"),
    threads: workflow.cores
    shell:
        """
        /usr/bin/time -vo {log.time} ../palss sketch -@{threads} -d0.1 {input.gbz} > {output.skt}
        """


rule palss_fmd:
    input:
        gbz=GBZ,
    output:
        fmd=pjoin(WD, "pangenome.fmd"),
    log:
        time=pjoin(WD, "times", "palss-fmd.time"),
    threads: workflow.cores
    shell:
        """
        /usr/bin/time -vo {log.time} sh -c 'LD_LIBRARY_PATH="$PWD/../lib" ../build/gbwtgraph-prefix/src/gbwtgraph/bin/gbz_extract {input.gbz} | ../build/rb3-prefix/src/rb3/ropebwt3 build -t {threads} -Ld - > {output.fmd}'
        """


rule palss_search:
    input:
        gbz=GBZ,
        skt=rules.palss_sketch.output.skt,
        fmd=rules.palss_fmd.output.fmd,
        reads=rules.hifiasm_ec.output.fa,
    output:
        sfs=pjoin(WD, "specific_strings.txt"),
    threads: workflow.cores
    log:
        time=pjoin(WD, "times", "palss-search.time"),
    shell:
        """
        /usr/bin/time -vo {log.time} ../palss sfs -@{threads} {input.gbz} {input.skt} {input.fmd} {input.reads} > {output.sfs}
        """


rule palss_align:
    input:
        gbz=GBZ,
        sfs=rules.palss_search.output.sfs,
    output:
        gaf=pjoin(WD, "consensus.gaf"),
    log:
        time=pjoin(WD, "times", "palss-align.time"),
    threads: workflow.cores
    shell:
        """
        /usr/bin/time -vo {log.time} ../palss align -@{threads} {input.gbz} {input.sfs} > {output.gaf}
        """


rule palss_augment:
    input:
        gfa=GFA,
        gaf=rules.palss_align.output.gaf,
    output:
        gfa=pjoin(WD, "pangenome-augmented.gfa"),
        gaf=pjoin(WD, "anchored-consensus.gaf"),
    params:
        wd=pjoin(WD, "augment.wd"),
        log=pjoin(WD, "augment.log"),
    conda:
        "./envs/graphaligner.yaml"
    log:
        time=pjoin(WD, "times", "palss-augment.time"),
    threads: workflow.cores
    shell:
        """
        /usr/bin/time -vo {log.time} bash ../scripts/augment.sh {input.gfa} {input.gaf} 2 {params.wd} {threads} > {output.gfa} 2> {params.log}
        mv {params.wd}/resulting_consensus.gaf {output.gaf}
        """


# TODO: refine


### MINIMAP2 ###########################
rule minimap2:
    input:
        fa=FA,
        fq=FQ,
    output:
        bam=pjoin(WD, "minimap2.bam"),
    conda:
        "./envs/minimap2.yaml"
    threads: workflow.cores
    shell:
        """
        minimap2 -ax map-hifi --MD --eqx -t{threads} {input.fa} {input.fq} | samtools view -bS | samtools sort > {output.bam}
        samtools index {output.bam}
        """


### GRAPHALIGNER AGAINST FULL ###########################
rule ga_full:
    input:
        gfa=GFA,
        fq=FQ,
    output:
        gaf=pjoin(WD, "graphaligner-full.gaf"),
    conda:
        "./envs/graphaligner.yaml"
    threads: workflow.cores
    shell:
        """
        GraphAligner --graph {input.gfa} --reads {input.fq} --alignments-out {output.gaf} --preset vg --threads {threads}
        """


# vg giraffe --parameter-preset hifi --haplotype-sampling --gbz-name {input.gbz} --fastq-in {input.fq} > {output.gam}


### GRAPHALIGNER AGAINST PALSS ###########################
rule ga_palss:
    input:
        gfa=rules.palss_augment.output.gfa,
        fq=FQ,
    output:
        gaf=pjoin(WD, "graphaligner-palss.gaf"),
    threads: workflow.cores
    shell:
        """
        GraphAligner --graph {input.gfa} --reads {input.fq} --alignments-out {output.gaf} --preset vg --threads {threads}
        """
