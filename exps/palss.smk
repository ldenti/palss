# from snakemake.utils import min_version
# min_version("6.4.1")

from os.path import join as pjoin

##### config file #####
# configfile: "config/config.yml"

Ks = [31]
Ds = [0.1, 0.5, 1.0]
Ws = [1,2,3]

GBZ = config["gbz"]
GFA = config["gfa"]
FQ = config["fq"]
WD=config["wd"]

wildcard_constraints:
    k=r"\d+",
    w=r"\d+",
    d=r"\d\.\d+"

rule run:
    input:
        expand(pjoin(WD, "pangenome-augmented.k{k}.d{d}.w{w}.gfa"), k=Ks, d=Ds, w=Ws)


rule palss_ec:
    input:
        fq=FQ,
    output:
        fa=pjoin(WD, "reads.ec.fa"),
    params:
        prefix=pjoin(WD, "reads"),
    threads: workflow.cores
    conda:
        "./envs/hifiasm.yml"
    log:
        time=pjoin(WD, "times", "ec.time"),
    shell:
        """
        /usr/bin/time -vo {log.time} hifiasm -t{threads} --write-ec --bin-only {input.fq} -o {params.prefix}
        """


rule palss_sketch:
    input:
        gbz=GBZ,
    output:
        skt=pjoin(WD, "pangenome.k{k}.d{d}.skt"),
    log:
        time=pjoin(WD, "times", "sketch.k{k}.d{d}.time"),
    shell:
        """
        /usr/bin/time -vo {log.time} ../palss sketch -d{wildcards.d} -k{wildcards.k} {input.gbz} > {output.skt}
        """


rule palss_fmd:
    input:
        gbz=GBZ,
    output:
        fmd=pjoin(WD, "pangenome.fmd")
    threads: workflow.cores
    log:
        time=pjoin(WD, "times", "fmd.time"),
    shell:
        """
        /usr/bin/time -vo {log.time} sh -c 'LD_LIBRARY_PATH="$PWD/../lib" ../build/gbwtgraph-prefix/src/gbwtgraph/bin/gbz_extract {input.gbz} | ../build/rb3-prefix/src/rb3/ropebwt3 build -t {threads} -m 2G -Ld - > {output.fmd}'
        """


rule palss_search:
    input:
        gbz=GBZ,
        skt=rules.palss_sketch.output.skt,
        fmd=rules.palss_fmd.output.fmd,
        fa=rules.palss_ec.output.fa,
    output:
        sfs=pjoin(WD, "specific_strings.k{k}.d{d}.txt"),
    threads: workflow.cores
    log:
        time=pjoin(WD, "times", "search.k{k}.d{d}.time"),
    shell:
        """
        /usr/bin/time -vo {log.time} ../palss sfs -@{threads} {input.gbz} {input.skt} {input.fmd} {input.fa} > {output.sfs}
        """

rule palss_align:
    input:
        gbz=GBZ,
        sfs=rules.palss_search.output.sfs,
    output:
        gaf=pjoin(WD, "consensus.k{k}.d{d}.gaf"),
    log:
        time=pjoin(WD, "times", "align.k{k}.d{d}.time"),
    shell:
        """
        /usr/bin/time -vo {log.time} ../palss align {input.gbz} {input.sfs} > {output.gaf}
        """


rule palss_augment:
    input:
        gfa=GFA,
        gaf=rules.palss_align.output.gaf,
    output:
        gfa=pjoin(WD, "pangenome-augmented.k{k}.d{d}.gfa")
    # conda:
    #     "./envs/vg.yml"
    log:
        time=pjoin(WD, "times", "augment.k{k}.d{d}.time"),
    shell:
        """
        /usr/bin/time -vo {log.time} vg augment --include-paths --min-coverage 1 --gaf {input.gfa} {input.gaf} | vg view - > ./example/graph.augmented.gfa
        """

rule palss_clean:
    input:
        gfa=rules.palss_augment.output.gfa,
        gaf=rules.palss_align.output.gaf,
    output:
        gfa=pjoin(WD, "pangenome-augmented.k{k}.d{d}.w{w}.gfa")
    log:
        time=pjoin(WD, "times", "clean.k{k}.d{d}.w{w}.time"),
    shell:
        """
        /usr/bin/time -vo {log.time} python3 ../clean_augmented_gfa.py {input.gfa} {input.gaf} {wildcards.w} > {output.gfa}
        """