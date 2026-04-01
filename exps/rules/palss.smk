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
        /usr/bin/time -vo {log.time} sh -c 'LD_LIBRARY_PATH="$PWD/../lib" ../build/gbwtgraph-prefix/src/gbwtgraph/bin/gbz_extract {input.gbz} | ../build/rb3-prefix/src/rb3/ropebwt3 build -t {threads} -Ld - > {output.fmd}'
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


# ====================================================================
# ====================================================================
HIFIASM_BIN = "/home/ldenti/code/hifiasm/hifiasm"
# ====================================================================
# ====================================================================


rule palss_refine_1:
    input:
        gfa=rules.palss_augment.output.gfa,
        fa=pjoin(WD, sample + "-reads.ec.fa"),
        sfs=rules.palss_search.output.sfs,
    output:
        fa=pjoin(
            WD,
            "n{n}",
            "palss-{t}",
            "unanchored_contigs.d{d}.w{w}.fa",
        ),
        gaf=pjoin(
            WD,
            "n{n}",
            "palss-{t}",
            "unanchored_contigs.d{d}.w{w}.gaf",
        ),
    params:
        wd=pjoin(WD, "n{n}", "palss-{t}", "augment-refine.d{d}.w{w}.wd"),
        log=pjoin(WD, "n{n}", "palss-{t}", "augment-refine.d{d}.w{w}.log"),
    conda:
        "../envs/refine.yaml"
    threads: workflow.cores
    log:
        time=pjoin(WD, "times", "n{n}", "palss-{t}", "refine_1.d{d}.w{w}.time"),
    shell:
        """
        mkdir -p {params.wd}
        grep -P "^1|^2|^3" {input.sfs} | cut -f2 | sort -u > {params.wd}/reads_with_unanchored.list
        python3 ../scripts/subfa.py {input.fa} {params.wd}/reads_with_unanchored.list > {params.wd}/reads_with_unanchored.fa
        /usr/bin/time -vo {log.time} {HIFIASM_BIN} -r0 -o {params.wd}/reads_with_unanchored -t{threads} {params.wd}/reads_with_unanchored.fa
        awk '/^S/{{print ">"$2;print $3}}' {params.wd}/reads_with_unanchored.bp.p_utg.gfa > {output.fa}
        # awk '/^S/{{print ">"$2;print $3}}' {params.wd}/reads_with_unanchored.bp.p_ctg.gfa > {output.fa}
        # awk '/^S/{{print ">"$2;print $3}}' {params.wd}/reads_with_unanchored.bp.hap1.p_ctg.gfa > {output.fa} # {params.wd}/reads_with_unanchored.bp.haps.p_ctg.fa
        # awk '/^S/{{print ">"$2;print $3}}' {params.wd}/reads_with_unanchored.bp.hap2.p_ctg.gfa >> {output.fa} # {params.wd}/reads_with_unanchored.bp.haps.p_ctg.fa

        #/usr/bin/time -vao {log.time} minimap2 -ax asm10 --MD --eqx -Y -t{threads} {FA} {params.wd}/reads_with_unanchored.bp.haps.p_ctg.fa | samtools view -bS | samtools sort > {params.wd}/reads_with_unanchored.bp.haps.p_ctg.bam
        # samtools view  -F 2048 {params.wd}/reads_with_unanchored.bp.haps.p_ctg.bam | cut -f1 | sort | uniq -c | grep " 1 " | egrep -o "h.+" | grep -A1 -f - --no-group-separator {params.wd}/reads_with_unanchored.bp.haps.p_ctg.fa > {output.fa}
        # samtools view -q50 {params.wd}/reads_with_unanchored.bp.haps.p_ctg.bam | cut -f1 | sort -u | grep -A1 -f - --no-group-separator {params.wd}/reads_with_unanchored.bp.haps.p_ctg.fa > {output.fa}
        /usr/bin/time  -vao {log.time} minigraph -x asm -c {input.gfa} {output.fa} -t{threads} > {output.gaf}
        """


rule palss_refine_2:
    input:
        gbz=pjoin(WD, "n{n}", "pangenome-{t}.gbz"),
        skt=rules.palss_sketch.output.skt,
        fmd=rules.palss_fmd.output.fmd,
        fa=rules.palss_refine_1.output.fa,
    output:
        sfs=pjoin(
            WD,
            "n{n}",
            "palss-{t}",
            "unanchored_contigs.d{d}.w{w}.fa.sfs",
        ),
    threads: workflow.cores
    log:
        time=pjoin(WD, "times", "n{n}", "palss-{t}", "refine_2.d{d}.w{w}.time"),
    shell:
        """
        /usr/bin/time -vo {log.time} ../palss sfs -@{threads} -n {input.gbz} {input.skt} {input.fmd} {input.fa} > {output.sfs}
        """


rule palss_refine_3:
    input:
        gfa=rules.palss_augment.output.gfa,
        fa=rules.palss_refine_1.output.fa,
        gaf=rules.palss_refine_1.output.gaf,
        sfs=rules.palss_refine_2.output.sfs,
    output:
        gaf=pjoin(
            WD,
            "n{n}",
            "palss-{t}",
            "unanchored_contigs.d{d}.w{w}.c{c}.m{m}.gaf",
        ),
    threads: max(1, workflow.cores / 4)
    conda:
        "../envs/biopython.yaml"
    log:
        time=pjoin(
            WD, "times", "n{n}", "palss-{t}", "refine_3.d{d}.w{w}.c{c}.m{m}.time"
        ),
        log=pjoin(
            WD,
            "n{n}",
            "palss-{t}",
            "unanchored_contigs.d{d}.w{w}.c{c}.m{m}.log",
        ),
    shell:
        """
        /usr/bin/time -vo {log.time}  python3 ../scripts/refine.py -c {wildcards.c} -m {wildcards.m} {input.gfa} {input.fa} {input.gaf} {input.sfs} > {output.gaf} 2> {log.log}
        """


rule palss_refine_4:
    input:
        gfa=rules.palss_augment.output.gfa,
        gaf=rules.palss_refine_3.output.gaf,
    output:
        gfa=pjoin(
            WD,
            "n{n}",
            "palss-{t}",
            "pangenome-augmented.d{d}.w{w}.c{c}.m{m}.gfa",
        ),
    threads: max(1, workflow.cores / 4)
    log:
        time=pjoin(
            WD, "times", "n{n}", "palss-{t}", "refine_4.d{d}.w{w}.c{c}.m{m}.time"
        ),
    shell:
        """
        /usr/bin/time -vo {log.time} sh -c "vg augment --min-coverage 1 --gaf {input.gfa} {input.gaf} | vg mod --unchop -" > {output.gfa}
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


rule gaf2fa_2:
    input:
        gaf=rules.palss_refine_3.output.gaf,
    output:
        fa=pjoin(
            WD,
            "n{n}",
            "palss-{t}",
            "pangenome-augmented.d{d}.w{w}.c{c}.m{m}.fa",
        ),
    shell:
        """
        cut -f1,17 {input.gaf} | sed "s/^/>/" | sed "s/\\tqs:Z:/\\n/g" > {output.fa}
        """
