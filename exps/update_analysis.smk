from os.path import join as pjoin
import gzip
import random

seed = 42
random.seed(seed)
K = 27  # kmer size
L = 512  # max vertex size
coverage = 7.5  # coverage per haplotype

Ns = [1, 2, 4, 8, 16, 32]

FA = config["fa"]
VCF = config["vcf"]
WD = config["wd"]
REALFQ = config["fq"]  # for sampling based simulation

sample = ""
SAMPLES = {}
for line in gzip.open(VCF, mode="rt"):
    if line.startswith("#CHROM"):
        line = line.strip("\n").split("\t")[10:]  # remove hg38 from the list
        nh = 1 + 1 + (len(line) - 1) * 2
        random.shuffle(line)
        sample = line[0]
        line = line[1:]
        random.shuffle(line)
        for n in Ns:
            SAMPLES[str(n)] = line[:n]
        break

print(sample, SAMPLES)


wildcard_constraints:
    x="(1out)|(full)",


rule run:
    input:
        expand(
            pjoin(WD, "{n}", "alignments-{x}-augmented.k27.w{w}.pkl"),
            n=SAMPLES.keys(),
            x=["1out", "full"],
            w=[2], # , 3, 5],
        ),
        expand(
            pjoin(WD, "{n}", "alignments-{x}.pkl"),
            n=SAMPLES.keys(),
            x=["1out", "full"],
        ),


# ======================= #
# === READ SIMULATION === #
# ======================= #
rule get_sample_vcf:
    input:
        vcf=VCF,
    output:
        vcf=pjoin(WD, f"{sample}.vcf.gz"),
    conda:
        "./envs/bcftools.yml"
    shell:
        """
        bcftools view -s {sample} {input.vcf} | bcftools view -c 1 -Oz > {output.vcf}
        sleep 5
        tabix -p vcf {output.vcf}
        """


# rule get_haplotype:
#     input:
#         fa=FA,
#         vcf=rules.get_sample_vcf.output.vcf,
#     output:
#         fa=pjoin(WD, sample + ".hap{h}.fasta"),
#     conda:
#         "./envs/bcftools.yml"
#     shell:
#         """
#         bcftools consensus -s {sample} -H {wildcards.h} --fasta-ref {input.fa} {input.vcf} > {output.fa}
#         samtools faidx {output.fa}
#         """


rule vg_construct_sample:
    input:
        fa=FA,
        vcf=rules.get_sample_vcf.output.vcf,
    output:
        vg=pjoin(WD, sample, "pangenome.vg"),
    params:
        prefix=pjoin(WD, sample, "pangenome"),
    threads: workflow.cores / 2
    conda:
        "./envs/vg.yml"
    shell:
        """
        vg construct --threads {threads} --reference {input.fa} --alt-paths --node-max {L} --vcf {input.vcf} > {params.prefix}.walts.vg
        vg gbwt --num-jobs {threads} --discard-overlaps --vcf-input {input.vcf} --xg-name {params.prefix}.walts.vg --output {params.prefix}.haplotypes.gbwt
        vg paths --drop-paths --variant-paths -x {params.prefix}.walts.vg > {params.prefix}.ref.vg
        vg paths --extract-gam --gbwt {params.prefix}.haplotypes.gbwt -x {params.prefix}.ref.vg > {params.prefix}.haplotypes.gam
        vg augment --label-paths {params.prefix}.ref.vg {params.prefix}.haplotypes.gam | vg mod --remove-non-path - > {output.vg}
        """


rule get_haplotype:
    input:
        vg=rules.vg_construct_sample.output.vg,
    output:
        fa=pjoin(WD, sample + ".hap{h}.fasta"),
    params:
        h=lambda wildcards: int(wildcards.h) - 1,
    conda:
        "./envs/vg.yml"
    shell:
        """
        vg paths --extract-fasta --xg {input.vg} --paths-by "{sample}#{params.h}" > {output.fa}
        """


rule pbsim3:
    input:
        fa=rules.get_haplotype.output.fa,
        fq=REALFQ,
    output:
        fq=pjoin(WD, "pbsim3", "hap{h}_0001.fastq"),
    params:
        oprefix=pjoin(WD, "pbsim3", "hap{h}"),
        cov=coverage,
    threads: workflow.cores / 2
    conda:
        "./envs/pbsim3.yml"
    shell:
        """
        pbsim --id-prefix "S{wildcards.h}_" --prefix {params.oprefix} --strategy wgs --method sample --sample {input.fq} --depth {params.cov} --genome {input.fa}
        """


rule combine:
    input:
        expand(pjoin(WD, "pbsim3", "hap{h}_0001.fastq"), h=[1, 2]),
    output:
        fq=pjoin(WD, sample + ".fq"),
    params:
        oprefix=pjoin(WD, "pbsim3"),
    shell:
        """
        cat {params.oprefix}/*.fastq > {output.fq}
        """


# ============================== #
# === PANGENOME CONSTRUCTION === #
# ============================== #


rule get_fullvcf:
    input:
        vcf=VCF,
    output:
        vcf=pjoin(WD, "{n}", "variations-full.vcf.gz"),
    params:
        idxs=lambda wildcards: ",".join(SAMPLES[wildcards.n]) + "," + sample,
    conda:
        "./envs/bcftools.yml"
    shell:
        """
        bcftools view -Ou -s {params.idxs} {input.vcf} | bcftools view -Oz -c1 > {output.vcf}
        sleep 3
        tabix -p vcf {output.vcf}
        """


rule get_reducedvcf:
    input:
        vcf=VCF,
    output:
        vcf=pjoin(WD, "{n}", "variations-1out.vcf.gz"),
    params:
        idxs=lambda wildcards: ",".join(SAMPLES[wildcards.n]),
    conda:
        "./envs/bcftools.yml"
    shell:
        """
        bcftools view -Ou -s {params.idxs} {input.vcf} | bcftools view -Oz -c1 > {output.vcf}
        sleep 5
        tabix -p vcf {output.vcf}
        """


# rule get_reducedvcf:
#    input:
#        vcf=pjoin(WD, "{n}", "variations-full.vcf.gz"),
#        vcf2=pjoin(WD, f"{sample}.vcf.gz"),
#    output:
#        vcf=pjoin(WD, "{n}", "variations-1out.vcf.gz"),
#    params:
#        idxs = lambda wildcards: ",".join(SAMPLES[wildcards.n]),
#    conda:
#        "./envs/bcftools.yml"
#    shell:
#        """
#        python3 scripts/remove_sample.py {input.vcf} {input.vcf2} {sample} | bcftools view -Oz -s ^{sample} > {output.vcf}
#        sleep 1
#        tabix -p vcf {output.vcf}
#        """


rule vg_construct:
    input:
        fa=FA,
        vcf=pjoin(WD, "{n}", "variations-{x}.vcf.gz"),
    output:
        vg=pjoin(WD, "{n}", "vg-{x}", "pangenome.walts.vg"),
    threads: workflow.cores / 2
    conda:
        "./envs/vg.yml"
    shell:
        """
        vg construct --threads {threads} --reference {input.fa} --alt-paths --node-max {L} --vcf {input.vcf} > {output.vg}
        """


rule vg_droppaths:
    input:
        vg=rules.vg_construct.output.vg,
    output:
        vg=pjoin(WD, "{n}", "vg-{x}", "pangenome.ref.vg"),
    threads: workflow.cores / 2
    conda:
        "./envs/vg.yml"
    shell:
        """
        vg paths --drop-paths --variant-paths -x {input.vg} > {output.vg}
        """


rule vg_gbwt:
    input:
        vcf=pjoin(WD, "{n}", "variations-{x}.vcf.gz"),
        vg=rules.vg_construct.output.vg,
    output:
        gbwt=pjoin(WD, "{n}", "vg-{x}", "haplotypes.gbwt"),
    threads: workflow.cores
    conda:
        "./envs/vg.yml"
    shell:
        """
        vg gbwt --num-jobs {threads} --discard-overlaps --vcf-input {input.vcf} --xg-name {input.vg} --output {output.gbwt}
        """


rule vg_extractgam:
    input:
        gbwt=rules.vg_gbwt.output.gbwt,
        vg=rules.vg_droppaths.output.vg,
    output:
        gam=pjoin(WD, "{n}", "vg-{x}", "haplotypes.gam"),
    threads: workflow.cores / 2
    conda:
        "./envs/vg.yml"
    shell:
        """
        vg paths --extract-gam --gbwt {input.gbwt} -x {input.vg} > {output.gam}
        """


rule vg_augment:
    input:
        vg=rules.vg_droppaths.output.vg,
        gam=rules.vg_extractgam.output.gam,
    output:
        vg=pjoin(WD, "{n}", "pangenome-{x}.vg"),
    threads: workflow.cores / 2
    conda:
        "./envs/vg.yml"
    shell:
        """
        vg augment --label-paths {input.vg} {input.gam} | vg mod --remove-non-path - > {output.vg}
        """


rule vg_view:
    input:
        vg=rules.vg_augment.output.vg,
    output:
        gfa=pjoin(WD, "{n}", "pangenome-{x}.gfa"),
    threads: 1
    conda:
        "./envs/vg.yml"
    shell:
        """
        vg view {input.vg} > {output.gfa}
        """


# ============== #
# === METHOD === #
# ============== #


rule remove_ns:
    input:
        fq=rules.combine.output.fq,
    output:
        fq=pjoin(WD, sample + ".clean.fq"),
    threads: workflow.cores / 4
    conda:
        "./envs/biopython.yml"
    shell:
        """
        python3 scripts/remove_n.py {input.fq} > {output.fq}
        """


rule hifiasm:
    input:
        fq=rules.remove_ns.output.fq,
    output:
        fa=pjoin(WD, sample + ".ec.fa"),
    params:
        prefix=pjoin(WD, sample),
    threads: workflow.cores
    conda:
        "./envs/hifiasm.yml"
    log:
        time=pjoin(WD, "hifiasm.time"),
    shell:
        """
        /usr/bin/time -vo {log.time} hifiasm -t{threads} --write-ec --bin-only {input.fq} -o {params.prefix}
        """


rule get_paths:
    input:
        vg=pjoin(WD, "{n}", "pangenome-{x}.vg"),
    output:
        fa=pjoin(WD, "{n}", "pangenome-{x}.paths.fa"),
    threads: workflow.cores / 2
    conda:
        "./envs/vg.yml"
    shell:
        """
        vg paths --extract-fasta --xg {input.vg} > {output.fa}
        """


rule index_paths:
    input:
        fa=rules.get_paths.output.fa,
    output:
        fmd=pjoin(WD, "{n}", "pangenome-{x}.paths.fa.fmd"),
    log:
        time=pjoin(WD, "{n}", "times", "ropebwt3-{x}.time"),
    threads: workflow.cores / 2
    shell:
        """
        /usr/bin/time -vo {log.time} ../build/rb3-prefix/src/rb3/ropebwt3 build -m 3G -t {threads} -d {input.fa} -o {output.fmd}
        """


rule sketch:
    input:
        gfa=rules.vg_view.output.gfa,
        fmd=rules.index_paths.output.fmd,
    output:
        skt=pjoin(WD, "{n}", "pangenome-{x}.k{k}.skt"),
    params:
        nh=lambda wildcards: (
            int(wildcards.n) * 2 + 1
            if wildcards.x == "1out"
            else int(wildcards.n) * 2 + 2 + 1
        ),
    log:
        time=pjoin(WD, "{n}", "times", "sketch-{x}.k{k}.time"),
        log=pjoin(WD, "{n}", "logs", "sketch-{x}.k{k}.log"),
    threads: workflow.cores / 2
    shell:
        """
        /usr/bin/time -vo {log.time} ../pansv sketch -g{params.nh} -k{wildcards.k} -v25000 {input.gfa} {input.fmd} > {output.skt} 2> {log.log}
        """


rule search:
    input:
        gfa=rules.vg_view.output.gfa,
        skt=rules.sketch.output.skt,
        fmd=rules.index_paths.output.fmd,
        fa=rules.hifiasm.output.fa,
    output:
        txt=pjoin(WD, "{n}", "specific-{x}.k{k}.txt"),
    log:
        time=pjoin(WD, "{n}", "times", "search-{x}.k{k}.time"),
        log=pjoin(WD, "{n}", "logs", "search-{x}.k{k}.log"),
    threads: workflow.cores / 2
    shell:
        """
        /usr/bin/time -vo {log.time} ../pansv search -k{wildcards.k} {input.gfa} {input.skt} {input.fmd} {input.fa} > {output.txt} 2> {log.log}
        """


rule call:
    input:
        gfa=rules.vg_view.output.gfa,
        skt=rules.sketch.output.skt,
        txt=rules.search.output.txt,
        fa=rules.hifiasm.output.fa,
    output:
        gaf=pjoin(WD, "{n}", "specific-{x}.k{k}.w{w}.gaf"),
    log:
        time=pjoin(WD, "{n}", "times", "call-{x}.k{k}.w{w}.time"),
        log=pjoin(WD, "{n}", "logs", "call-{x}.k{k}.w{w}.log"),
    threads: workflow.cores / 2
    shell:
        """
        /usr/bin/time -vo {log.time} ../pansv call -w{wildcards.w} -k{wildcards.k} {input.gfa} {input.skt} {input.txt} {input.fa} > {output.gaf} 2> {log.log}
        """


rule augment:
    input:
        gfa=rules.vg_view.output.gfa,
        gaf=rules.call.output.gaf,
    output:
        gfa=pjoin(WD, "{n}", "pangenome-{x}-augmented.k{k}.w{w}.gfa"),
    log:
        time=pjoin(WD, "{n}", "times", "augment-{x}.k{k}.w{w}.time"),
    threads: workflow.cores / 2
    conda:
        "./envs/vg.yml"
    shell:
        """
        vg augment --min-coverage 1 --gaf {input.gfa} {input.gaf} > {output.gfa}
        """


# ================ #
# === ANALYSIS === #
# ================ #
rule minigraph:
    input:
        gfa=pjoin(WD, "{n}", "pangenome-{x}.gfa"),
        fq=rules.combine.output.fq,
    output:
        gaf=pjoin(WD, "{n}", "alignments-{x}.gaf"),
    threads: workflow.cores / 2
    conda:
        "./envs/minigraph.yml"
    shell:
        """
        minigraph -t{threads} -cx lr {input.gfa} {input.fq} > {output.gaf}
        """


rule minigraph_augmented:
    input:
        gfa=pjoin(WD, "{n}", "pangenome-{x}-augmented.k{k}.w{w}.gfa"),
        fq=rules.combine.output.fq,
    output:
        gaf=pjoin(WD, "{n}", "alignments-{x}-augmented.k{k}.w{w}.gaf"),
    threads: workflow.cores / 2
    conda:
        "./envs/minigraph.yml"
    shell:
        """
        minigraph -t{threads} -cx lr {input.gfa} {input.fq} > {output.gaf}
        """


rule analyze:
    input:
        gfa=pjoin(WD, "{n}", "pangenome-{x}.gfa"),
        gaf=pjoin(WD, "{n}", "alignments-{x}.gaf"),
    output:
        pkl=pjoin(WD, "{n}", "alignments-{x}.pkl"),
    shell:
        """
        python3 scripts/analyze.py {input.gfa} {input.gaf} -o {output.pkl}
        """


rule analyze_augmented:
    input:
        gfa=pjoin(WD, "{n}", "pangenome-{x}-augmented.k{k}.w{w}.gfa"),
        gaf=pjoin(WD, "{n}", "alignments-{x}-augmented.k{k}.w{w}.gaf"),
    output:
        pkl=pjoin(WD, "{n}", "alignments-{x}-augmented.k{k}.w{w}.pkl"),
    shell:
        """
        python3 scripts/analyze.py {input.gfa} {input.gaf} -o {output.pkl}
        """
