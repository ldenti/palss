from os.path import join as pjoin
import random

random.seed(23)
import gzip

FA = config["fa"]
VCF = config["vcf"]
WD = config["wd"]

nths = workflow.cores

Ns = [2]  # [1, 2]  # , 4] # , 8, 16, 32]
Ks = [27]  # [23, 27, 31]

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


rule run:
    input:
        expand(pjoin(WD, "{n}", "k{k}", "reference-anchors.bed"), n=Ns, k=Ks),
        # expand(pjoin(WD, "k{k}", "sample.nuk.png"), k=Ks),


# ============================== #
# === PANGENOME CONSTRUCTION === #
# ============================== #


rule get_vcf:
    input:
        vcf=VCF,
    output:
        vcf=pjoin(WD, "{n}", "variations.vcf.gz"),
    params:
        idxs=lambda wildcards: ",".join(SAMPLES[wildcards.n]),
    conda:
        "./envs/bcftools.yml"
    shell:
        """
        bcftools view -Ou -s {params.idxs} {input.vcf} | bcftools view -Oz -c1 > {output.vcf}
        sleep 3
        tabix -p vcf {output.vcf}
        """


rule vg_construct:
    input:
        fa=FA,
        vcf=rules.get_vcf.output.vcf,
    output:
        vg=pjoin(WD, "{n}", "vg", "pangenome.walts.vg"),
    threads: workflow.cores
    conda:
        "./envs/vg.yml"
    shell:
        """
        vg construct --threads {threads} --reference {input.fa} --alt-paths --node-max 512 --vcf {input.vcf} > {output.vg}
        """


rule vg_droppaths:
    input:
        vg=rules.vg_construct.output.vg,
    output:
        vg=pjoin(WD, "{n}", "vg", "pangenome.ref.vg"),
    threads: workflow.cores / 2
    conda:
        "./envs/vg.yml"
    shell:
        """
        vg paths --drop-paths --variant-paths -x {input.vg} > {output.vg}
        """


rule vg_gbwt:
    input:
        vcf=rules.get_vcf.output.vcf,
        vg=rules.vg_construct.output.vg,
    output:
        gbwt=pjoin(WD, "{n}", "vg", "haplotypes.gbwt"),
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
        gam=pjoin(WD, "{n}", "vg", "haplotypes.gam"),
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
        vg=pjoin(WD, "{n}", "pangenome.vg"),
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
        gfa=pjoin(WD, "{n}", "pangenome.gfa"),
    threads: 1
    conda:
        "./envs/vg.yml"
    shell:
        """
        vg view {input.vg} > {output.gfa}
        """


rule get_paths:
    input:
        gfa=rules.vg_view.output.gfa,
    output:
        fa=pjoin(WD, "{n}", "paths.fa"),
    shell:
        """
        vg paths --extract-fasta --xg {input.gfa} > {output.fa}
        """


rule rb3_index:
    input:
        fa=rules.get_paths.output.fa,
    output:
        fmd=pjoin(WD, "{n}", "pangenome.paths.fa.fmd"),
    log:
        time=pjoin(WD, "times", "{n}", "rb3-build.time"),
    threads: workflow.cores
    shell:
        """
        /usr/bin/time -vo {log.time} ../build/rb3-prefix/src/rb3/ropebwt3 build -m 2G -t {threads} -d {input.fa} > {output.fmd}
        """


rule sketch:
    input:
        gfa=rules.vg_view.output.gfa,
        fmd=rules.rb3_index.output.fmd,
    output:
        skt=pjoin(WD, "{n}", "k{k}", "pangenome.skt"),
    params:
        nh=lambda wildcards: int(wildcards.n) * 2 + 1,
    log:
        time=pjoin(WD, "times", "{n}", "k{k}", "sketch.time"),
        log=pjoin(WD, "times", "{n}", "k{k}", "sketch.log"),
    threads: workflow.cores
    shell:
        """
        /usr/bin/time -vo {log.time} ../pansv sketch -g{params.nh} -k {wildcards.k} -v25000 -@{threads} {input.gfa} {input.fmd} > {output.skt} 2> {log.log}
        """


rule kan_ref:
    input:
        skt=rules.sketch.output.skt,
        fa=FA,
    output:
        bed=pjoin(WD, "{n}", "k{k}", "reference-anchors.bed"),
    log:
        time=pjoin(WD, "times", "{n}", "k{k}", "kan-reference.time"),
    shell:
        """
        /usr/bin/time -vo {log.time} ../pansv kan -k{wildcards.k} {input.skt} {input.fa} > {output.bed}
        """


# rule kan_py:
#     input:
#         bed=rules.kan_ref.output.bed,
#         fai=FA + ".fai",
#     output:
#         pjoin(WD, "k{k}", "reference.kan.png"),
#     params:
#         prefix=pjoin(WD, "k{k}", "reference.kan"),
#     conda: "workflow/envs/seaborn.yml"
#     shell:
#         """
#         python3 ../scripts/kan_hist.py {input.bed} {input.fai} -o {params.prefix}
#         """


# rule kan_reads:
#     input:
#         skt=rules.sketch.output.skt,
#         fq=FQ,
#     output:
#         nuk=pjoin(WD, "k{k}", "sample.nuk"),
#     log:
#         time = pjoin(WD, "TIMES", "k{k}", "kan-reads.time"),
#     shell:
#         """
#         /usr/bin/time -vo {log.time} ../pansv chreads -k{wildcards.k} {input.skt} {input.fq} > {output.nuk}
#         """
# rule kan_reads_py:
#     input:
#         nuk=rules.kan_reads.output.nuk,
#     output:
#         png=pjoin(WD, "k{k}", "sample.nuk.png"),
#         txt=pjoin(WD, "k{k}", "sample.nuk.txt"),
#     conda: "workflow/envs/seaborn.yml"
#     shell:
#         """
#         python3 ../scripts/ran_hist.py {input.nuk} {output.png} > {output.txt}
#         """
