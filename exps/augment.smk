from os.path import join as pjoin
import gzip
import random

seed = 42 # 23
random.seed(seed)
K = 27  # kmer size
coverage = 7.5  # coverage per haplotype

Ns = [1, 8, 32]  # , 2, 4, 8, 16, 32]

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
    x="(1out)|(full)", # |(mgcactus)",


rule run:
    input:
        expand(
            pjoin(WD, "{n}", "alignments-{x}-augmented.k27.w{w}.s{s}.gaf"),
            n=SAMPLES.keys(),
            x=["1out", "full"],
            w=[2], # , 3, 5],
            s=[0],
        ),
        expand(
            pjoin(WD, "{n}", "alignments-{x}.gaf"),
            n=SAMPLES.keys(),
            x=["1out", "full"], #, "mgcactus"],
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
        tabix -p vcf {output.vcf}
        """


# we cannot use bcftools consensus since it filters variations differently from vg construct


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
        vg construct --threads {threads} --reference {input.fa} --alt-paths --node-max 512 --vcf {input.vcf} > {params.prefix}-walts.vg
        vg gbwt --num-jobs {threads} --discard-overlaps --vcf-input {input.vcf} --xg-name {params.prefix}-walts.vg --output {params.prefix}-haplotypes.gbwt
        vg paths --drop-paths --variant-paths -x {params.prefix}-walts.vg > {params.prefix}-ref.vg
        vg paths --extract-gam --gbwt {params.prefix}-haplotypes.gbwt -x {params.prefix}-ref.vg > {params.prefix}-haplotypes.gam
        vg augment --label-paths {params.prefix}-ref.vg {params.prefix}-haplotypes.gam | vg mod --remove-non-path - > {output.vg}
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
    conda:
        "./envs/biopython.yml"
    shell:
        """
        cat {params.oprefix}/*.fastq | python3 ./scripts/remove_n.py > {output.fq}
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


rule build_pangenome:
    input:
        fa=FA,
        vcf=pjoin(WD, "{n}", "variations-{x}.vcf.gz"),
    output:
        vg=pjoin(WD, "{n}", "pangenome-{x}.vg"),
        gfa=pjoin(WD, "{n}", "pangenome-{x}.gfa"),
    params:
        prefix=pjoin(WD, "{n}", "vg-{x}", "pangenome"),
    threads: workflow.cores
    conda:
        "./envs/vg.yml"
    shell:
        """
        vg construct --threads {threads} --reference {input.fa} --alt-paths --node-max 512 --vcf {input.vcf} > {params.prefix}-walts.vg
        vg paths --drop-paths --variant-paths -x {params.prefix}-walts.vg > {params.prefix}-ref.vg
        vg gbwt --num-jobs {threads} --discard-overlaps --vcf-input {input.vcf} --xg-name {params.prefix}-walts.vg --output {params.prefix}-haplotypes.gbwt
        vg paths --extract-gam --gbwt {params.prefix}-haplotypes.gbwt -x {params.prefix}-ref.vg > {params.prefix}-haplotypes.gam
        vg augment --label-paths {params.prefix}-ref.vg {params.prefix}-haplotypes.gam | vg mod --remove-non-path - > {output.vg}
        vg view {output.vg} > {output.gfa}
        """


# ============== #
# === METHOD === #
# ============== #
rule hifiasm_ec:
    input:
        fq=rules.combine.output.fq,
    output:
        fa=pjoin(WD, sample + ".ec.fa"),
    params:
        prefix=pjoin(WD, sample),
    threads: workflow.cores
    conda:
        "./envs/hifiasm.yml"
    log:
        time=pjoin(WD, "hifiasm_ec.time"),
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
        gfa=pjoin(WD, "{n}", "pangenome-{x}.gfa"),
        fmd=rules.index_paths.output.fmd,
    output:
        skt=pjoin(WD, "{n}", "pangenome-{x}.gfa.k{k}.skt"),
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
        /usr/bin/time -vo {log.time} ../palss sketch -g{params.nh} -k{wildcards.k} {input.gfa} {input.fmd} > {output.skt} 2> {log.log}
        """


rule search:
    input:
        gfa=pjoin(WD, "{n}", "pangenome-{x}.gfa"),
        skt=rules.sketch.output.skt,
        fmd=rules.index_paths.output.fmd,
        fa=rules.hifiasm_ec.output.fa,
    output:
        fa=pjoin(WD, "{n}", "clusters-{x}.k{k}.w{w}.fa"),
    log:
        time=pjoin(WD, "{n}", "times", "search-{x}.k{k}.w{w}.time"),
        log=pjoin(WD, "{n}", "logs", "search-{x}.k{k}.w{w}.log"),
    threads: workflow.cores / 2
    shell:
        """
        /usr/bin/time -vo {log.time} ../palss search -@{threads} -k{wildcards.k} -w{wildcards.w} {input.gfa} {input.skt} {input.fmd} {input.fa} > {output.fa} 2> {log.log}
        """

rule GraphAligner:
    input:
        gfa=pjoin(WD, "{n}", "pangenome-{x}.gfa"),
        fa=rules.search.output.fa,
    output:
        gaf=pjoin(WD, "{n}", "clusters-{x}.k{k}.w{w}.gaf"),
    log:
        time=pjoin(WD, "{n}", "times", "graphaligner-{x}.k{k}.w{w}.time"),
    conda:
        "./envs/graphaligner.yml"
    threads: workflow.cores
    shell:
        """
        /usr/bin/time -vo {log.time} GraphAligner --graph {input.gfa} --reads {input.fa} --alignments-out {output.gaf} --preset vg --threads {threads}
        """

rule clean_graphaligner:
    input:
        gaf=rules.GraphAligner.output.gaf,
        fa=rules.search.output.fa,
    output:
        txt=pjoin(WD, "{n}", "goodclusters-{x}.k{k}.w{w}.txt"),
    log:
        time=pjoin(WD, "{n}", "times", "filtergaf-{x}.k{k}.w{w}.time"),
    conda:
        "./envs/biopython.yml"
    shell:
        """
        /usr/bin/time -vo {log.time} python3 ../scripts/filter_gaf.py {input.gaf} {input.fa} > {output.txt}
        """

rule realign:
    input:
        gfa=pjoin(WD, "{n}", "pangenome-{x}.gfa"),
        txt=rules.clean_graphaligner.output.txt,
    output:
        gaf=pjoin(WD, "{n}", "goodclusters-{x}.k{k}.w{w}.s{s}.gaf"),
    log:
        time=pjoin(WD, "{n}", "times", "call-{x}.k{k}.w{w}.s{s}.time"),
        log=pjoin(WD, "{n}", "logs", "call-{x}.k{k}.w{w}.s{s}.log"),
    threads: workflow.cores / 2
    shell:
        """
        /usr/bin/time -vo {log.time} ../palss realign -s{wildcards.s} {input.gfa} {input.txt} > {output.gaf} 2> {log.log}
        """

rule augment:
    input:
        gfa=pjoin(WD, "{n}", "pangenome-{x}.gfa"),
        gaf=rules.realign.output.gaf,
    output:
        gfa=pjoin(WD, "{n}", "pangenome-{x}-augmented.k{k}.w{w}.s{s}.gfa"),
    log:
        time=pjoin(WD, "{n}", "times", "augment-{x}.k{k}.w{w}.s{s}.time"),
    threads: workflow.cores / 2
    conda:
        "./envs/vg.yml"
    shell:
        """
        vg augment --min-coverage 1 --gaf {input.gfa} {input.gaf} > {output.gfa}.tmp
        python3 ../scripts/clean_augmented_gfa.py {output.gfa}.tmp > {output.gfa}
        """

# ======================== #
# === MINIGRAPH-CACTUS === #
# ======================== #
rule hifiasm:
    input:
        fq=rules.combine.output.fq,
    output:
        fa1=pjoin(WD, sample + ".asm.bp.hap1.p_ctg.fa"),
        fa2=pjoin(WD, sample + ".asm.bp.hap2.p_ctg.fa"),
    params:
        prefix=pjoin(WD, sample + ".asm"),
    threads: workflow.cores
    conda:
        "./envs/hifiasm.yml"
    log:
        time=pjoin(WD, "hifiasm.time"),
    shell:
        """
        /usr/bin/time -vo {log.time} hifiasm -o {params.prefix} -t{threads} {input.fq}
        awk '/^S/{{print ">"$2;print $3}}' {params.prefix}.bp.hap1.p_ctg.gfa > {output.fa1}
        awk '/^S/{{print ">"$2;print $3}}' {params.prefix}.bp.hap2.p_ctg.gfa > {output.fa2}
        """

rule install_minigraphcactus:
    output:
        pjoin(WD, "mgc-src", "cactus_env", "bin", "activate"),
    params:
        d=pjoin(WD, "mgc-src"),
    shell:
        """
        rm -rf {params.d}
        git clone --recursive https://github.com/ComparativeGenomicsToolkit/cactus.git {params.d}
        cd {params.d}
        virtualenv -p python3 cactus_env
        echo "export PATH=$(pwd)/bin:\$PATH" >> cactus_env/bin/activate
        echo "export PYTHONPATH=$(pwd)/lib:\$PYTHONPATH" >> cactus_env/bin/activate
        set +u; source cactus_env/bin/activate; set -u
        python3 -m pip install -U setuptools pip wheel
        python3 -m pip install -U .
        python3 -m pip install -U -r ./toil-requirement.txt
        """

rule minigraphcactus:
    input:
        ref=FA,
        vcf=pjoin(WD, "{n}", "variations-1out.vcf.gz"),
        fa1=rules.hifiasm.output.fa1,
        fa2=rules.hifiasm.output.fa2,
        venv=pjoin(WD, "mgc-src", "cactus_env", "bin", "activate"),
    output:
        gfa=pjoin(WD, "{n}", "pangenome-mgcactus.gfa"),
    params:
        prefix = pjoin(WD, "{n}", "mgcactus"),
    threads: workflow.cores
    log:
        time = pjoin(WD, "{n}", "times", "mgcactus.time"),
        log = pjoin(WD, "{n}", "logs", "mgcactus.log"),
    shell:
        """
        set +u; source {input.venv}; set -u
        /usr/bin/time -vo {log.time} bash run_mgcactus.sh {input.ref} {input.vcf} {sample} {input.fa1} {input.fa2} {params.prefix} {threads} > {output.gfa}
        """

# ================ #
# === ANALYSIS === #
# ================ #
rule graphaligner:
    input:
        gfa=pjoin(WD, "{n}", "pangenome-{x}.gfa"),
        fq=rules.combine.output.fq,
    output:
        gaf=pjoin(WD, "{n}", "alignments-{x}.gaf"),
    threads: workflow.cores / 2
    conda:
        "./envs/graphaligner.yml"
    shell:
        """
        GraphAligner --graph {input.gfa} --reads {input.fq} --alignments-out {output.gaf} --preset vg --threads {threads}
        """

rule graphaligner_augmented:
    input:
        gfa=pjoin(WD, "{n}", "pangenome-{x}-augmented.k{k}.w{w}.s{s}.gfa"),
        fq=rules.combine.output.fq,
    output:
        gaf=pjoin(WD, "{n}", "alignments-{x}-augmented.k{k}.w{w}.s{s}.gaf"),
    threads: workflow.cores / 2
    conda:
        "./envs/graphaligner.yml"
    shell:
        """
        GraphAligner --graph {input.gfa} --reads {input.fq} --alignments-out {output.gaf} --preset vg --threads {threads}
        """