"""
  * pangenome: pjoin(WD, "n{n}", "pangenome-{t}.gbz")
  * reads: pjoin(WD, sample + "-reads.ec.fa"),
"""


rule hifiasm:
    input:
        fq=pjoin(WD, "cov{cov}", sample + "-reads.fq.gz"),
    output:
        fa1=pjoin(WD, "cov{cov}", sample + ".asm.bp.hap1.p_ctg.fa"),
        fa2=pjoin(WD, "cov{cov}", sample + ".asm.bp.hap2.p_ctg.fa"),
        fa=pjoin(WD, "cov{cov}", sample + ".asm.bp.haps.p_ctg.fa"),
    params:
        prefix=pjoin(WD, "cov{cov}", sample + ".asm"),
    threads: workflow.cores
    conda:
        "../envs/hifiasm.yaml"
    log:
        time=pjoin(WD, "times", "cov{cov}", "hifiasm.time"),
    shell:
        """
        /usr/bin/time -vo {log.time} hifiasm -o {params.prefix} -t{threads} {input.fq}
        awk '/^S/{{print ">"$2;print $3}}' {params.prefix}.bp.hap1.p_ctg.gfa > {output.fa1}
        awk '/^S/{{print ">"$2;print $3}}' {params.prefix}.bp.hap2.p_ctg.gfa > {output.fa2}
        cat {output.fa1} {output.fa2} > {output.fa}
        """


# rule align_mgc_contigs:
#     input:
#         fa=FA,
#         faq=pjoin(WD, "cov{cov}",sample + ".asm.bp.hap{h}.p_ctg.fa"),
#     output:
#         bam=pjoin(WD,"cov{cov}", sample + ".asm.bp.hap{h}.p_ctg.bam"),
#     conda:
#         "../envs/minimap2.yaml"
#     threads: workflow.cores
#     shell:
#         """
#         minimap2 -t{threads} --MD -ax asm5 --eqx {input.fa} {input.faq} | samtools view -bS | samtools sort > {output.bam}
#         samtools index {output.bam}
#         """


# rule install_minigraphcactus:
#     output:
#         pjoin(WD, "mgc-src", "cactus_env", "bin", "activate"),
#     params:
#         d=pjoin(WD, "mgc-src"),
#     # conda:
#     #     "../envs/mgc.yaml"
#     shell:
#         """
#         rm -rf {params.d}
#         git clone --recursive https://github.com/ComparativeGenomicsToolkit/cactus.git {params.d}
#         cd {params.d}
#         git checkout v3.1.4
#         virtualenv -p python3 cactus_env
#         echo "export PATH=$(pwd)/bin:\$PATH" >> cactus_env/bin/activate
#         echo "export PYTHONPATH=$(pwd)/lib:\$PYTHONPATH" >> cactus_env/bin/activate
#         set +u; source cactus_env/bin/activate; set -u
#         python3 -m pip install -U setuptools pip wheel
#         python3 -m pip install -U .
#         python3 -m pip install -U -r ./toil-requirement.txt
#         sed -i "s/base_singularity_call += \['-u', /base_singularity_call += \[/g" {params.d}/cactus_env/lib/python3.12/site-packages/cactus/shared/common.py
#         """


rule minigraphcactus:
    input:
        ref=FA,
        gbz=pjoin(WD, "n{n}", "pangenome-oneout.gbz"),
        fa1=rules.hifiasm.output.fa1,
        fa2=rules.hifiasm.output.fa2,
        venv=cactus_activate,
    output:
        gfa=pjoin(WD, "n{n}", "cov{cov}", "pangenome-mgcactus.gfa"),
    params:
        prefix=pjoin("/scratch2", "luca-palss", WD[1:], "n{n}", "cov{cov}", "mgcactus"),
    threads: workflow.cores
    # conda:
    #     "../envs/mgc.yaml"
    log:
        time=pjoin(WD, "times", "cov{cov}", "n{n}", "mgcactus.time"),
    shell:
        """
        set +u; source {input.venv}; set -u
        /usr/bin/time -vo {log.time} bash ./utils/run_mgcactus.sh {input.ref} {input.gbz} {sample} {input.fa1} {input.fa2} {params.prefix} {threads} | vg mod --unchop - > {output.gfa}
        """
