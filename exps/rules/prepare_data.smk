# sample is the new sample
# kept_samples is a map n->list of samples (with reference)

rule dump_samples:
    output:
        full=pjoin(WD, "n{n}", "samples-full.list"),
        oneout=pjoin(WD, "n{n}", "samples-oneout.list"),
        new=pjoin(WD, "n{n}", "new.list"),
    params:
        ss=lambda wildcards: kept_samples[wildcards.n],
    shell:
        """
        echo {sample} > {output.new}
        echo {params.ss} | tr " " "\\n" > {output.oneout}
        cat {output.oneout} {output.new} > {output.full}
        """

rule extract_haplotype:
    input:
        gbz=GBZ,
    output:
        fa=pjoin(WD, sample + "-hap{h}.fa"),
    shell:
        """
        vg paths --paths-by "{sample}#{wildcards.h}" --extract-fasta --xg {input.gbz} > {output.fa}
        """


rule extract_subgraph:
    input:
        gbz=GBZ,
        txt=pjoin(WD, "n{n}", "samples-{t}.list"),
    output:
        gfa=pjoin(WD, "n{n}", "pangenome-{t}.gfa"),
    shell:
        """
        ./utils/extract_subgraph {input.gbz} {input.txt} | vg mod --unchop - > {output.gfa}
        """

rule gfa2gbz:
    input:
        gfa=pjoin(WD, "n{n}", "pangenome-{t}.gfa"),
    output:
        gbz=pjoin(WD, "n{n}", "pangenome-{t}.gbz"),
    shell:
        """
        vg gbwt --gbz-format -g {output.gbz} --xg-name {input.gfa} --index-paths
        """





rule pbsim3:
    input:
        fa=rules.extract_haplotype.output.fa,
        fq=REALFQ,
    output:
        fq=pjoin(WD, "pbsim3", "hap{h}_0001.fq.gz"),
    params:
        oprefix=pjoin(WD, "pbsim3", "hap{h}"),
        cov=coverage,
    threads: workflow.cores / 2
    conda:
        "../envs/pbsim3.yaml"
    shell:
        """
        pbsim --id-prefix "S{wildcards.h}_" --prefix {params.oprefix} --strategy wgs --method sample --sample {input.fq} --depth {params.cov} --genome {input.fa} 2> {params.oprefix}.log
        """

rule combine:
    input:
        pjoin(WD, "pbsim3", "hap1_0001.fq.gz"),
        pjoin(WD, "pbsim3", "hap2_0001.fq.gz"),
    output:
        fq=pjoin(WD, sample + "-reads.fq.gz"),
    params:
        oprefix=pjoin(WD, "pbsim3"),
    threads:
        workflow.cores
    shell:
        """
        i=0
        for fq in {params.oprefix}/*.fq.gz
        do
            python3 ./utils/remove_n.py $fq | gzip -c > $fq.clean.gz &
            i=$((i+1))
            if (( i % {threads} == 0 )); then
                wait
            fi
        done
        wait
        #
        cat {params.oprefix}/*.clean.gz > {output.fq}
        rm {params.oprefix}/*.clean.gz
        """

rule hifiasm_ec:
    input:
        fq=pjoin(WD, sample + "-reads.fq.gz"),
    output:
        fa=pjoin(WD, sample + "-reads.ec.fa"),
    params:
        prefix=pjoin(WD, sample + "-reads"),
    threads: workflow.cores
    conda:
        "../envs/hifiasm.yaml"
    log:
        time=pjoin(WD, "times", "hifiasm-ec.time"),
        log=pjoin(WD, "logs", "hifiasm-ec.log"),
    shell:
        """
        /usr/bin/time -vo {log.time} hifiasm -t{threads} --write-ec --bin-only {input.fq} -o {params.prefix} 2> {log.log}
        """

        
rule align_hap:
    input:
        fa=FA,
        faq=pjoin(WD, sample + "-hap{h}.fa"),
    output:
        bam=pjoin(WD, sample + "-hap{h}.bam"),
    conda:
        "../envs/minimap2.yaml"
    threads:
        workflow.cores
    shell:
        """
        minimap2 -t{threads} --MD -ax asm5 --eqx {input.fa} {input.faq} | samtools view -bS | samtools sort > {output.bam}
        samtools index {output.bam}
        """
        
rule align_reads:
    input:
        fa=FA,
        fq=rules.combine.output.fq,
    output:
        bam=pjoin(WD, sample + "-reads.bam"),
    conda:
        "../envs/minimap2.yaml"
    threads:
        workflow.cores
    shell:
        """
        minimap2 -t{threads} --MD -ax map-hifi --eqx {input.fa} {input.fq} | samtools view -bS | samtools sort > {output.bam}
        samtools index {output.bam}
        """

rule align_corrected_reads:
    input:
        fa=FA,
        fq=rules.hifiasm_ec.output.fa,
    output:
        bam=pjoin(WD, sample + "-reads.ec.bam"),
    conda:
        "../envs/minimap2.yaml"
    threads:
        workflow.cores
    shell:
        """
        minimap2 -t{threads} --MD -ax map-hifi --eqx {input.fa} {input.fq} | samtools view -bS | samtools sort > {output.bam}
        samtools index {output.bam}
        """
