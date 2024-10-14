rule graphaligner:
    input:
        gfa=pjoin(WD, sample, "pansv-l{l}", "{graph}.gfa"),
        fq=pjoin(WD, sample, "hifi.fq"),
    output:
        gaf=pjoin(WD, sample, "graphaligner-l{l}", "{graph}.gaf"),
    threads: workflow.cores
    conda:
        "../envs/graphaligner.yml"
    log:
        time=pjoin(WD, sample, "TIMES", "graphaligner-l{l}", "{graph}.time"),
    shell:
        """
        /usr/bin/time -vo {log.time} GraphAligner -g {input.gfa} -f {input.fq} -a {output.gaf} -x vg --threads {threads}
        """


# rule graphchainer:
#     input:
#         gfa=pjoin(WD, sample, "pansv-l{l}", "{graph}.gfa"),
#         fq=pjoin(WD, sample, "hifi.fq"),
#     output:
#         gaf=pjoin(WD, sample, "graphchainer-l{l}", "{graph}.gaf"),
#     threads: workflow.cores
#     conda:
#         "."
#     shell:
#         """
#         GraphChainer -g {input.gfa} -f {input.fq} -a {output.gaf} --threads {threads}
#         """
