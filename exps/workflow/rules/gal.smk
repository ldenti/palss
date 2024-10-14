rule graphaligner:
    input:
        gfa=pjoin(WD, SAMPLE, "pansv-l{l}", "{graph}.gfa"),
        fq=pjoin(WD, SAMPLE, "hifi.fq"),
    output:
        gaf=pjoin(WD, SAMPLE, "graphaligner-l{l}", "{graph}.gaf"),
    threads: workflow.cores
    conda:
        "../envs/graphaligner.yml"
    log:
        time=pjoin(WD, SAMPLE, "TIMES", "graphaligner-l{l}", "{graph}.time"),
    shell:
        """
        /usr/bin/time -vo {log.time} GraphAligner -g {input.gfa} -f {input.fq} -a {output.gaf} -x vg --threads {threads}
        """


# rule graphchainer:
#     input:
#         gfa=pjoin(WD, SAMPLE, "pansv-l{l}", "{graph}.gfa"),
#         fq=pjoin(WD, SAMPLE, "hifi.fq"),
#     output:
#         gaf=pjoin(WD, SAMPLE, "graphchainer-l{l}", "{graph}.gaf"),
#     threads: workflow.cores
#     conda:
#         "."
#     shell:
#         """
#         GraphChainer -g {input.gfa} -f {input.fq} -a {output.gaf} --threads {threads}
#         """
