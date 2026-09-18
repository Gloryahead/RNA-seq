"""
Rules: Alternative splicing with SplAdder (Part 8)
Environment: spladder_env (04_spladder.yml)
Input: BAM files from STAR alignment (rule star_align)
"""

OUTDIR = config["outdir"]
LOGDIR = config["logdir"]


rule spladder_build:
    input:
        bams = expand(f"{OUTDIR}/bam/{{sample}}_Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        gtf  = config["genome"]["gtf"],
    output:
        graph = directory(f"{OUTDIR}/splicing/spladder/graph"),
    log:   f"{LOGDIR}/spladder/build.log"
    threads: 8
    resources: mem_mb=32000, runtime=180
    params:
        activate = mamba_activate("spladder_env"),
        bam_list = lambda wc, input: ",".join(input.bams),
    shell:
        """
        set -eo pipefail
        {params.activate}
        spladder build \
            --bams {params.bam_list} \
            --gtf  {input.gtf} \
            --outdir {output.graph} \
            --threads {threads} \
            --event-types exon_skip,intron_retention,alt_3prime,alt_5prime,mult_exon_skip \
            --readlen {config[library][read_length]} \
            2>{log}
        """


rule spladder_test:
    input:
        graph = f"{OUTDIR}/splicing/spladder/graph",
    output:
        tsv = f"{OUTDIR}/splicing/spladder/{{sample}}_events.tsv",
    log:   f"{LOGDIR}/spladder/test_{{sample}}.log"
    threads: 4
    resources: mem_mb=16000, runtime=120
    params:
        activate  = mamba_activate("spladder_env"),
        treatment = lambda wc: config["comparisons"][0][0],
        control   = lambda wc: config["comparisons"][0][1],
    shell:
        """
        set -eo pipefail
        {params.activate}
        spladder test \
            --outdir {input.graph} \
            --conditionA {params.treatment} \
            --conditionB {params.control} \
            --output {output.tsv} \
            2>{log}
        """
