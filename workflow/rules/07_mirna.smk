"""
Rules: miRNA-seq and UMI-based miRNA-seq (Parts 15, 15.2)
Environment: mirna_env (09_mirna.yml)
Tools: miRDeep2 (Bowtie v1), umi_tools (Part 15.2; install manually in mirna_env)
Input: single-end FASTQ
"""

OUTDIR = config["outdir"]
LOGDIR = config["logdir"]

MIRNA_ADAPTER = config["mirna"]["adapter"]
MIRNA_MIN_LEN = config["mirna"]["min_length"]
MIRNA_MAX_LEN = config["mirna"]["max_length"]


rule umi_extract:
    input:
        r1 = get_r1,
    output:
        r1_umi = temp(f"{OUTDIR}/trimmed/{{sample}}_umi_extracted.fastq.gz"),
    log:   f"{LOGDIR}/umi/{{sample}}_extract.log"
    resources: mem_mb=4000, runtime=30
    params:
        activate = mamba_activate("mirna_env"),
        pattern  = config["umi"]["pattern"],
    shell:
        """
        set -eo pipefail
        {params.activate}
        if [ "{config[umi][enabled]}" = "True" ]; then
            umi_tools extract \
                --stdin      {input.r1} \
                --stdout     {output.r1_umi} \
                --extract-method=string \
                --bc-pattern {params.pattern} \
                --log        {log}
        else
            cp {input.r1} {output.r1_umi}
        fi
        """


rule trim_mirna:
    input:
        r1 = f"{OUTDIR}/trimmed/{{sample}}_umi_extracted.fastq.gz",
    output:
        r1 = temp(f"{OUTDIR}/trimmed/{{sample}}_mirna_trimmed.fastq.gz"),
    log:   f"{LOGDIR}/trim_mirna/{{sample}}.log"
    threads: 4
    resources: mem_mb=4000, runtime=30
    params:
        activate = mamba_activate("mirna_env"),
    shell:
        """
        set -eo pipefail
        {params.activate}
        cutadapt \
            -a {MIRNA_ADAPTER} \
            --minimum-length {MIRNA_MIN_LEN} \
            --maximum-length {MIRNA_MAX_LEN} \
            --cores {threads} \
            -o {output.r1} \
            {input.r1} \
            2>{log}
        """


rule mirdeep2_mapper:
    input:
        r1    = f"{OUTDIR}/trimmed/{{sample}}_mirna_trimmed.fastq.gz",
        fasta = config["genome"]["fasta"],
    output:
        reads        = f"{OUTDIR}/mirna/{{sample}}_collapsed.fa",
        arf          = f"{OUTDIR}/mirna/{{sample}}_reads_vs_genome.arf",
        bowtie_index = temp(directory(f"{OUTDIR}/mirna/bowtie_idx_{{sample}}")),
    log:   f"{LOGDIR}/mirdeep2/{{sample}}_mapper.log"
    threads: 8
    resources: mem_mb=16000, runtime=60
    params:
        activate = mamba_activate("mirna_env"),
    shell:
        """
        set -eo pipefail
        {params.activate}
        mkdir -p {output.bowtie_index}
        bowtie-build {input.fasta} {output.bowtie_index}/genome 2>>{log}
        mapper.pl \
            <(zcat {input.r1}) \
            -e -h -j -l {MIRNA_MIN_LEN} \
            -m -q \
            -p {output.bowtie_index}/genome \
            -s {output.reads} \
            -t {output.arf} \
            -o {threads} \
            2>>{log}
        """


rule mirdeep2_quantifier:
    input:
        reads   = f"{OUTDIR}/mirna/{{sample}}_collapsed.fa",
        arf     = f"{OUTDIR}/mirna/{{sample}}_reads_vs_genome.arf",
        mature  = config["mirna"]["mirbase_mature"],
        hairpin = config["mirna"]["mirbase_hairpin"],
        fasta   = config["genome"]["fasta"],
    output:
        html = f"{OUTDIR}/mirna/{{sample}}_miRDeep2_results.html",
        tsv  = f"{OUTDIR}/mirna/{{sample}}_miRDeep2_results.tsv",
    log:   f"{LOGDIR}/mirdeep2/{{sample}}_quantifier.log"
    threads: 4
    resources: mem_mb=16000, runtime=90
    params:
        activate = mamba_activate("mirna_env"),
        outdir   = f"{OUTDIR}/mirna",
    shell:
        """
        set -eo pipefail
        {params.activate}
        cd {params.outdir}
        miRDeep2.pl \
            {input.reads} \
            {input.fasta} \
            {input.arf} \
            {input.mature} \
            none \
            {input.hairpin} \
            -P -v \
            2>{log}
        mv result_*.html {wildcards.sample}_miRDeep2_results.html 2>>{log} || true
        mv result_*.tsv  {wildcards.sample}_miRDeep2_results.tsv  2>>{log} || true
        """
