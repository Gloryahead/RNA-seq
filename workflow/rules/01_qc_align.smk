"""
Rules: QC, trimming, STAR alignment, featureCounts
Tutorials: Parts 1, 2
Environment: rnaseq_env (01_rnaseq_core.yml)
"""

OUTDIR = config["outdir"]
LOGDIR = config["logdir"]

# ── FastQC on raw reads ───────────────────────────────────────────────
rule fastqc_raw:
    input:
        r1 = get_r1,
        r2 = get_r2,
    output:
        html_r1 = f"{OUTDIR}/qc/fastqc/{{sample}}_R1_fastqc.html",
        html_r2 = f"{OUTDIR}/qc/fastqc/{{sample}}_R2_fastqc.html",
        zip_r1  = f"{OUTDIR}/qc/fastqc/{{sample}}_R1_fastqc.zip",
        zip_r2  = f"{OUTDIR}/qc/fastqc/{{sample}}_R2_fastqc.zip",
    log:   f"{LOGDIR}/fastqc/{{sample}}.log"
    threads: 64
    resources:
        cpus_per_task   = 64,
        runtime         = 15,
        slurm_partition = "standard",
    params:
        activate = mamba_activate("rnaseq_env"),
        outdir   = f"{OUTDIR}/qc/fastqc",
    shell:
        """
        set -eo pipefail
        {params.activate}
        fastqc --threads {threads} --outdir {params.outdir} {input.r1} {input.r2} 2>{log}
        base_r1=$(basename {input.r1} .fastq.gz)
        base_r2=$(basename {input.r2} .fastq.gz)
        mv {params.outdir}/${{base_r1}}_fastqc.html {output.html_r1}
        mv {params.outdir}/${{base_r2}}_fastqc.html {output.html_r2}
        mv {params.outdir}/${{base_r1}}_fastqc.zip  {output.zip_r1}
        mv {params.outdir}/${{base_r2}}_fastqc.zip  {output.zip_r2}
        """


# ── Trim Galore (adapter trim + embedded FastQC) ─────────────────────
rule trim_galore:
    input:
        r1 = get_r1,
        r2 = get_r2,
    output:
        r1     = temp(f"{OUTDIR}/trimmed/{{sample}}_val_1.fq.gz"),
        r2     = temp(f"{OUTDIR}/trimmed/{{sample}}_val_2.fq.gz"),
        report = f"{OUTDIR}/qc/trimming/{{sample}}_trimming_report.txt",
    log:   f"{LOGDIR}/trim_galore/{{sample}}.log"
    threads: 64
    resources:
        cpus_per_task   = 64,
        runtime         = 30,
        slurm_partition = "standard",
    params:
        activate   = mamba_activate("rnaseq_env"),
        quality    = config["trimming"]["quality"],
        min_length = config["trimming"]["min_length"],
        outdir     = f"{OUTDIR}/trimmed",
    shell:
        """
        set -eo pipefail
        {params.activate}
        trim_galore \
            --quality {params.quality} \
            --length {params.min_length} \
            --fastqc \
            --paired \
            --cores {threads} \
            --output_dir {params.outdir} \
            {input.r1} {input.r2} \
            2>{log}
        base_r1=$(basename {input.r1} .fastq.gz)
        base_r2=$(basename {input.r2} .fastq.gz)
        mv {params.outdir}/${{base_r1}}_val_1.fq.gz {output.r1}
        mv {params.outdir}/${{base_r2}}_val_2.fq.gz {output.r2}
        cat {params.outdir}/${{base_r1}}.fastq.gz_trimming_report.txt \
            {params.outdir}/${{base_r2}}.fastq.gz_trimming_report.txt \
            > {output.report} 2>>{log}
        """


# ── STAR alignment ────────────────────────────────────────────────────
rule star_align:
    input:
        r1    = f"{OUTDIR}/trimmed/{{sample}}_val_1.fq.gz",
        r2    = f"{OUTDIR}/trimmed/{{sample}}_val_2.fq.gz",
        index = config["genome"]["star_index"],
    output:
        bam      = f"{OUTDIR}/bam/{{sample}}_Aligned.sortedByCoord.out.bam",
        log      = f"{OUTDIR}/bam/{{sample}}_Log.final.out",
        chimeric = f"{OUTDIR}/bam/{{sample}}_Chimeric.out.junction",
    log:   f"{LOGDIR}/star/{{sample}}.log"
    threads: 64
    resources:
        cpus_per_task   = 64,
        runtime         = 90,
        slurm_partition = "standard",
    params:
        activate = mamba_activate("rnaseq_env"),
        extra    = config["star"]["extra"],
        prefix   = f"{OUTDIR}/bam/{{sample}}_",
        gtf      = config["genome"]["gtf"],
    shell:
        """
        set -eo pipefail
        {params.activate}
        STAR \
            --genomeDir {input.index} \
            --readFilesIn {input.r1} {input.r2} \
            --readFilesCommand zcat \
            --runThreadN {threads} \
            --outFileNamePrefix {params.prefix} \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI AS NM \
            --chimOutType Junctions \
            --chimSegmentMin 10 \
            {params.extra} \
            2>{log}
        """


# ── Index BAM ────────────────────────────────────────────────────────
rule samtools_index:
    input:  f"{OUTDIR}/bam/{{sample}}_Aligned.sortedByCoord.out.bam"
    output: f"{OUTDIR}/bam/{{sample}}_Aligned.sortedByCoord.out.bam.bai"
    log:    f"{LOGDIR}/samtools/{{sample}}_index.log"
    threads: 64
    resources: cpus_per_task=64, runtime=15
    params:
        activate = mamba_activate("rnaseq_env"),
    shell:
        """
        set -eo pipefail
        {params.activate}
        samtools index -@ {threads} {input} 2>{log}
        """


# ── featureCounts (gene-level count matrix) ───────────────────────────
rule featurecounts:
    input:
        bams = expand(f"{OUTDIR}/bam/{{sample}}_Aligned.sortedByCoord.out.bam",
                      sample=SAMPLES),
        gtf  = config["genome"]["gtf"],
    output:
        counts  = f"{OUTDIR}/counts/counts_raw.tsv",
        summary = f"{OUTDIR}/counts/counts_raw.tsv.summary",
    log:   f"{LOGDIR}/featurecounts/featurecounts.log"
    threads: 64
    resources:
        cpus_per_task   = 64,
        runtime         = 30,
        slurm_partition = "standard",
    params:
        activate    = mamba_activate("rnaseq_env"),
        strand      = config["library"]["strandedness"],
        feature     = config["featurecounts"]["feature_type"],
        attr        = config["featurecounts"]["attribute"],
        paired_flag = "-p" if config["library"]["paired"] else "",
    shell:
        """
        set -eo pipefail
        {params.activate}
        featureCounts \
            -T {threads} \
            -t {params.feature} \
            -g {params.attr} \
            -s {params.strand} \
            {params.paired_flag} \
            -a {input.gtf} \
            -o {output.counts} \
            {input.bams} \
            2>{log}
        """


# ── Salmon (transcript-level quantification) ──────────────────────────
rule salmon_quant:
    input:
        r1    = f"{OUTDIR}/trimmed/{{sample}}_val_1.fq.gz",
        r2    = f"{OUTDIR}/trimmed/{{sample}}_val_2.fq.gz",
        index = config["genome"]["salmon_index"],
    output:
        quant = f"{OUTDIR}/salmon/{{sample}}/quant.sf",
        dir   = directory(f"{OUTDIR}/salmon/{{sample}}"),
    log:   f"{LOGDIR}/salmon/{{sample}}.log"
    threads: 64
    resources: cpus_per_task=64, runtime=30
    params:
        activate = mamba_activate("rnaseq_env"),
        lib_type = config["salmon"]["lib_type"],
    shell:
        """
        set -eo pipefail
        {params.activate}
        salmon quant \
            --index {input.index} \
            --libType {params.lib_type} \
            -1 {input.r1} -2 {input.r2} \
            --threads {threads} \
            --output {output.dir} \
            --validateMappings \
            --gcBias \
            2>{log}
        """


# ── MultiQC (aggregate all QC) ────────────────────────────────────────
rule multiqc:
    input:
        expand(f"{OUTDIR}/qc/fastqc/{{sample}}_R1_fastqc.zip", sample=SAMPLES),
        expand(f"{OUTDIR}/bam/{{sample}}_Log.final.out",        sample=SAMPLES),
        f"{OUTDIR}/counts/counts_raw.tsv.summary",
    output: f"{OUTDIR}/qc/multiqc_report.html"
    log:    f"{LOGDIR}/multiqc/multiqc.log"
    resources: runtime=30
    params:
        activate = mamba_activate("rnaseq_env"),
        indir    = OUTDIR,
        outdir   = f"{OUTDIR}/qc",
    shell:
        """
        set -eo pipefail
        {params.activate}
        multiqc {params.indir} \
            --outdir {params.outdir} \
            --filename multiqc_report.html \
            --force \
            2>{log}
        """
