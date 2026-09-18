"""
Rules: circRNA detection (Parts 13, 13.2)
  - CIRCexplorer2 (Part 13): uses STAR chimeric junctions → circrna_env
  - CIRI3 (Part 13.2): uses BWA alignment → ciri3_env
Input: trimmed FASTQs from trim_galore
"""

OUTDIR = config["outdir"]
LOGDIR = config["logdir"]


# ── CIRCexplorer2 (Part 13) ───────────────────────────────────────────
if config["circrna"]["run_circexplorer2"]:

    rule circexplorer2_parse:
        input:
            chimeric = f"{OUTDIR}/bam/{{sample}}_Chimeric.out.junction",
        output:
            bed = f"{OUTDIR}/circrna/circexplorer2/{{sample}}_junction.bed",
        log:   f"{LOGDIR}/circexplorer2/{{sample}}_parse.log"
        resources: mem_mb=8000, runtime=30
        params:
            activate = mamba_activate("circrna_env"),
        shell:
            """
            set -eo pipefail
            {params.activate}
            CIRCexplorer2 parse \
                -t STAR \
                {input.chimeric} \
                > {output.bed} \
                2>{log}
            """

    rule circexplorer2_annotate:
        input:
            bed    = f"{OUTDIR}/circrna/circexplorer2/{{sample}}_junction.bed",
            fasta  = config["genome"]["fasta"],
            refann = config["circrna"]["annotation_refseq"],
        output:
            txt = f"{OUTDIR}/circrna/circexplorer2/{{sample}}_circularRNA_known.txt",
        log:   f"{LOGDIR}/circexplorer2/{{sample}}_annotate.log"
        resources: mem_mb=8000, runtime=30
        params:
            activate = mamba_activate("circrna_env"),
        shell:
            """
            set -eo pipefail
            {params.activate}
            CIRCexplorer2 annotate \
                -r {input.refann} \
                -g {input.fasta} \
                {input.bed} \
                -o {output.txt} \
                2>{log}
            """


# ── CIRI3 (Part 13.2) ─────────────────────────────────────────────────
if config["circrna"]["run_ciri3"]:

    rule bwa_align_ciri3:
        input:
            r1    = f"{OUTDIR}/trimmed/{{sample}}_val_1.fq.gz",
            r2    = f"{OUTDIR}/trimmed/{{sample}}_val_2.fq.gz",
            fasta = config["genome"]["fasta"],
        output:
            bam = f"{OUTDIR}/bam_bwa/{{sample}}.bam",
        log:   f"{LOGDIR}/bwa/{{sample}}.log"
        threads: 16
        resources: mem_mb=32000, runtime=120
        params:
            activate = mamba_activate("ciri3_env"),
        shell:
            """
            set -eo pipefail
            {params.activate}
            bwa mem -T 19 -t {threads} \
                {input.fasta} {input.r1} {input.r2} \
                | samtools sort -@ 4 -o {output.bam} \
                2>{log}
            samtools index {output.bam}
            """

    rule ciri3_detect:
        input:
            bam   = f"{OUTDIR}/bam_bwa/{{sample}}.bam",
            fasta = config["genome"]["fasta"],
            gtf   = config["genome"]["gtf"],
        output:
            txt = f"{OUTDIR}/circrna/ciri3/{{sample}}_ciri3.txt",
        log:   f"{LOGDIR}/ciri3/{{sample}}.log"
        threads: 8
        resources: mem_mb=16000, runtime=90
        params:
            activate = mamba_activate("ciri3_env"),
            jar      = config["circrna"]["ciri3_jar"],
        shell:
            """
            set -eo pipefail
            {params.activate}
            java -jar {params.jar} \
                -bam  {input.bam} \
                -ref  {input.fasta} \
                -gtf  {input.gtf} \
                -out  {output.txt} \
                -nthread {threads} \
                2>{log}
            """

    rule ciri3_merge:
        input:
            txts = expand(f"{OUTDIR}/circrna/ciri3/{{sample}}_ciri3.txt", sample=SAMPLES),
        output:
            bsj = f"{OUTDIR}/circrna/ciri3/all_samples.BSJ_Matrix.txt",
            fsj = f"{OUTDIR}/circrna/ciri3/all_samples.FSJ_Matrix.txt",
        log:   f"{LOGDIR}/ciri3/merge.log"
        params:
            activate = mamba_activate("ciri3_env"),
            jar      = config["circrna"]["ciri3_jar"],
            indir    = f"{OUTDIR}/circrna/ciri3",
        shell:
            """
            set -eo pipefail
            {params.activate}
            java -jar {params.jar} merge \
                -indir  {params.indir} \
                -outbsj {output.bsj} \
                -outfsj {output.fsj} \
                2>{log}
            """
