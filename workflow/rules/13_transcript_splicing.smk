"""
Rules: Transcript-level quantification (StringTie) and differential
       alternative splicing (rMATS) from STAR BAMs.

StringTie workflow:
  1. stringtie_assemble  — guided assembly per sample → sample GTF
  2. stringtie_merge     — merge all sample GTFs → consensus reference GTF
  3. stringtie_quant     — re-estimate expression against merged GTF (-e -B)
     Output: gene_abund.tab (gene TPM) + t_data.ctab (Ballgown input)

rMATS workflow:
  rmats — one job per comparison; 5 event types (SE, A3SS, A5SS, MXE, RI)
  Output: {event}.MATS.JC.txt (junction counts) per comparison

Environment: transcript_splicing_env (17_transcript_splicing.yml)
"""

OUTDIR = config["outdir"]
LOGDIR = config["logdir"]


def _group_bams(group):
    """Return sorted BAM paths for all samples in the given group."""
    return sorted(
        f"{OUTDIR}/bam/{s}_Aligned.sortedByCoord.out.bam"
        for s in SAMPLES if samples_df.loc[s, "group"] == group
    )


def _strand_stringtie(strand):
    """Map featureCounts strandedness code to StringTie flag."""
    return {"1": "--rf", "2": "--fr"}.get(str(strand), "")


def _libtype_rmats(strand):
    """Map featureCounts strandedness code to rMATS libType."""
    return {"1": "fr-firststrand", "2": "fr-secondstrand"}.get(
        str(strand), "fr-unstranded"
    )


# ── StringTie: guided assembly (per sample) ───────────────────────────
rule stringtie_assemble:
    input:
        bam = f"{OUTDIR}/bam/{{sample}}_Aligned.sortedByCoord.out.bam",
        bai = f"{OUTDIR}/bam/{{sample}}_Aligned.sortedByCoord.out.bam.bai",
        gtf = config["genome"]["gtf"],
    output:
        gtf = f"{OUTDIR}/stringtie/assembly/{{sample}}.gtf",
    log:   f"{LOGDIR}/stringtie/{{sample}}_assemble.log"
    threads: 8
    resources:
        cpus_per_task   = 8,
        runtime         = 60,
        slurm_partition = "standard",
    params:
        activate = mamba_activate("transcript_splicing_env"),
        strand   = _strand_stringtie(config["library"]["strandedness"]),
        outdir   = f"{OUTDIR}/stringtie/assembly",
    shell:
        """
        set -eo pipefail
        {params.activate}
        mkdir -p {params.outdir}
        stringtie {input.bam} \
            -G {input.gtf} \
            {params.strand} \
            -o {output.gtf} \
            -p {threads} \
            2>{log}
        """


# ── StringTie: merge all sample GTFs into one consensus reference ─────
rule stringtie_merge:
    input:
        sample_gtfs = expand(f"{OUTDIR}/stringtie/assembly/{{sample}}.gtf", sample=SAMPLES),
        ref_gtf     = config["genome"]["gtf"],
    output:
        merged = f"{OUTDIR}/stringtie/merged.gtf",
    log:   f"{LOGDIR}/stringtie/merge.log"
    threads: 4
    resources:
        cpus_per_task = 4,
        runtime       = 30,
    params:
        activate = mamba_activate("transcript_splicing_env"),
    shell:
        """
        set -eo pipefail
        {params.activate}
        stringtie --merge \
            -G {input.ref_gtf} \
            -o {output.merged} \
            {input.sample_gtfs} \
            2>{log}
        """


# ── StringTie: re-estimate expression against merged reference ────────
# -e: restrict to transcripts in guide GTF (expression mode)
# -B: write Ballgown ctab files for downstream R analysis
rule stringtie_quant:
    input:
        bam    = f"{OUTDIR}/bam/{{sample}}_Aligned.sortedByCoord.out.bam",
        bai    = f"{OUTDIR}/bam/{{sample}}_Aligned.sortedByCoord.out.bam.bai",
        merged = f"{OUTDIR}/stringtie/merged.gtf",
    output:
        gtf   = f"{OUTDIR}/stringtie/quant/{{sample}}/quant.gtf",
        abund = f"{OUTDIR}/stringtie/quant/{{sample}}/gene_abund.tab",
        ctab  = f"{OUTDIR}/stringtie/quant/{{sample}}/t_data.ctab",
    log:   f"{LOGDIR}/stringtie/{{sample}}_quant.log"
    threads: 8
    resources:
        cpus_per_task   = 8,
        runtime         = 60,
        slurm_partition = "standard",
    params:
        activate = mamba_activate("transcript_splicing_env"),
        strand   = _strand_stringtie(config["library"]["strandedness"]),
        outdir   = f"{OUTDIR}/stringtie/quant/{{sample}}",
    shell:
        """
        set -eo pipefail
        {params.activate}
        mkdir -p {params.outdir}
        stringtie {input.bam} \
            -G {input.merged} \
            -e -B \
            {params.strand} \
            -o {output.gtf} \
            -A {output.abund} \
            -p {threads} \
            2>{log}
        """


# ── rMATS: differential alternative splicing per comparison ──────────
# Outputs 5 event types: SE (skipped exon), A3SS/A5SS (alt splice sites),
# MXE (mutually exclusive exons), RI (retained intron)
# JC = junction-spanning reads only; also produces JCEC (junction + exon body)
rule rmats:
    input:
        b1_bams = lambda wc: _group_bams(get_treatment(wc.comp)),
        b2_bams = lambda wc: _group_bams(get_control(wc.comp)),
        gtf     = config["genome"]["gtf"],
    output:
        se  = f"{OUTDIR}/rmats/{{comp}}/SE.MATS.JC.txt",
        a3  = f"{OUTDIR}/rmats/{{comp}}/A3SS.MATS.JC.txt",
        a5  = f"{OUTDIR}/rmats/{{comp}}/A5SS.MATS.JC.txt",
        mxe = f"{OUTDIR}/rmats/{{comp}}/MXE.MATS.JC.txt",
        ri  = f"{OUTDIR}/rmats/{{comp}}/RI.MATS.JC.txt",
    log:   f"{LOGDIR}/rmats/{{comp}}.log"
    threads: 8
    resources:
        cpus_per_task   = 8,
        runtime         = 120,
        slurm_partition = "standard",
    params:
        activate  = mamba_activate("transcript_splicing_env"),
        outdir    = f"{OUTDIR}/rmats/{{comp}}",
        tmpdir    = f"{OUTDIR}/rmats/{{comp}}_tmp",
        read_len  = config["library"]["read_length"],
        libtype   = _libtype_rmats(config["library"]["strandedness"]),
        read_type = "paired" if config["library"]["paired"] else "single",
    shell:
        """
        set -eo pipefail
        {params.activate}
        mkdir -p {params.outdir} {params.tmpdir}
        b1=$(echo {input.b1_bams} | sed 's/ /,/g')
        b2=$(echo {input.b2_bams} | sed 's/ /,/g')
        rmats.py \
            --b1 $b1 \
            --b2 $b2 \
            --gtf {input.gtf} \
            -t {params.read_type} \
            --readLength {params.read_len} \
            --libType {params.libtype} \
            --od {params.outdir} \
            --tmp {params.tmpdir} \
            --nthread {threads} \
            2>{log}
        """
