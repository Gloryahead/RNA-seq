"""
Rules: RNA editing detection (A-to-I / adenosine-to-inosine) via JACUSA2.

Workflow:
  1. jacusa2_call1   — per-sample site detection (call-1, single-condition)
  2. jacusa2_call2   — per-comparison differential editing (call-2, two-condition)
  3. rna_edit_analysis — R downstream: filtering, A>G signature plots, summary

Environments:
  - JACUSA2 rules: rnaediting_env  (environments/05_reditools2.yml)
  - R analysis:    rnaseq_r_env    (environments/02_rnaseq_r.yml)
    (R packages are NOT installed in rnaediting_env — see env YAML comment)
"""

OUTDIR = config["outdir"]
LOGDIR = config["logdir"]


def _editing_group_bams(group):
    """Return sorted BAM paths for all samples in the given group."""
    return sorted(
        f"{OUTDIR}/bam/{s}_Aligned.sortedByCoord.out.bam"
        for s in SAMPLES if samples_df.loc[s, "group"] == group
    )


def _comp_call1_tsvs(wc):
    """Return call-1 TSVs for every sample belonging to this comparison."""
    trt = get_treatment(wc.comp)
    ctl = get_control(wc.comp)
    comp_samples = [
        s for s in SAMPLES
        if samples_df.loc[s, "group"] in (trt, ctl)
    ]
    return [f"{OUTDIR}/rna_editing/call1/{s}_sites.tsv" for s in comp_samples]


# ── JACUSA2 call-1: per-sample editing site detection ────────────────
rule jacusa2_call1:
    input:
        bam = f"{OUTDIR}/bam/{{sample}}_Aligned.sortedByCoord.out.bam",
        bai = f"{OUTDIR}/bam/{{sample}}_Aligned.sortedByCoord.out.bam.bai",
    output:
        tsv = f"{OUTDIR}/rna_editing/call1/{{sample}}_sites.tsv",
    log:   f"{LOGDIR}/rna_editing/{{sample}}_call1.log"
    threads: 8
    resources:
        cpus_per_task   = 8,
        runtime         = 60,
        slurm_partition = "standard",
    params:
        activate = mamba_activate("rnaediting_env"),
        outdir   = f"{OUTDIR}/rna_editing/call1",
    shell:
        """
        set -eo pipefail
        {params.activate}
        mkdir -p {params.outdir}
        JACUSA2 call-1 \
            -r {output.tsv} \
            -p {threads} \
            -m 20 \
            -q 20 \
            -c 10 \
            {input.bam} \
            2>{log}
        """


# ── JACUSA2 call-2: per-comparison differential editing ──────────────
rule jacusa2_call2:
    input:
        trt_bams = lambda wc: _editing_group_bams(get_treatment(wc.comp)),
        ctl_bams = lambda wc: _editing_group_bams(get_control(wc.comp)),
    output:
        tsv = f"{OUTDIR}/rna_editing/call2/{{comp}}_differential.tsv",
    log:   f"{LOGDIR}/rna_editing/{{comp}}_call2.log"
    threads: 16
    resources:
        cpus_per_task   = 16,
        runtime         = 120,
        slurm_partition = "standard",
    params:
        activate = mamba_activate("rnaediting_env"),
        outdir   = f"{OUTDIR}/rna_editing/call2",
    shell:
        """
        set -eo pipefail
        {params.activate}
        mkdir -p {params.outdir}
        trt=$(echo {input.trt_bams} | sed 's/ /,/g')
        ctl=$(echo {input.ctl_bams} | sed 's/ /,/g')
        JACUSA2 call-2 \
            -r {output.tsv} \
            -p {threads} \
            -m 20 \
            -q 20 \
            -c 10 \
            $trt \
            $ctl \
            2>{log}
        """


# ── RNA editing R analysis: filter, A>G signature, summary ───────────
rule rna_edit_analysis:
    input:
        call2    = f"{OUTDIR}/rna_editing/call2/{{comp}}_differential.tsv",
        call1    = _comp_call1_tsvs,
        samples  = config["samples"],
    output:
        table   = f"{OUTDIR}/rna_editing/{{comp}}_editing_filtered.tsv",
        plots   = f"{OUTDIR}/rna_editing/{{comp}}_editing_plots.pdf",
        summary = f"{OUTDIR}/rna_editing/{{comp}}_editing_summary.tsv",
    log:   f"{LOGDIR}/rna_editing/{{comp}}_analysis.log"
    threads: 4
    resources:
        cpus_per_task   = 4,
        runtime         = 30,
        slurm_partition = "standard",
    params:
        activate      = mamba_activate("rnaseq_r_env"),
        script        = "workflow/scripts/rna_editing.R",
        call1_dir     = f"{OUTDIR}/rna_editing/call1",
        treatment     = lambda wc: get_treatment(wc.comp),
        control       = lambda wc: get_control(wc.comp),
        min_coverage  = config["rna_editing"]["min_coverage"],
        min_frequency = config["rna_editing"]["min_frequency"],
        min_edited    = config["rna_editing"]["min_edited_reads"],
        known_snps    = config["rna_editing"]["known_snps"],
    shell:
        """
        set -eo pipefail
        {params.activate}
        Rscript {params.script} \
            --call2      {input.call2} \
            --call1_dir  {params.call1_dir} \
            --samples    {input.samples} \
            --treatment  {params.treatment} \
            --control    {params.control} \
            --min_coverage   {params.min_coverage} \
            --min_frequency  {params.min_frequency} \
            --min_edited_reads {params.min_edited} \
            --known_snps "{params.known_snps}" \
            --out_table   {output.table} \
            --out_plots   {output.plots} \
            --out_summary {output.summary} \
            2>{log}
        """
