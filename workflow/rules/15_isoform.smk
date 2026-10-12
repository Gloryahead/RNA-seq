"""
Rules: Isoform-level differential transcript usage (DTU) analysis.

Workflow:
  isoform_analysis — per comparison:
    1. Import Salmon transcript-level quant.sf files via tximport
    2. Run DRIMSeq DTU pipeline (dmDSdata → dmFilter → dmPrecision → dmFit → dmTest)
    3. Write transcript- and gene-level DTU results + proportion plots

Salmon quant.sf files are produced by the salmon_quant rule in 01_qc_align.smk
and are consumed here without re-running Salmon.

Environment: isoform_env  (environments/18_isoform_dtu.yml)
  Bioconductor packages (tximport, DRIMSeq, GenomicFeatures, BiocParallel,
  AnnotationDbi) installed via:
    Rscript environments/r_scripts/install_isoform_r.R
"""

OUTDIR = config["outdir"]
LOGDIR = config["logdir"]


def _comp_quant_files(wc):
    """Return quant.sf paths for all samples belonging to this comparison."""
    trt = get_treatment(wc.comp)
    ctl = get_control(wc.comp)
    comp_samples = [
        s for s in SAMPLES
        if samples_df.loc[s, "group"] in (trt, ctl)
    ]
    return [f"{OUTDIR}/salmon/{s}/quant.sf" for s in comp_samples]


# ── Isoform DTU: tximport + DRIMSeq per comparison ───────────────────
rule isoform_analysis:
    input:
        quant_files = _comp_quant_files,
        gtf         = config["genome"]["gtf"],
        samples     = config["samples"],
    output:
        dtu_results = f"{OUTDIR}/isoform/{{comp}}_dtu_results.tsv",
        dtu_plots   = f"{OUTDIR}/isoform/{{comp}}_dtu_plots.pdf",
    log:   f"{LOGDIR}/isoform/{{comp}}_isoform.log"
    threads: 8
    resources:
        cpus_per_task   = 8,
        runtime         = 120,
        mem_mb          = 32000,
        slurm_partition = "standard",
    params:
        activate          = mamba_activate("isoform_env"),
        script            = "workflow/scripts/isoform_analysis.R",
        treatment         = lambda wc: get_treatment(wc.comp),
        control           = lambda wc: get_control(wc.comp),
        organism          = config["organism"],
        min_feature_expr  = config["isoform"]["min_feature_expr"],
        min_gene_expr     = config["isoform"]["min_gene_expr"],
        min_feature_prop  = config["isoform"]["min_feature_prop"],
        outdir            = f"{OUTDIR}/isoform",
    shell:
        """
        set -eo pipefail
        {params.activate}
        mkdir -p {params.outdir}
        Rscript {params.script} \
            --quant_files     {input.quant_files} \
            --gtf             {input.gtf} \
            --samples         {input.samples} \
            --treatment       {params.treatment} \
            --control         {params.control} \
            --organism        {params.organism} \
            --min_feature_expr {params.min_feature_expr} \
            --min_gene_expr    {params.min_gene_expr} \
            --min_feature_prop {params.min_feature_prop} \
            --threads          {threads} \
            --out_dtu          {output.dtu_results} \
            --out_plots        {output.dtu_plots} \
            2>{log}
        """
