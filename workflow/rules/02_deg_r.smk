"""
Rules: Differential expression analysis
Tutorials: Parts 3, 4, 19, 20
Container: rnaseq_r_env (02_rnaseq_r.yml / containers/rnaseq_r.sif)
Supports DESeq2, edgeR, and limma-voom (set in config.yaml deg.method)
"""

OUTDIR = config["outdir"]
LOGDIR = config["logdir"]


rule deg_analysis:
    """
    Runs DEG analysis for one comparison (treatment vs control).
    Outputs: DEG table, PCA plot, MA plot, sample-distance heatmap.
    """
    input:
        counts  = f"{OUTDIR}/counts/counts_raw.tsv",
        samples = config["samples"],
    output:
        degs    = f"{OUTDIR}/deg/{{comp}}_DEG_results.tsv",
        pca     = f"{OUTDIR}/figures/{{comp}}_pca.pdf",
        ma      = f"{OUTDIR}/figures/{{comp}}_ma_plot.pdf",
        heatmap = f"{OUTDIR}/figures/{{comp}}_sample_distance.pdf",
        rds     = f"{OUTDIR}/deg/{{comp}}_dge_object.rds",
    log:   f"{LOGDIR}/deg/{{comp}}.log"
    conda: "../../environments/02_rnaseq_r.yml"
    container: "file://containers/rnaseq_r.sif"
    threads: 4
    resources:
        mem_mb   = 16000,
        runtime  = 60,
        slurm_partition = "standard",
    params:
        method      = config["deg"]["method"],
        min_count   = config["deg"]["min_count"],
        min_samples = config["deg"]["min_samples"],
        lfc         = config["deg"]["lfc_threshold"],
        alpha       = config["deg"]["alpha"],
        organism    = config["organism"],
        treatment   = lambda wc: wc.comp.split("_vs_")[0],
        control     = lambda wc: wc.comp.split("_vs_")[1],
        outdir      = f"{OUTDIR}",
        script      = "workflow/scripts/deg_analysis.R",
    shell:
        """
        Rscript {params.script} \
            --counts    {input.counts} \
            --samples   {input.samples} \
            --treatment {params.treatment} \
            --control   {params.control} \
            --method    {params.method} \
            --min_count {params.min_count} \
            --min_samples {params.min_samples} \
            --lfc       {params.lfc} \
            --alpha     {params.alpha} \
            --organism  {params.organism} \
            --out_degs  {output.degs} \
            --out_pca   {output.pca} \
            --out_ma    {output.ma} \
            --out_heatmap {output.heatmap} \
            --out_rds   {output.rds} \
            2>{log}
        """


rule publication_figures:
    """
    Publication-ready figures from DEG results (Part 4):
    EnhancedVolcano, pheatmap of top DEGs, PCA biplot.
    """
    input:
        degs = f"{OUTDIR}/deg/{{comp}}_DEG_results.tsv",
        rds  = f"{OUTDIR}/deg/{{comp}}_dge_object.rds",
    output:
        volcano = f"{OUTDIR}/figures/{{comp}}_volcano.pdf",
        heatmap = f"{OUTDIR}/figures/{{comp}}_top50_heatmap.pdf",
    log:   f"{LOGDIR}/figures/{{comp}}.log"
    conda: "../../environments/02_rnaseq_r.yml"
    container: "file://containers/rnaseq_r.sif"
    threads: 1
    resources: mem_mb=8000, runtime=30
    params:
        alpha      = config["deg"]["alpha"],
        lfc        = config["deg"]["lfc_threshold"],
        treatment  = lambda wc: wc.comp.split("_vs_")[0],
        control    = lambda wc: wc.comp.split("_vs_")[1],
    shell:
        """
        Rscript workflow/scripts/figures.R \
            --degs     {input.degs} \
            --rds      {input.rds} \
            --treatment {params.treatment} \
            --control   {params.control} \
            --alpha    {params.alpha} \
            --lfc      {params.lfc} \
            --out_volcano {output.volcano} \
            --out_heatmap {output.heatmap} \
            2>{log}
        """


rule normalization_report:
    """
    Normalization comparison report (Part 19):
    CPM, TPM, RPKM, TMM, VST — side-by-side density plots.
    """
    input:
        counts = f"{OUTDIR}/counts/counts_raw.tsv",
        gtf    = config["genome"]["gtf"],
    output:
        report = f"{OUTDIR}/qc/normalization_comparison.pdf",
    log:   f"{LOGDIR}/normalization/normalization.log"
    conda: "../../environments/02_rnaseq_r.yml"
    container: "file://containers/rnaseq_r.sif"
    threads: 1
    resources: mem_mb=8000, runtime=30
    shell:
        """
        Rscript workflow/scripts/normalization.R \
            --counts {input.counts} \
            --gtf    {input.gtf} \
            --output {output.report} \
            2>{log}
        """
