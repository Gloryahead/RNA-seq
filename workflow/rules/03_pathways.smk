"""
Rules: Pathway enrichment analysis
Tutorial: Part 5
Container: rnaseq_r_env (02_rnaseq_r.yml / containers/rnaseq_r.sif)
Methods: GO ORA, KEGG ORA, GSEA (MSigDB optional)
"""

OUTDIR = config["outdir"]
LOGDIR = config["logdir"]


rule pathway_enrichment:
    """
    ORA and GSEA from a DEG results table.
    Outputs: GO/KEGG/GSEA tables, dot plots, enrichment plots.
    """
    input:
        degs = f"{OUTDIR}/deg/{{comp}}_DEG_results.tsv",
    output:
        go_table   = f"{OUTDIR}/pathways/{{comp}}_GO_results.tsv"   if config["pathways"]["run_go"]   else [],
        kegg_table = f"{OUTDIR}/pathways/{{comp}}_KEGG_results.tsv" if config["pathways"]["run_kegg"] else [],
        gsea_table = f"{OUTDIR}/pathways/{{comp}}_GSEA_results.tsv" if config["pathways"]["run_gsea"] else [],
        go_plot    = f"{OUTDIR}/pathways/{{comp}}_GO_dotplot.pdf"    if config["pathways"]["run_go"]   else [],
        kegg_plot  = f"{OUTDIR}/pathways/{{comp}}_KEGG_dotplot.pdf" if config["pathways"]["run_kegg"] else [],
    log:   f"{LOGDIR}/pathways/{{comp}}.log"
    conda: "../../environments/02_rnaseq_r.yml"
    container: "file://containers/rnaseq_r.sif"
    threads: 2
    resources:
        mem_mb   = 12000,
        runtime  = 60,
        slurm_partition = "standard",
    params:
        organism     = config["organism"],
        pval         = config["pathways"]["pvalue_cutoff"],
        qval         = config["pathways"]["qvalue_cutoff"],
        msigdb_gmt   = config["pathways"]["msigdb_gmt"],
        run_go       = config["pathways"]["run_go"],
        run_kegg     = config["pathways"]["run_kegg"],
        run_gsea     = config["pathways"]["run_gsea"],
        outdir       = f"{OUTDIR}/pathways",
        label        = lambda wc: wc.comp,
    shell:
        """
        Rscript workflow/scripts/pathway_analysis.R \
            --degs       {input.degs} \
            --organism   {params.organism} \
            --pval       {params.pval} \
            --qval       {params.qval} \
            --msigdb_gmt "{params.msigdb_gmt}" \
            --run_go     {params.run_go} \
            --run_kegg   {params.run_kegg} \
            --run_gsea   {params.run_gsea} \
            --outdir     {params.outdir} \
            --label      {params.label} \
            2>{log}
        """
