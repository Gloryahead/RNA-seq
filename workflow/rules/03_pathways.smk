"""
Rules: Pathway enrichment analysis
Tutorial: Part 5
Environment: rnaseq_r_env (02_rnaseq_r.yml)
Methods: GO ORA (ALL/UP/DOWN × ALL/BP/CC/MF), KEGG ORA, GSEA optional
Output structure: {outdir}/pathways/{comp}/GO-ALL|GO-UP|GO-DOWN/
"""

OUTDIR = config["outdir"]
LOGDIR = config["logdir"]


rule pathway_enrichment:
    input:
        degs = f"{OUTDIR}/deg/{{comp}}_DEG_results.tsv",
    output:
        go_all_enrich  = f"{OUTDIR}/pathways/{{comp}}/GO-ALL/output-GO-Enrichment.txt"  if config["pathways"]["run_go"]   else [],
        go_up_enrich   = f"{OUTDIR}/pathways/{{comp}}/GO-UP/output-GO-Enrichment.txt"   if config["pathways"]["run_go"]   else [],
        go_down_enrich = f"{OUTDIR}/pathways/{{comp}}/GO-DOWN/output-GO-Enrichment.txt" if config["pathways"]["run_go"]   else [],
        kegg_all       = f"{OUTDIR}/pathways/{{comp}}/GO-ALL/output-KEGG.txt"           if config["pathways"]["run_kegg"] else [],
        kegg_up        = f"{OUTDIR}/pathways/{{comp}}/GO-UP/output-KEGG.txt"            if config["pathways"]["run_kegg"] else [],
        kegg_down      = f"{OUTDIR}/pathways/{{comp}}/GO-DOWN/output-KEGG.txt"          if config["pathways"]["run_kegg"] else [],
    log:   f"{LOGDIR}/pathways/{{comp}}.log"
    threads: 2
    resources:
        mem_mb   = 12000,
        runtime  = 60,
        slurm_partition = "standard",
    params:
        activate   = mamba_activate("rnaseq_r_env"),
        organism   = config["organism"],
        pval       = config["pathways"]["pvalue_cutoff"],
        qval       = config["pathways"]["qvalue_cutoff"],
        msigdb_gmt = config["pathways"]["msigdb_gmt"],
        run_go     = config["pathways"]["run_go"],
        run_kegg   = config["pathways"]["run_kegg"],
        run_gsea   = config["pathways"]["run_gsea"],
        outdir     = lambda wc: f"{OUTDIR}/pathways/{wc.comp}",
        label      = lambda wc: wc.comp,
    shell:
        """
        set -eo pipefail
        {params.activate}
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
