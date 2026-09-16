"""
Rules: WGCNA co-expression networks + GENIE3 GRN + Master Regulator analysis
Tutorials: WGCNA, GENIE3, Master Regulator tutorials
Container: wgcna_env / genie3_masterreg_env
Input: normalized count matrix from DEG pipeline
"""

OUTDIR = config["outdir"]
LOGDIR = config["logdir"]


if config["networks"]["run_wgcna"]:

    rule wgcna:
        """
        WGCNA gene co-expression network: soft threshold → module detection
        → module-trait correlation → hub gene identification → STRINGdb overlay.
        """
        input:
            counts  = f"{OUTDIR}/counts/counts_raw.tsv",
            samples = config["samples"],
        output:
            modules = f"{OUTDIR}/networks/wgcna/module_assignments.tsv",
            hubs    = f"{OUTDIR}/networks/wgcna/hub_genes.tsv",
            plot    = f"{OUTDIR}/networks/wgcna/module_trait_correlations.pdf",
        log:   f"{LOGDIR}/networks/wgcna.log"
        conda: "../../environments/15_wgcna_r.yml"
        container: "file://containers/wgcna.sif"
        threads: config["networks"]["wgcna_threads"]
        resources: mem_mb=32000, runtime=240
        params:
            organism = config["organism"],
            outdir   = f"{OUTDIR}/networks/wgcna",
        shell:
            """
            Rscript workflow/scripts/wgcna.R \
                --counts    {input.counts} \
                --samples   {input.samples} \
                --threads   {threads} \
                --organism  {params.organism} \
                --outdir    {params.outdir} \
                2>{log}
            """


if config["networks"]["run_genie3"]:

    rule genie3:
        """
        GENIE3 gene regulatory network: random forest importance → adjacency matrix
        → top regulatory links → DoRothEA TF regulon overlap.
        """
        input:
            counts  = f"{OUTDIR}/counts/counts_raw.tsv",
            samples = config["samples"],
        output:
            links   = f"{OUTDIR}/networks/genie3/regulatory_links.tsv",
            network = f"{OUTDIR}/networks/genie3/network_plot.pdf",
        log:   f"{LOGDIR}/networks/genie3.log"
        conda: "../../environments/16_genie3_masterreg_r.yml"
        container: "file://containers/genie3_masterreg.sif"
        threads: 8
        resources: mem_mb=32000, runtime=300
        params:
            organism = config["organism"],
            outdir   = f"{OUTDIR}/networks/genie3",
        shell:
            """
            Rscript workflow/scripts/genie3.R \
                --counts   {input.counts} \
                --samples  {input.samples} \
                --organism {params.organism} \
                --outdir   {params.outdir} \
                2>{log}
            """


if config["networks"]["run_master_regulator"]:

    rule master_regulator:
        """
        Master regulator analysis: RegEnrich + RTN.
        Identifies top transcription factors driving observed expression changes.
        """
        input:
            counts  = f"{OUTDIR}/counts/counts_raw.tsv",
            samples = config["samples"],
            degs    = expand(f"{OUTDIR}/deg/{{comp}}_DEG_results.tsv", comp=COMPARISONS),
        output:
            mr_table = f"{OUTDIR}/networks/masterreg/master_regulators.tsv",
            plot     = f"{OUTDIR}/networks/masterreg/master_regulator_network.pdf",
        log:   f"{LOGDIR}/networks/masterreg.log"
        conda: "../../environments/16_genie3_masterreg_r.yml"
        container: "file://containers/genie3_masterreg.sif"
        threads: 4
        resources: mem_mb=24000, runtime=180
        params:
            organism = config["organism"],
            outdir   = f"{OUTDIR}/networks/masterreg",
        shell:
            """
            Rscript workflow/scripts/master_regulator.R \
                --counts   {input.counts} \
                --samples  {input.samples} \
                --degs     {input.degs[0]} \
                --organism {params.organism} \
                --outdir   {params.outdir} \
                2>{log}
            """
