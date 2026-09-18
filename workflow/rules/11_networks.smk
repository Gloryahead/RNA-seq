"""
Rules: WGCNA co-expression networks + GENIE3 GRN + Master Regulator analysis
Environments: wgcna_env (15_wgcna_r.yml) / genie3_masterreg_env (16_genie3_masterreg_r.yml)
Input: normalized count matrix from DEG pipeline
"""

OUTDIR = config["outdir"]
LOGDIR = config["logdir"]


if config["networks"]["run_wgcna"]:

    rule wgcna:
        input:
            counts  = f"{OUTDIR}/counts/counts_raw.tsv",
            samples = config["samples"],
        output:
            modules = f"{OUTDIR}/networks/wgcna/module_assignments.tsv",
            hubs    = f"{OUTDIR}/networks/wgcna/hub_genes.tsv",
            plot    = f"{OUTDIR}/networks/wgcna/module_trait_correlations.pdf",
        log:   f"{LOGDIR}/networks/wgcna.log"
        threads: config["networks"]["wgcna_threads"]
        resources: mem_mb=32000, runtime=240
        params:
            activate = mamba_activate("wgcna_env"),
            organism = config["organism"],
            outdir   = f"{OUTDIR}/networks/wgcna",
        shell:
            """
            set -eo pipefail
            {params.activate}
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
        input:
            counts  = f"{OUTDIR}/counts/counts_raw.tsv",
            samples = config["samples"],
        output:
            links   = f"{OUTDIR}/networks/genie3/regulatory_links.tsv",
            network = f"{OUTDIR}/networks/genie3/network_plot.pdf",
        log:   f"{LOGDIR}/networks/genie3.log"
        threads: 8
        resources: mem_mb=32000, runtime=300
        params:
            activate = mamba_activate("genie3_masterreg_env"),
            organism = config["organism"],
            outdir   = f"{OUTDIR}/networks/genie3",
        shell:
            """
            set -eo pipefail
            {params.activate}
            Rscript workflow/scripts/genie3.R \
                --counts   {input.counts} \
                --samples  {input.samples} \
                --organism {params.organism} \
                --outdir   {params.outdir} \
                2>{log}
            """


if config["networks"]["run_master_regulator"]:

    rule master_regulator:
        input:
            counts  = f"{OUTDIR}/counts/counts_raw.tsv",
            samples = config["samples"],
            degs    = expand(f"{OUTDIR}/deg/{{comp}}_DEG_results.tsv", comp=COMPARISONS),
        output:
            mr_table = f"{OUTDIR}/networks/masterreg/master_regulators.tsv",
            plot     = f"{OUTDIR}/networks/masterreg/master_regulator_network.pdf",
        log:   f"{LOGDIR}/networks/masterreg.log"
        threads: 4
        resources: mem_mb=24000, runtime=180
        params:
            activate = mamba_activate("genie3_masterreg_env"),
            organism = config["organism"],
            outdir   = f"{OUTDIR}/networks/masterreg",
        shell:
            """
            set -eo pipefail
            {params.activate}
            Rscript workflow/scripts/master_regulator.R \
                --counts   {input.counts} \
                --samples  {input.samples} \
                --degs     {input.degs[0]} \
                --organism {params.organism} \
                --outdir   {params.outdir} \
                2>{log}
            """
