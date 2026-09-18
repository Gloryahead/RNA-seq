"""
Rules: Cancer subtype classification + GSVA + deconvolution (Parts 5–6, 11)
  - PAM50 breast cancer subtyping
  - GSVA gene set variation analysis
  - CIBERSORTx immune deconvolution (web-only; prep script only)
Environment: cancer_subtype_env (03_cancer_subtype_r.yml)
"""

OUTDIR = config["outdir"]
LOGDIR = config["logdir"]


if config["cancer_subtype"]["run_pam50"]:

    rule pam50_subtype:
        input:
            counts  = f"{OUTDIR}/counts/counts_raw.tsv",
            samples = config["samples"],
        output:
            subtypes = f"{OUTDIR}/cancer_subtype/pam50_subtypes.tsv",
            plot     = f"{OUTDIR}/cancer_subtype/pam50_heatmap.pdf",
        log:   f"{LOGDIR}/cancer_subtype/pam50.log"
        resources: mem_mb=8000, runtime=30
        params:
            activate = mamba_activate("cancer_subtype_env"),
            outdir   = f"{OUTDIR}/cancer_subtype",
        shell:
            """
            set -eo pipefail
            {params.activate}
            Rscript workflow/scripts/cancer_subtype.R \
                --mode    pam50 \
                --counts  {input.counts} \
                --samples {input.samples} \
                --outdir  {params.outdir} \
                2>{log}
            """


if config["cancer_subtype"]["run_gsva"]:

    rule gsva_analysis:
        input:
            counts  = f"{OUTDIR}/counts/counts_raw.tsv",
            samples = config["samples"],
        output:
            scores  = f"{OUTDIR}/cancer_subtype/gsva_scores.tsv",
            heatmap = f"{OUTDIR}/cancer_subtype/gsva_heatmap.pdf",
        log:   f"{LOGDIR}/cancer_subtype/gsva.log"
        threads: 4
        resources: mem_mb=16000, runtime=60
        params:
            activate = mamba_activate("cancer_subtype_env"),
            organism = config["organism"],
            outdir   = f"{OUTDIR}/cancer_subtype",
        shell:
            """
            set -eo pipefail
            {params.activate}
            Rscript workflow/scripts/cancer_subtype.R \
                --mode    gsva \
                --counts  {input.counts} \
                --samples {input.samples} \
                --organism {params.organism} \
                --outdir  {params.outdir} \
                2>{log}
            """


if config["cancer_subtype"]["run_cibersortx_prep"]:

    rule cibersortx_prep:
        """Prepare mixture matrix for CIBERSORTx upload (cibersortx.stanford.edu)."""
        input:
            counts  = f"{OUTDIR}/counts/counts_raw.tsv",
            samples = config["samples"],
        output:
            mixture = f"{OUTDIR}/cancer_subtype/cibersortx_mixture.tsv",
        log:   f"{LOGDIR}/cancer_subtype/cibersortx_prep.log"
        resources: mem_mb=8000, runtime=15
        params:
            activate = mamba_activate("cancer_subtype_env"),
            outdir   = f"{OUTDIR}/cancer_subtype",
        shell:
            """
            set -eo pipefail
            {params.activate}
            Rscript workflow/scripts/cancer_subtype.R \
                --mode    cibersortx_prep \
                --counts  {input.counts} \
                --samples {input.samples} \
                --outdir  {params.outdir} \
                2>{log}
            echo "Upload {output.mixture} to https://cibersortx.stanford.edu" >> {log}
            """
