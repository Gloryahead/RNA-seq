"""
Rules: Cancer subtype classification + GSVA + deconvolution (Parts 5–6, 11)
  - PAM50 breast cancer subtyping (Part 5/6)
  - GSVA gene set variation analysis (Part 11)
  - CIBERSORTx immune deconvolution (Part 11 — web-only; prep script only)
Container: cancer_subtype_env (03_cancer_subtype_r.yml / containers/cancer_subtype.sif)
"""

OUTDIR = config["outdir"]
LOGDIR = config["logdir"]


if config["cancer_subtype"]["run_pam50"]:

    rule pam50_subtype:
        """
        PAM50 breast cancer intrinsic subtype calling (Luminal A/B, HER2, Basal, Normal).
        Requires normalized expression matrix (log2 TPM or microarray RMA).
        """
        input:
            counts  = f"{OUTDIR}/counts/counts_raw.tsv",
            samples = config["samples"],
        output:
            subtypes = f"{OUTDIR}/cancer_subtype/pam50_subtypes.tsv",
            plot     = f"{OUTDIR}/cancer_subtype/pam50_heatmap.pdf",
        log:   f"{LOGDIR}/cancer_subtype/pam50.log"
        conda: "../../environments/03_cancer_subtype_r.yml"
        container: "file://containers/cancer_subtype.sif"
        resources: mem_mb=8000, runtime=30
        params:
            outdir = f"{OUTDIR}/cancer_subtype",
        shell:
            """
            Rscript workflow/scripts/cancer_subtype.R \
                --mode    pam50 \
                --counts  {input.counts} \
                --samples {input.samples} \
                --outdir  {params.outdir} \
                2>{log}
            """


if config["cancer_subtype"]["run_gsva"]:

    rule gsva_analysis:
        """
        GSVA (Gene Set Variation Analysis) — sample-level pathway scoring.
        Uses MSigDB Hallmark + KEGG gene sets.
        """
        input:
            counts  = f"{OUTDIR}/counts/counts_raw.tsv",
            samples = config["samples"],
        output:
            scores = f"{OUTDIR}/cancer_subtype/gsva_scores.tsv",
            heatmap= f"{OUTDIR}/cancer_subtype/gsva_heatmap.pdf",
        log:   f"{LOGDIR}/cancer_subtype/gsva.log"
        conda: "../../environments/03_cancer_subtype_r.yml"
        container: "file://containers/cancer_subtype.sif"
        threads: 4
        resources: mem_mb=16000, runtime=60
        params:
            organism = config["organism"],
            outdir   = f"{OUTDIR}/cancer_subtype",
        shell:
            """
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
        """
        Prepare mixture matrix for CIBERSORTx upload.
        CIBERSORTx itself is web-only (cibersortx.stanford.edu) — this rule
        writes the input matrix in the required format.
        Upload <outdir>/cibersortx_mixture.tsv to the CIBERSORTx portal.
        """
        input:
            counts  = f"{OUTDIR}/counts/counts_raw.tsv",
            samples = config["samples"],
        output:
            mixture = f"{OUTDIR}/cancer_subtype/cibersortx_mixture.tsv",
        log:   f"{LOGDIR}/cancer_subtype/cibersortx_prep.log"
        conda: "../../environments/03_cancer_subtype_r.yml"
        container: "file://containers/cancer_subtype.sif"
        resources: mem_mb=8000, runtime=15
        params:
            outdir = f"{OUTDIR}/cancer_subtype",
        shell:
            """
            Rscript workflow/scripts/cancer_subtype.R \
                --mode    cibersortx_prep \
                --counts  {input.counts} \
                --samples {input.samples} \
                --outdir  {params.outdir} \
                2>{log}
            echo "Upload {output.mixture} to https://cibersortx.stanford.edu" >> {log}
            """
