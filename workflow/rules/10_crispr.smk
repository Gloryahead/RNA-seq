"""
Rules: CRISPR screen analysis with MAGeCK (Part 20)
Container: mageck_env (14_mageck.yml / containers/mageck.sif)
Input: FASTQ files from sequencing of sgRNA libraries
"""

OUTDIR = config["outdir"]
LOGDIR = config["logdir"]


rule mageck_count:
    """Count sgRNA reads per sample from FASTQ."""
    input:
        r1      = get_r1,
        library = config["crispr"]["library"],
    output:
        count  = f"{OUTDIR}/crispr/counts/{{sample}}.count.txt",
        log_   = f"{OUTDIR}/crispr/counts/{{sample}}.count_normalized.txt",
    log:   f"{LOGDIR}/mageck/count_{{sample}}.log"
    conda: "../../environments/14_mageck.yml"
    container: "file://containers/mageck.sif"
    threads: 4
    resources: mem_mb=8000, runtime=60
    params:
        outpfx    = f"{OUTDIR}/crispr/counts/{{sample}}",
        pam       = config["crispr"].get("pam", "TTTN"),
        sgrna_len = config["crispr"].get("sgrna_length", 20),
    shell:
        """
        mageck count \
            -l {input.library} \
            --fastq {input.r1} \
            -n {params.outpfx} \
            --sgrna-len {params.sgrna_len} \
            --pam {params.pam} \
            2>{log}
        """


rule mageck_merge_counts:
    """Merge per-sample count files into a single count matrix."""
    input:
        counts = expand(f"{OUTDIR}/crispr/counts/{{sample}}.count.txt", sample=SAMPLES),
    output:
        matrix = f"{OUTDIR}/crispr/all_samples.count_matrix.txt",
    log: f"{LOGDIR}/mageck/merge_counts.log"
    run:
        import pandas as pd
        frames = []
        for f, s in zip(input.counts, SAMPLES):
            df = pd.read_csv(f, sep="\t", index_col=0, usecols=[0, 1])
            df.columns = [s]
            frames.append(df)
        merged = pd.concat(frames, axis=1)
        merged.to_csv(output.matrix, sep="\t")


rule mageck_test:
    """MAGeCK MLE/RRA test for enriched/depleted sgRNAs per comparison."""
    input:
        matrix = f"{OUTDIR}/crispr/all_samples.count_matrix.txt",
    output:
        gene_summary   = f"{OUTDIR}/crispr/results/{{comp}}.gene_summary.txt",
        sgrna_summary  = f"{OUTDIR}/crispr/results/{{comp}}.sgrna_summary.txt",
    log:   f"{LOGDIR}/mageck/test_{{comp}}.log"
    conda: "../../environments/14_mageck.yml"
    container: "file://containers/mageck.sif"
    threads: 4
    resources: mem_mb=8000, runtime=60
    params:
        outpfx    = f"{OUTDIR}/crispr/results/{{comp}}",
        treatment = lambda wc: ",".join(
            samples_df[samples_df.group == wc.comp.split("_vs_")[0]].index.tolist()
        ),
        control   = lambda wc: ",".join(
            samples_df[samples_df.group == wc.comp.split("_vs_")[1]].index.tolist()
        ),
        norm      = config["crispr"].get("normalization", "median"),
    shell:
        """
        mageck test \
            -k {input.matrix} \
            -t {params.treatment} \
            -c {params.control} \
            -n {params.outpfx} \
            --norm-method {params.norm} \
            2>{log}
        """


rule mageck_plot:
    """MAGeCK R-based visualization: volcano, GSEA enrichment plots."""
    input:
        gene_summary = expand(
            f"{OUTDIR}/crispr/results/{{comp}}.gene_summary.txt",
            comp=COMPARISONS,
        ),
    output:
        plots = directory(f"{OUTDIR}/crispr/figures"),
    log:   f"{LOGDIR}/mageck/plot.log"
    conda: "../../environments/14_mageck.yml"
    container: "file://containers/mageck.sif"
    resources: mem_mb=8000, runtime=30
    params:
        indir    = f"{OUTDIR}/crispr/results",
        outdir   = f"{OUTDIR}/crispr/figures",
        organism = config["organism"],
    shell:
        """
        Rscript workflow/scripts/mageck_plots.R \
            --indir    {params.indir} \
            --outdir   {params.outdir} \
            --organism {params.organism} \
            2>{log}
        """
