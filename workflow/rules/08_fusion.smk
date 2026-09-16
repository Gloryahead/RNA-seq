"""
Rules: Fusion gene detection (Parts 16–17)
  - STAR-Fusion (Part 16): uses STAR chimeric junctions
  - FusionCatcher (Part 17): end-to-end pipeline from FASTQ
Container: starfusion_env / fusioncatcher_env
"""

OUTDIR = config["outdir"]
LOGDIR = config["logdir"]


if config["fusion"]["run_starfusion"]:

    rule star_fusion:
        """
        STAR-Fusion from chimeric junctions produced by STAR.
        Requires CTAT genome library (~30 GB); set fusion.ctat_lib in config.
        Note: Docker/Apptainer preferred over conda for production use.
        """
        input:
            chimeric = f"{OUTDIR}/bam/{{sample}}_Chimeric.out.junction",
        output:
            tsv = f"{OUTDIR}/fusion/starfusion/{{sample}}/star-fusion.fusion_predictions.tsv",
        log:   f"{LOGDIR}/starfusion/{{sample}}.log"
        conda: "../../environments/10_starfusion.yml"
        container: "file://containers/starfusion.sif"
        threads: 8
        resources: mem_mb=64000, runtime=120
        params:
            ctat_lib = config["fusion"]["ctat_lib"],
            outdir   = f"{OUTDIR}/fusion/starfusion/{{sample}}",
        shell:
            """
            STAR-Fusion \
                --genome_lib_dir {params.ctat_lib} \
                -J {input.chimeric} \
                --output_dir {params.outdir} \
                --CPU {threads} \
                --examine_coding_effect \
                2>{log}
            """

    rule merge_starfusion:
        """Merge per-sample STAR-Fusion predictions into a cohort TSV."""
        input:
            tsvs = expand(
                f"{OUTDIR}/fusion/starfusion/{{sample}}/star-fusion.fusion_predictions.tsv",
                sample=SAMPLES,
            ),
        output:
            merged = f"{OUTDIR}/fusion/starfusion/all_samples_fusions.tsv",
        log: f"{LOGDIR}/starfusion/merge.log"
        run:
            import pandas as pd
            frames = []
            for tsv, s in zip(input.tsvs, SAMPLES):
                df = pd.read_csv(tsv, sep="\t", comment="#")
                df.insert(0, "sample", s)
                frames.append(df)
            pd.concat(frames).to_csv(output.merged, sep="\t", index=False)


if config["fusion"]["run_fusioncatcher"]:

    rule fusioncatcher:
        """
        FusionCatcher end-to-end detection from raw paired-end FASTQs.
        Human data only for production; set fusion.fusioncatcher_data in config.
        """
        input:
            r1 = get_r1,
            r2 = get_r2,
        output:
            tsv = f"{OUTDIR}/fusion/fusioncatcher/{{sample}}/final-list_candidate-fusion-genes.txt",
        log:   f"{LOGDIR}/fusioncatcher/{{sample}}.log"
        conda: "../../environments/11_fusioncatcher.yml"
        container: "file://containers/fusioncatcher.sif"
        threads: 16
        resources: mem_mb=64000, runtime=240
        params:
            data   = config["fusion"]["fusioncatcher_db"],
            outdir = f"{OUTDIR}/fusion/fusioncatcher/{{sample}}",
        shell:
            """
            fusioncatcher \
                -d {params.data} \
                -i {input.r1},{input.r2} \
                -o {params.outdir} \
                -p {threads} \
                2>{log}
            """
