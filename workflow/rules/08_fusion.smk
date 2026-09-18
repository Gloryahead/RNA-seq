"""
Rules: Fusion gene detection (Parts 16–17)
  - STAR-Fusion (Part 16): runs via Apptainer SIF (fusiongene_env is a stub)
  - FusionCatcher (Part 17): fusioncatcher_env (pip install manually)

STAR-Fusion SIF: apptainer pull docker://trinityctat/starfusion:latest
CTAT library:    bash scripts/00b_download_starfusion_refs.sh
"""

OUTDIR = config["outdir"]
LOGDIR = config["logdir"]

# Apptainer SIF path (set once here; used in the STAR-Fusion shell block)
_STARFUSION_SIF = "/xdisk/haining/maarowosegbe/apptainer_images/starfusion_latest.sif"
_BIND_PATH      = "/xdisk/haining/maarowosegbe:/xdisk/haining/maarowosegbe"


if config["fusion"]["run_starfusion"]:

    rule star_fusion:
        """
        STAR-Fusion from chimeric junctions produced by STAR.
        Runs via Apptainer (arriba/STAR-Fusion cannot be installed via conda on UA HPC).
        Requires CTAT genome library; set fusion.ctat_lib in config.yaml.
        Download: bash scripts/00b_download_starfusion_refs.sh
        """
        input:
            chimeric = f"{OUTDIR}/bam/{{sample}}_Chimeric.out.junction",
        output:
            tsv = f"{OUTDIR}/fusion/starfusion/{{sample}}/star-fusion.fusion_predictions.tsv",
        log:   f"{LOGDIR}/starfusion/{{sample}}.log"
        threads: 8
        resources: mem_mb=64000, runtime=120, slurm_partition="standard"
        params:
            ctat_lib = config["fusion"]["ctat_lib"],
            outdir   = f"{OUTDIR}/fusion/starfusion/{{sample}}",
            sif      = _STARFUSION_SIF,
            bind     = _BIND_PATH,
        shell:
            """
            set -eo pipefail
            mkdir -p {params.outdir}
            apptainer exec \
                --bind {params.bind} \
                {params.sif} \
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
        Install manually: micromamba activate fusioncatcher_env && pip install git+...
        Set fusion.fusioncatcher_db in config.yaml.
        """
        input:
            r1 = get_r1,
            r2 = get_r2,
        output:
            tsv = f"{OUTDIR}/fusion/fusioncatcher/{{sample}}/final-list_candidate-fusion-genes.txt",
        log:   f"{LOGDIR}/fusioncatcher/{{sample}}.log"
        threads: 16
        resources: mem_mb=64000, runtime=240, slurm_partition="standard"
        params:
            activate = mamba_activate("fusioncatcher_env"),
            data     = config["fusion"]["fusioncatcher_db"],
            outdir   = f"{OUTDIR}/fusion/fusioncatcher/{{sample}}",
        shell:
            """
            set -eo pipefail
            {params.activate}
            fusioncatcher \
                -d {params.data} \
                -i {input.r1},{input.r2} \
                -o {params.outdir} \
                -p {threads} \
                2>{log}
            """
