"""
Rules: Viral RNA detection (Parts 18–19)
  - esViritu (Part 18): esviritu_env (pip install manually)
  - VIRTUS2  (Part 19): virtus2_env + cwltool
"""

OUTDIR = config["outdir"]
LOGDIR = config["logdir"]


if config["viral"]["run_esviritu"]:

    rule esviritu_detect:
        """
        esViritu: virus-aware alignment against human + viral reference.
        Install manually: micromamba activate esviritu_env
                          pip install git+https://github.com/mtisza1/EsViritu.git
        Database: esviritu download_db -o ~/esviritu_db
        """
        input:
            r1 = get_r1,
            r2 = get_r2,
        output:
            tsv = f"{OUTDIR}/viral/esviritu/{{sample}}/virus_abundance.tsv",
        log:   f"{LOGDIR}/esviritu/{{sample}}.log"
        threads: 8
        resources: mem_mb=32000, runtime=120
        params:
            activate = mamba_activate("esviritu_env"),
            db       = config["viral"]["esviritu_db"],
            outdir   = f"{OUTDIR}/viral/esviritu/{{sample}}",
        shell:
            """
            set -eo pipefail
            {params.activate}
            esviritu detect \
                -1 {input.r1} \
                -2 {input.r2} \
                --db  {params.db} \
                -o    {params.outdir} \
                -t    {threads} \
                2>{log}
            """

    rule esviritu_merge:
        input:
            tsvs = expand(
                f"{OUTDIR}/viral/esviritu/{{sample}}/virus_abundance.tsv",
                sample=SAMPLES,
            ),
        output:
            merged = f"{OUTDIR}/viral/esviritu/all_samples_virus_abundance.tsv",
        log: f"{LOGDIR}/esviritu/merge.log"
        run:
            import pandas as pd
            frames = []
            for tsv, s in zip(input.tsvs, SAMPLES):
                df = pd.read_csv(tsv, sep="\t")
                df.insert(0, "sample", s)
                frames.append(df)
            pd.concat(frames).to_csv(output.merged, sep="\t", index=False)


if config["viral"]["run_virtus2"]:

    rule virtus2:
        """
        VIRTUS2 CWL workflow for viral transcript detection.
        Requires cwltool (in virtus2_env) and Singularity/Apptainer accessible to cwltool.
        Clone yyoshiaki/VIRTUS2 and set viral.virtus2_dir in config.
        """
        input:
            r1 = get_r1,
            r2 = get_r2,
        output:
            tsv = f"{OUTDIR}/viral/virtus2/{{sample}}/output.tsv",
        log:   f"{LOGDIR}/virtus2/{{sample}}.log"
        threads: 8
        resources: mem_mb=32000, runtime=180
        params:
            activate   = mamba_activate("virtus2_env"),
            cwl        = config["viral"]["virtus2_dir"],
            outdir     = f"{OUTDIR}/viral/virtus2/{{sample}}",
            star_human = config["genome"]["star_index"],
        shell:
            """
            set -eo pipefail
            {params.activate}
            mkdir -p {params.outdir}
            cwltool \
                {params.cwl}/workflow/VIRTUS.PE.cwl \
                --fastq  {input.r1} \
                --fastq2 {input.r2} \
                --genomeDir {params.star_human} \
                --outFileNamePrefix {params.outdir}/ \
                --nthreads {threads} \
                2>{log}
            """
