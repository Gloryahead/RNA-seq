# NGS101 RNA-seq Pipeline

A fully reproducible, modular RNA-seq pipeline covering all 21 parts of the [NGS101 tutorial series](https://ngs101.com). Runs on an HPC with SLURM via Snakemake, with every tool pinned in conda environments and optionally containerized with Apptainer (Singularity).

---

## What this pipeline does

| Module | Tutorial | Tools |
|--------|----------|-------|
| QC + alignment | Parts 1–3 | FastQC · Trim Galore · STAR · featureCounts · Salmon · MultiQC |
| DEG analysis | Parts 4, 10 | DESeq2 · edgeR · limma-voom · ggplot2 · pheatmap |
| Pathway enrichment | Parts 4, 7 | clusterProfiler · KEGG · MSigDB GSEA |
| Cancer subtype / GSVA | Parts 5–6, 11 | PAM50 · GSVA · CIBERSORTx prep |
| Alternative splicing | Parts 8–9 | SplAdder · SUPPA2 |
| RNA editing | Part 10 | REDItools2 |
| circRNA | Parts 13, 13.2 | CIRCexplorer2 · CIRI3 |
| miRNA / UMI-miRNA | Parts 15, 15.2 | miRDeep2 · umi_tools |
| Fusion genes | Parts 16–17 | STAR-Fusion · FusionCatcher |
| Viral detection | Parts 18–19 | esViritu · VIRTUS2 |
| CRISPR screens | Part 20 | MAGeCK |
| WGCNA | WGCNA tutorial | WGCNA · STRINGdb |
| Gene networks | GENIE3 tutorial | GENIE3 · DoRothEA |
| Master regulators | MR tutorial | RegEnrich · RTN |

---

## Repository layout

```
.
├── config/
│   ├── config.yaml          # ALL pipeline parameters — edit this first
│   └── samples.tsv          # sample metadata (sample, group, fastq paths, batch)
├── containers/
│   ├── rnaseq_core.def      # Apptainer definition (core alignment env)
│   ├── rnaseq_r.def         # Apptainer definition (R/Bioconductor env)
│   └── build_all_sifs.sh    # Build all SIF containers at once
├── environments/
│   ├── 01_rnaseq_core.yml   # conda YAML: STAR, Salmon, featureCounts, …
│   ├── 02_rnaseq_r.yml      # conda YAML: DESeq2, edgeR, limma, ggplot2, …
│   ├── …                    # (16 total YAMLs, one per analysis module)
│   ├── r_scripts/           # R install scripts for packages not in bioconda
│   └── setup_all_envs.sh    # Create all conda envs in one command
├── scripts/                 # Standalone bash scripts (step-by-step manual run)
│   ├── 00_download_refs.sh  # Download genome FASTA + GTF + miRBase
│   ├── 01_build_indices.sh  # Build STAR, Salmon, Bowtie v1, BWA indices
│   ├── 02_download_sra.sh   # Download FASTQs from NCBI SRA
│   ├── 03_qc_trim.sh        # FastQC + Trim Galore per sample
│   ├── 04_align_star.sh     # STAR alignment → sorted BAM + chimeric junctions
│   ├── 05_quantify.sh       # featureCounts (genes) + Salmon (transcripts)
│   ├── 06_multiqc.sh        # Aggregate all QC into one MultiQC report
│   ├── 07_convert_formats.sh # BAM→bigWig, BAM→BED, GTF→BED12, SAM→BAM
│   └── run_pipeline.sh      # Launch full Snakemake pipeline on SLURM
└── workflow/
    ├── Snakefile             # Master DAG
    ├── profiles/slurm/       # SLURM executor config
    ├── rules/               # 12 Snakemake rule files (one per module)
    └── scripts/             # R analysis scripts called by Snakemake rules
```

---

## Step-by-step setup

### Step 0 — Clone and configure

```bash
git clone https://github.com/<your-org>/ngs101-pipeline.git
cd ngs101-pipeline
```

Open `config/config.yaml` and set:
- `genome.star_index`, `genome.fasta`, `genome.gtf` — paths to your reference files
- `genome.salmon_index`
- `samples` — path to your `config/samples.tsv`
- `comparisons` — list of `[treatment, control]` pairs for DEG

Fill in `config/samples.tsv`:

```tsv
sample	group	fastq_r1	fastq_r2	batch
KO_rep1	KO	/data/KO_rep1_R1.fastq.gz	/data/KO_rep1_R2.fastq.gz	1
WT_rep1	WT	/data/WT_rep1_R1.fastq.gz	/data/WT_rep1_R2.fastq.gz	1
```

---

### Step 1 — Install Miniforge (if not already on HPC)

```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p $HOME/miniforge3
source $HOME/miniforge3/etc/profile.d/conda.sh
conda install -n base -c conda-forge mamba snakemake -y
```

---

### Step 2 — Build conda environments

This creates all 16 environments. Takes ~1–2 hours the first time.

```bash
bash environments/setup_all_envs.sh
```

To create a specific environment only:

```bash
mamba env create -f environments/01_rnaseq_core.yml
```

---

### Step 3 — Download reference genome

```bash
# Human (hg38 / GRCh38)
bash scripts/00_download_refs.sh

# Mouse (mm39)
bash scripts/00_download_refs.sh mouse

# Rat (mRatBN7.2)
bash scripts/00_download_refs.sh rat
```

Downloads to `ref/`: genome FASTA, GTF, transcript FASTA, miRBase sequences, refFlat.

---

### Step 4 — Build genome indices

Needs to run once. STAR index takes ~32 GB RAM and ~2 hours for hg38.

```bash
THREADS=16 bash scripts/01_build_indices.sh
```

Individual indices:

```bash
bash scripts/01_build_indices.sh --star    # STAR only (~32 GB RAM)
bash scripts/01_build_indices.sh --salmon  # Salmon only (<8 GB RAM)
bash scripts/01_build_indices.sh --bowtie  # Bowtie v1 for miRDeep2
bash scripts/01_build_indices.sh --bwa     # BWA for CIRI3 circRNA
```

---

### Step 5 — Download FASTQs (if using SRA data)

Create a text file with one SRA accession per line:

```bash
# data/sra_accessions.txt
SRR1234567
SRR7654321
```

Then run:

```bash
bash scripts/02_download_sra.sh data/sra_accessions.txt
```

If you already have FASTQ files, skip this step and put the paths in `config/samples.tsv` directly.

---

### Step 6 — Run the pipeline

**Option A: Full Snakemake pipeline on SLURM (recommended)**

```bash
# Dry run first — shows what will run without executing
bash scripts/run_pipeline.sh --dry-run

# Full run with conda environments
bash scripts/run_pipeline.sh

# Full run with Apptainer containers (build SIFs first — see Step 7)
bash scripts/run_pipeline.sh --singularity
```

**Option B: Run step by step (manually, useful for debugging)**

```bash
THREADS=16 bash scripts/03_qc_trim.sh       # FastQC + Trim Galore
THREADS=16 bash scripts/04_align_star.sh    # STAR alignment
THREADS=16 bash scripts/05_quantify.sh      # featureCounts + Salmon
bash scripts/06_multiqc.sh                  # MultiQC QC report
```

Then run individual R scripts:

```bash
# DEG analysis
conda run -n rnaseq_r_env Rscript workflow/scripts/deg_analysis.R \
  --counts  results/counts/counts_raw.tsv \
  --samples config/samples.tsv \
  --comparison KO_vs_WT \
  --method  DESeq2 \
  --outdir  results/deg/KO_vs_WT

# Pathway enrichment
conda run -n rnaseq_r_env Rscript workflow/scripts/pathway_analysis.R \
  --deg     results/deg/KO_vs_WT/KO_vs_WT_DEG_results.tsv \
  --organism human \
  --outdir  results/pathways/KO_vs_WT

# WGCNA
conda run -n wgcna_env Rscript workflow/scripts/wgcna.R \
  --counts  results/counts/counts_raw.tsv \
  --samples config/samples.tsv \
  --threads 8 \
  --outdir  results/networks/wgcna
```

---

### Step 7 — (Optional) Build Apptainer containers

Build SIF files **before** running the pipeline in singularity mode.  
Requires `apptainer >= 1.2` (usually `module load apptainer` on HPC).

```bash
# Build all containers (~20–30 GB total, takes 1–2 hours)
bash containers/build_all_sifs.sh

# Or build one at a time
apptainer build containers/rnaseq_core.sif containers/rnaseq_core.def
```

Then launch with:

```bash
bash scripts/run_pipeline.sh --singularity
```

---

### Step 8 — Format conversions (as needed)

```bash
# BAM → bigWig coverage track (for IGV / UCSC browser)
bash scripts/07_convert_formats.sh bam2bw

# BAM → BED (for BEDTools operations)
bash scripts/07_convert_formats.sh bam2bed

# GTF → BED12 (for RSeQC / infer_experiment.py)
bash scripts/07_convert_formats.sh gtf2bed

# BAM → FASTQ (recover reads)
bash scripts/07_convert_formats.sh bam2fastq
```

---

## Known limitations and manual steps

| Tool | Issue | Workaround |
|------|-------|------------|
| REDItools2 | Requires Python 2.7 (EOL) | Isolated in `05_reditools2.yml`; downstream analysis in `rnaseq_r_env` |
| STAR-Fusion | Docker officially preferred | conda env works on HPC; use Apptainer SIF for production |
| CIBERSORTx | Web-only tool | `12_cancer_subtype.smk` writes the mixture matrix for upload at cibersortx.stanford.edu |
| miRDeep2 | Requires Bowtie v1, not v2 | Isolated in `09_mirna.yml`; `01_build_indices.sh --bowtie` builds the correct index |
| CIRI3 | Java JAR not in conda | Download from GitHub releases; set `circrna.ciri3_jar` in config.yaml |
| STAR-Fusion CTAT library | ~30 GB download | Download separately: `wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/` |
| FusionCatcher DB | ~15 GB per organism | Run `fusioncatcher-build -g homo_sapiens` after creating the fusioncatcher_env |
| Part 21 (batch correction) | Tutorial returned 404 | sva ComBat + limma::removeBatchEffect are already in `rnaseq_r_env`; see `normalization.R` |

---

## Configuration quick reference

All parameters live in `config/config.yaml`. Key flags for toggling modules:

```yaml
# Turn modules on/off
splicing:
  run_spladder: true
  run_suppa2: true
circrna:
  run_circexplorer2: true
  run_ciri3: false
mirna:
  enabled: true
fusion:
  run_starfusion: true
  run_fusioncatcher: false
viral:
  run_esviritu: true
crispr:
  enabled: false
networks:
  run_wgcna: true
  run_genie3: false
```

---

## Reproducibility

- All conda environments pin exact package versions
- `config/config.yaml` records every parameter used
- Snakemake logs every shell command with timestamps in `logs/`
- `--use-singularity` mode uses immutable SIF containers for full bit-for-bit reproducibility
- Snakemake `--report` generates an HTML provenance report after each run
