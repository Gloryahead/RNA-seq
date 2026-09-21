#!/bin/bash
# --------------------------------------------------------------
### PART 1: Request resources
# --------------------------------------------------------------
# This script runs the pipeline step-by-step (without Snakemake),
# useful for debugging, learning, or running a single module.
# Mirror of your eclip_pipeline pattern.

#SBATCH --job-name=ngs101_steps
#SBATCH --output=/xdisk/haining/maarowosegbe/RNA-seq/logs/slurm/%x-%j.out
#SBATCH --error=/xdisk/haining/maarowosegbe/RNA-seq/logs/slurm/%x-%j.err
#SBATCH --account=haining
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem-per-cpu=8gb      # 32 × 8gb = 256 GB — enough for STAR + DESeq2
#SBATCH --time=240:00:00

# --------------------------------------------------------------
### PART 2: Setup and run steps
# --------------------------------------------------------------

echo "========== Job started at: $(date) =========="
echo "  Node    : $(hostname)"
echo "  Job ID  : ${SLURM_JOB_ID}"

# Two-path pattern — code in home, data in xdisk
SCRIPT_BASE=/home/u11/maarowosegbe/RNA-seq
DATA_BASE=/xdisk/haining/maarowosegbe/RNA-seq

export THREADS=${SLURM_NTASKS}
export DATA_BASE
# Read organism from config so it propagates to R scripts
ORGANISM=$(grep '^organism:' "${SCRIPT_BASE}/config/config.yaml" | awk '{print $2}' | tr -d '"')

mkdir -p "${DATA_BASE}/logs/slurm"

export MAMBA_ROOT_PREFIX=/groups/haining/maarowosegbe/micromamba
export MAMBA_EXE=/opt/ohpc/pub/apps/micromamba/2.0.2-2/bin/micromamba
MC="${MAMBA_EXE}"

# ── Step 1: QC + trimming ─────────────────────────────────────────────
"${MC}" run -n rnaseq_env bash "${SCRIPT_BASE}/scripts/03_qc_trim.sh"

# ── Step 2: STAR alignment ────────────────────────────────────────────
"${MC}" run -n rnaseq_env bash "${SCRIPT_BASE}/scripts/04_align_star.sh"

# ── Step 3: Gene + transcript quantification ──────────────────────────
"${MC}" run -n rnaseq_env bash "${SCRIPT_BASE}/scripts/05_quantify.sh"

# ── Step 4: MultiQC ───────────────────────────────────────────────────
"${MC}" run -n rnaseq_env bash "${SCRIPT_BASE}/scripts/06_multiqc.sh"

# ── Step 5: DEG analysis (R) ─────────────────────────────────────────
"${MC}" run -n rnaseq_r_env \
  Rscript "${SCRIPT_BASE}/workflow/scripts/deg_analysis.R" \
    --counts   "${DATA_BASE}/results/counts/counts_raw.tsv" \
    --samples  "${SCRIPT_BASE}/config/samples.tsv" \
    --outdir   "${DATA_BASE}/results/deg" \
    --method   DESeq2 \
    2>"${DATA_BASE}/logs/slurm/deg_${SLURM_JOB_ID}.log"

# ── Step 6: Pathway enrichment (R) ───────────────────────────────────
"${MC}" run -n rnaseq_r_env \
  Rscript "${SCRIPT_BASE}/workflow/scripts/pathway_analysis.R" \
    --deg      "${DATA_BASE}/results/deg/$(ls ${DATA_BASE}/results/deg/*_DEG_results.tsv | head -1)" \
    --organism "${ORGANISM}" \
    --outdir   "${DATA_BASE}/results/pathways" \
    2>"${DATA_BASE}/logs/slurm/pathways_${SLURM_JOB_ID}.log"

# ── Optional: WGCNA (comment out if not needed) ───────────────────────
# "${MC}" run -n wgcna_env \
#   Rscript "${SCRIPT_BASE}/workflow/scripts/wgcna.R" \
#     --counts  "${DATA_BASE}/results/counts/counts_raw.tsv" \
#     --samples "${SCRIPT_BASE}/config/samples.tsv" \
#     --threads "${THREADS}" \
#     --outdir  "${DATA_BASE}/results/networks/wgcna"

echo ""
echo "========== Job finished at: $(date) =========="
