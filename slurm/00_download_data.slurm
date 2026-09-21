#!/bin/bash
# --------------------------------------------------------------
### PART 1: Request resources
# --------------------------------------------------------------

#SBATCH --job-name=ngs101_download_all
#SBATCH --output=/xdisk/haining/maarowosegbe/RNA-seq/logs/slurm/%x-%j.out
#SBATCH --error=/xdisk/haining/maarowosegbe/RNA-seq/logs/slurm/%x-%j.err
#SBATCH --account=haining
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=48:00:00        # SRA downloads can be slow; 48 h for all datasets

# --------------------------------------------------------------
### PART 2: Download all NGS101 tutorial datasets
# --------------------------------------------------------------

echo "========== Download job started at: $(date) =========="
echo "  Node: $(hostname)  |  Job: ${SLURM_JOB_ID}"

SCRIPT_BASE=/home/u11/maarowosegbe/RNA-seq
DATA_BASE=/xdisk/haining/maarowosegbe/RNA-seq

export SCRIPT_BASE DATA_BASE THREADS=${SLURM_NTASKS}
mkdir -p "${DATA_BASE}/logs/slurm" "${DATA_BASE}/data/fastq" "${DATA_BASE}/tmp"

export MAMBA_ROOT_PREFIX=/groups/haining/maarowosegbe/micromamba
export MAMBA_EXE=/opt/ohpc/pub/apps/micromamba/2.0.2-2/bin/micromamba
MC="${MAMBA_EXE}"

# ── Which datasets to download? ───────────────────────────────────────
# Uncomment the modules you want. "all" downloads everything.
MODULE="${1:-core}"    # core | mirna | circrna | fusion | wgcna | crispr | all

echo ""
echo "--- Downloading module: ${MODULE} ---"
"${MC}" run -n rnaseq_env \
  bash "${SCRIPT_BASE}/scripts/download_all_data.sh" "${MODULE}"

echo ""
echo "========== Download finished at: $(date) =========="
echo ""
echo "Files in ${DATA_BASE}/data/:"
find "${DATA_BASE}/data" -name "*.fastq.gz" | sort | sed 's/^/  /'
