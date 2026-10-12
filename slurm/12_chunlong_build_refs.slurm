#!/bin/bash
# --------------------------------------------------------------
# Download hg38 references and build STAR + Salmon indices
# for the Chunlong Ma RNA-seq project (9040-0001).
#
# Refs land in: /xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/ref/
# Run this ONCE before submitting 11_chunlong_rnaseq.slurm.
#
# Submit: sbatch slurm/12_chunlong_build_refs.slurm
# Wall time: ~5 h  (download ~30 min, STAR ~2 h, Salmon ~15 min, BWA skipped)
# --------------------------------------------------------------

#SBATCH --job-name=chunlong_build_refs
#SBATCH --output=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/logs/slurm/%x-%j.out
#SBATCH --error=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/logs/slurm/%x-%j.err
#SBATCH --account=haining
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=5120    # 32 × 5 GB = 160 GB; STAR genomeGenerate needs ~35 GB for hg38
#SBATCH --time=06:00:00

set -euo pipefail

echo "========== Build refs started at: $(date) =========="
echo "  Node    : $(hostname)"
echo "  Job ID  : ${SLURM_JOB_ID}"

SCRIPT_BASE=/home/u11/maarowosegbe/RNA-seq

# All ref files go under the Chunlong project directory, separate from mm10 refs
export DATA_BASE=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq
export THREADS=32

mkdir -p "${DATA_BASE}/logs/slurm" "${DATA_BASE}/ref"

# ── Step 1: download hg38 genome, GTF, and transcript FASTA ──────────
echo ""
echo "--- Downloading hg38 references (Ensembl 112) ---"
bash "${SCRIPT_BASE}/scripts/00_download_refs.sh" human "${THREADS}"

# ── Step 2: build STAR index (~2 h, ~35 GB RAM for hg38) ─────────────
echo ""
echo "--- Building STAR index ---"
bash "${SCRIPT_BASE}/scripts/01_build_indices.sh" --star

# ── Step 3: build Salmon index (~15 min) ─────────────────────────────
echo ""
echo "--- Building Salmon index ---"
bash "${SCRIPT_BASE}/scripts/01_build_indices.sh" --salmon

echo ""
echo "========== Build refs finished at: $(date) =========="
echo ""
echo "Ref directory: ${DATA_BASE}/ref"
echo "Next step: sbatch slurm/11_chunlong_rnaseq.slurm"
