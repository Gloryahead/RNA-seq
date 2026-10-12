#!/bin/bash
# Recovery: generate transcriptome.fa and build Salmon index.
# Run this because job 24213014 completed STAR but gffread was not on PATH,
# so transcriptome.fa was never created and the Salmon index was skipped.
#
# Submit: sbatch slurm/13_chunlong_salmon_index.slurm

#SBATCH --job-name=chunlong_salmon_idx
#SBATCH --output=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/logs/slurm/%x-%j.out
#SBATCH --error=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/logs/slurm/%x-%j.err
#SBATCH --account=haining
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=5120    # 16 × 5 GB = 80 GB; Salmon needs ~8 GB
#SBATCH --time=01:00:00

set -euo pipefail

echo "========== Salmon index recovery started at: $(date) =========="
echo "  Node    : $(hostname)"
echo "  Job ID  : ${SLURM_JOB_ID}"

export MAMBA_ROOT_PREFIX=/groups/haining/maarowosegbe/micromamba
export MAMBA_EXE=/opt/ohpc/pub/apps/micromamba/2.0.2-2/bin/micromamba
MC="${MAMBA_EXE}"

REF_DIR=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/ref
THREADS=16

# ── Step 1: generate transcriptome.fa (gffread via rnaseq_env) ────────
echo ""
echo "--- Extracting transcript sequences ---"
"${MC}" run -n rnaseq_env \
  gffread "${REF_DIR}/annotation.gtf" \
    -g "${REF_DIR}/genome.fa" \
    -w "${REF_DIR}/transcriptome.fa"
echo "  -> ${REF_DIR}/transcriptome.fa"

# ── Step 2: build Salmon index ────────────────────────────────────────
echo ""
echo "--- Building Salmon index ---"
SALMON_DECOYS="${REF_DIR}/salmon_decoys.txt"
grep "^>" "${REF_DIR}/genome.fa" | cut -d ' ' -f 1 | tr -d '>' > "${SALMON_DECOYS}"

"${MC}" run -n rnaseq_env \
  salmon index \
    -t "${REF_DIR}/transcriptome.fa" \
    -d "${SALMON_DECOYS}" \
    -i "${REF_DIR}/salmon_index" \
    -p "${THREADS}"
echo "  -> ${REF_DIR}/salmon_index"

echo ""
echo "========== Salmon index finished at: $(date) =========="
echo ""
echo "Next step: sbatch slurm/11_chunlong_rnaseq.slurm"
