#!/bin/bash
# ──────────────────────────────────────────────────────────────────────
# bamCoverage — generate CPM-normalized bigWig files for IGV
# Chunlong Ma RNA-seq (9040-0001): Control, Disease, Remedy
#
# Samples:
#   S9040-0001-0001 = Control  (C3-4sU)
#   S9040-0001-0002 = Disease  (014-4sU)
#   S9040-0001-0003 = Remedy   (003-4sU)
#
# Output: one .CPM.bw file per sample in 06.BigWig/
# Load these in IGV instead of BAMs — y-axis is CPM, depth-normalized.
# ──────────────────────────────────────────────────────────────────────

#SBATCH --job-name=bigwig_chunlong
#SBATCH --output=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/logs/bigwig_%j.out
#SBATCH --error=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/logs/bigwig_%j.err
#SBATCH --account=haining
#SBATCH --partition=standard
#SBATCH --ntasks=8
#SBATCH --time=02:00:00

# ── Paths ──────────────────────────────────────────────────────────────
BAM_DIR=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/9040-0001-human-analysis/02.BamFiles
BW_DIR=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/9040-0001-human-analysis/06.BigWig
THREADS=${SLURM_NTASKS}

# ── Environment ────────────────────────────────────────────────────────
export MAMBA_ROOT_PREFIX=/groups/haining/maarowosegbe/micromamba
eval "$(/opt/ohpc/pub/apps/micromamba/2.0.2-2/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq_env

mkdir -p "${BW_DIR}"

echo "========== bamCoverage started at: $(date) =========="

# ── Generate CPM-normalized bigWig per sample ──────────────────────────
for SAMPLE in S9040-0001-0001 S9040-0001-0002 S9040-0001-0003; do
    BAM="${BAM_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam"
    BW="${BW_DIR}/${SAMPLE}.CPM.bw"

    if [ ! -f "${BAM}" ]; then
        echo "WARNING: BAM not found, skipping: ${BAM}"
        continue
    fi

    echo "--- Processing ${SAMPLE} ---"
    bamCoverage \
        --bam "${BAM}" \
        --outFileName "${BW}" \
        --normalizeUsing CPM \
        --binSize 10 \
        --numberOfProcessors "${THREADS}" \
        --skipNonCoveredRegions

    echo "  Saved: ${BW}"
done

echo "========== bamCoverage finished at: $(date) =========="
echo ""
echo "Load these files in IGV (all on the same scale, depth-normalized):"
ls "${BW_DIR}"/*.CPM.bw 2>/dev/null
