#!/bin/bash
# ──────────────────────────────────────────────────────────────────────
# SPLAdder visualization — MAPT gene (ENSG00000186868)
# Chunlong Ma RNA-seq (9040-0001), GRCh38/Ensembl 110
#
# Requires: 07_spladder_chunlong.slurm to have completed successfully.
# Runtime: ~5–15 minutes for a single gene.
#
# Output:
#   MAPT_splicing.pdf — coverage, events, and splice graph for MAPT
#   saved in VIZDIR below
# ──────────────────────────────────────────────────────────────────────

#SBATCH --job-name=spladder_viz
#SBATCH --output=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/logs/spladder_viz_%j.out
#SBATCH --error=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/logs/spladder_viz_%j.err
#SBATCH --account=haining
#SBATCH --partition=standard
#SBATCH --ntasks=4
#SBATCH --time=00:30:00

# ── Paths ──────────────────────────────────────────────────────────────
BAM_DIR=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/9040-0001-human-analysis/02.BamFiles
GRAPH_DIR=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/9040-0001-human-analysis/05.SPLAdder/graph
VIZDIR=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/9040-0001-human-analysis/05.SPLAdder/viz

CONTROL="${BAM_DIR}/S9040-0001-0001_Aligned.sortedByCoord.out.bam"
DISEASE="${BAM_DIR}/S9040-0001-0002_Aligned.sortedByCoord.out.bam"
REMEDY="${BAM_DIR}/S9040-0001-0003_Aligned.sortedByCoord.out.bam"

MAPT_GENE_ID="ENSG00000186868"

# ── Environment ────────────────────────────────────────────────────────
export MAMBA_ROOT_PREFIX=/groups/haining/maarowosegbe/micromamba
eval "$(/opt/ohpc/pub/apps/micromamba/2.0.2-2/bin/micromamba shell hook --shell bash)"
micromamba activate spladder_env

mkdir -p "${VIZDIR}"

# ── Guard: confirm build finished ─────────────────────────────────────
if [ ! -d "${GRAPH_DIR}" ]; then
    echo "ERROR: Graph directory not found: ${GRAPH_DIR}"
    echo "Run 07_spladder_chunlong.slurm first."
    exit 1
fi

echo "========== SPLAdder viz started at: $(date) =========="

# ── MAPT visualization ─────────────────────────────────────────────────
echo "--- Generating MAPT visualization ---"
spladder viz \
    --range gene "${MAPT_GENE_ID}" \
    --track coverage,segments \
    Control:"${CONTROL}" \
    Disease:"${DISEASE}" \
    Remedy:"${REMEDY}" \
    --track event any \
    --track splicegraph \
    -o "${GRAPH_DIR}" \
    -O "${VIZDIR}/MAPT_splicing" \
    --format pdf

echo "========== SPLAdder viz finished at: $(date) =========="
echo ""
echo "Output: ${VIZDIR}/MAPT_splicing.pdf"
