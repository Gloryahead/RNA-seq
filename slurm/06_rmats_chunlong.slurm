#!/bin/bash
# ──────────────────────────────────────────────────────────────────────
# rMATS differential alternative splicing — Chunlong Ma RNA-seq (9040-0001)
# 3 human samples, GRCh38/Ensembl 110, paired-end 151 bp
#
# Outputs (per comparison, per event type):
#   SE.MATS.JC.txt   — skipped exon
#   A3SS.MATS.JC.txt — alternative 3' splice site
#   A5SS.MATS.JC.txt — alternative 5' splice site
#   MXE.MATS.JC.txt  — mutually exclusive exon
#   RI.MATS.JC.txt   — retained intron
#
# ── BEFORE RUNNING: set B1_SAMPLES and B2_SAMPLES below ───────────────
# ──────────────────────────────────────────────────────────────────────

#SBATCH --job-name=rmats_chunlong
#SBATCH --output=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/logs/rmats_%j.out
#SBATCH --error=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/logs/rmats_%j.err
#SBATCH --account=haining
#SBATCH --partition=standard
#SBATCH --ntasks=8
#SBATCH --time=04:00:00

# ── Paths ──────────────────────────────────────────────────────────────
BAM_DIR=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/9040-0001-human-analysis/02.BamFiles
GTF=/xdisk/haining/maarowosegbe/refs/ensembl/Homo_sapiens.GRCh38.110.gtf
OUTDIR=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/9040-0001-human-analysis/04.rMATS
LOGDIR=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/logs
THREADS=${SLURM_NTASKS}
READ_LENGTH=151

# ── Sample grouping — EDIT THESE ──────────────────────────────────────
# List BAMs for each group, comma-separated (no spaces)
# Example: B1=treatment, B2=control
B1_BAMS="${BAM_DIR}/S9040-0001-XXXX_Aligned.sortedByCoord.out.bam,\
${BAM_DIR}/S9040-0001-XXXX_Aligned.sortedByCoord.out.bam"

B2_BAMS="${BAM_DIR}/S9040-0001-XXXX_Aligned.sortedByCoord.out.bam"
# ──────────────────────────────────────────────────────────────────────

TMPDIR="${OUTDIR}/tmp"

# ── Environment ────────────────────────────────────────────────────────
export MAMBA_ROOT_PREFIX=/groups/haining/maarowosegbe/micromamba
eval "$(/opt/ohpc/pub/apps/micromamba/2.0.2-2/bin/micromamba shell hook --shell bash)"
micromamba activate transcript_splicing_env

mkdir -p "${OUTDIR}" "${TMPDIR}" "${LOGDIR}"

echo "========== rMATS started at: $(date) =========="
echo "  Group 1 (b1): ${B1_BAMS}"
echo "  Group 2 (b2): ${B2_BAMS}"

rmats.py \
    --b1 "${B1_BAMS}" \
    --b2 "${B2_BAMS}" \
    --gtf "${GTF}" \
    -t paired \
    --readLength ${READ_LENGTH} \
    --libType fr-unstranded \
    --od "${OUTDIR}" \
    --tmp "${TMPDIR}" \
    --nthread ${THREADS} \
    2>"${LOGDIR}/rmats.log"

echo "========== rMATS finished at: $(date) =========="
echo ""
echo "Outputs in: ${OUTDIR}/"
echo "  SE.MATS.JC.txt   — skipped exon"
echo "  A3SS.MATS.JC.txt — alternative 3' splice site"
echo "  A5SS.MATS.JC.txt — alternative 5' splice site"
echo "  MXE.MATS.JC.txt  — mutually exclusive exon"
echo "  RI.MATS.JC.txt   — retained intron"
