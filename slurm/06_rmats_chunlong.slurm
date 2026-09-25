#!/bin/bash
# ──────────────────────────────────────────────────────────────────────
# rMATS differential alternative splicing — Chunlong Ma RNA-seq (9040-0001)
# 3 human samples, GRCh38/Ensembl 110, paired-end 151 bp
#
# Samples:
#   S9040-0001-0001 = Control  (C3-4sU)
#   S9040-0001-0002 = Disease  (014-4sU)
#   S9040-0001-0003 = Remedy   (003-4sU)
#
# Comparisons (n=1 per group — results are exploratory, no replicates):
#   1. Disease_vs_Control  (0002 vs 0001)
#   2. Remedy_vs_Control   (0003 vs 0001)
#   3. Remedy_vs_Disease   (0003 vs 0002)
#
# Outputs per comparison, per event type:
#   SE.MATS.JC.txt / A3SS / A5SS / MXE / RI
# ──────────────────────────────────────────────────────────────────────

#SBATCH --job-name=rmats_chunlong
#SBATCH --output=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/logs/rmats_%j.out
#SBATCH --error=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/logs/rmats_%j.err
#SBATCH --account=haining
#SBATCH --partition=standard
#SBATCH --ntasks=16
#SBATCH --time=06:00:00

# ── Paths ──────────────────────────────────────────────────────────────
BAM_DIR=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/9040-0001-human-analysis/02.BamFiles
GTF=/xdisk/haining/maarowosegbe/refs/ensembl/Homo_sapiens.GRCh38.110.gtf
OUTBASE=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/9040-0001-human-analysis/04.rMATS
LOGDIR=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/logs
THREADS=${SLURM_NTASKS}
READ_LENGTH=151

CONTROL="${BAM_DIR}/S9040-0001-0001_Aligned.sortedByCoord.out.bam"
DISEASE="${BAM_DIR}/S9040-0001-0002_Aligned.sortedByCoord.out.bam"
REMEDY="${BAM_DIR}/S9040-0001-0003_Aligned.sortedByCoord.out.bam"

# ── Environment ────────────────────────────────────────────────────────
export MAMBA_ROOT_PREFIX=/groups/haining/maarowosegbe/micromamba
eval "$(/opt/ohpc/pub/apps/micromamba/2.0.2-2/bin/micromamba shell hook --shell bash)"
micromamba activate transcript_splicing_env

mkdir -p "${LOGDIR}"

run_rmats() {
    local LABEL=$1
    local B1=$2
    local B2=$3
    local OUTDIR="${OUTBASE}/${LABEL}"
    local TMPDIR="${OUTBASE}/tmp_${LABEL}"
    mkdir -p "${OUTDIR}" "${TMPDIR}"
    echo "--- ${LABEL}: $(date) ---"
    rmats.py \
        --b1 "${B1}" \
        --b2 "${B2}" \
        --gtf "${GTF}" \
        -t paired \
        --readLength ${READ_LENGTH} \
        --libType fr-unstranded \
        --od "${OUTDIR}" \
        --tmp "${TMPDIR}" \
        --nthread ${THREADS} \
        2>"${LOGDIR}/rmats_${LABEL}.log"
    echo "  Done → ${OUTDIR}/"
}

echo "========== rMATS started at: $(date) =========="

run_rmats "Disease_vs_Control" "${DISEASE}" "${CONTROL}"
run_rmats "Remedy_vs_Control"  "${REMEDY}"  "${CONTROL}"
run_rmats "Remedy_vs_Disease"  "${REMEDY}"  "${DISEASE}"

echo "========== rMATS finished at: $(date) =========="
echo ""
echo "NOTE: n=1 per group — p-values are not reliable without replicates."
echo "Use inclusion level difference (IncLevelDifference) as the main metric."
echo ""
echo "Outputs:"
echo "  ${OUTBASE}/Disease_vs_Control/"
echo "  ${OUTBASE}/Remedy_vs_Control/"
echo "  ${OUTBASE}/Remedy_vs_Disease/"
