#!/bin/bash
# ──────────────────────────────────────────────────────────────────────
# SPLAdder alternative splicing — Chunlong Ma RNA-seq (9040-0001)
# 3 human samples, GRCh38/Ensembl 110, paired-end 151 bp
#
# Samples:
#   S9040-0001-0001 = Control  (C3-4sU)
#   S9040-0001-0002 = Disease  (014-4sU)
#   S9040-0001-0003 = Remedy   (003-4sU)
#
# Workflow:
#   Step 1: spladder build — augment splice graph from all BAMs
#   Step 2: spladder test  — differential splicing per comparison
#
# NOTE: n=1 per group — p-values are not reliable without replicates.
# Use effect size / delta-PSI as the main metric.
# ──────────────────────────────────────────────────────────────────────

#SBATCH --job-name=spladder_chunlong
#SBATCH --output=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/logs/spladder_%j.out
#SBATCH --error=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/logs/spladder_%j.err
#SBATCH --account=haining
#SBATCH --partition=standard
#SBATCH --ntasks=8
#SBATCH --time=06:00:00

# ── Paths ──────────────────────────────────────────────────────────────
BAM_DIR=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/9040-0001-human-analysis/02.BamFiles
GTF=/xdisk/haining/maarowosegbe/refs/ensembl/Homo_sapiens.GRCh38.110.gtf
OUTDIR=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/9040-0001-human-analysis/05.SPLAdder
LOGDIR=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/logs
THREADS=${SLURM_NTASKS}
READ_LENGTH=151

CONTROL="${BAM_DIR}/S9040-0001-0001_Aligned.sortedByCoord.out.bam"
DISEASE="${BAM_DIR}/S9040-0001-0002_Aligned.sortedByCoord.out.bam"
REMEDY="${BAM_DIR}/S9040-0001-0003_Aligned.sortedByCoord.out.bam"

# ── Environment ────────────────────────────────────────────────────────
export MAMBA_ROOT_PREFIX=/groups/haining/maarowosegbe/micromamba
eval "$(/opt/ohpc/pub/apps/micromamba/2.0.2-2/bin/micromamba shell hook --shell bash)"
micromamba activate spladder_env

mkdir -p "${OUTDIR}/graph" "${LOGDIR}"

echo "========== SPLAdder started at: $(date) =========="

# ── Step 1: build augmented splice graph from all BAMs ─────────────────
echo "--- Step 1: build splice graph ---"
spladder build \
    --bams "${CONTROL},${DISEASE},${REMEDY}" \
    --gtf  "${GTF}" \
    --outdir "${OUTDIR}/graph" \
    --threads ${THREADS} \
    --event-types exon_skip,intron_retention,alt_3prime,alt_5prime,mult_exon_skip \
    --readlen ${READ_LENGTH} \
    2>"${LOGDIR}/spladder_build.log"
echo "  Done: ${OUTDIR}/graph/"

# ── Step 2: test differential splicing per comparison ─────────────────
echo "--- Step 2: differential splicing tests ---"

run_test() {
    local LABEL=$1
    local COND_A=$2
    local COND_B=$3
    echo "  Testing ${LABEL}..."
    spladder test \
        --outdir "${OUTDIR}/graph" \
        --conditionA "${COND_A}" \
        --conditionB "${COND_B}" \
        --output "${OUTDIR}/${LABEL}_events.tsv" \
        2>"${LOGDIR}/spladder_test_${LABEL}.log"
    echo "  Done: ${OUTDIR}/${LABEL}_events.tsv"
}

run_test "Disease_vs_Control" "${DISEASE}" "${CONTROL}"
run_test "Remedy_vs_Control"  "${REMEDY}"  "${CONTROL}"
run_test "Remedy_vs_Disease"  "${REMEDY}"  "${DISEASE}"

echo "========== SPLAdder finished at: $(date) =========="
echo ""
echo "Outputs:"
echo "  Splice graph : ${OUTDIR}/graph/"
echo "  Disease_vs_Control : ${OUTDIR}/Disease_vs_Control_events.tsv"
echo "  Remedy_vs_Control  : ${OUTDIR}/Remedy_vs_Control_events.tsv"
echo "  Remedy_vs_Disease  : ${OUTDIR}/Remedy_vs_Disease_events.tsv"
