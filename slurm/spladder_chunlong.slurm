#!/bin/bash
# ──────────────────────────────────────────────────────────────────────
# SPLAdder alternative splicing — Chunlong Ma RNA-seq (9040-0001)
# 3 human samples, GRCh38/Ensembl 110, paired-end 151 bp
#
#   S9040-0001-0001 = Control  (C3-4sU)
#   S9040-0001-0002 = Disease  (014-4sU)
#   S9040-0001-0003 = Remedy   (003-4sU)
#
# NOTE: n=1 per group — p-values are not reliable without replicates.
# Use effect size / delta-PSI as the main metric.
#
# RUN HISTORY:
#   24013309 — OOM: --ntasks=32 spread over 10 nodes; SPLAdder got 4 CPUs/20 GB
#   24021154 — build OK (peak 32 GB); Disease_vs_Control OK; then crashed
#              because spladder test defaults to mutex_exons, which build
#              was not asked to make -> test now gets the same --event-types
#   This run — build skipped (done); remaining two comparisons only
# ──────────────────────────────────────────────────────────────────────

#SBATCH --job-name=spladder_chunlong
#SBATCH --output=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/logs/spladder_%j.out
#SBATCH --error=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/logs/spladder_%j.err
#SBATCH --account=haining
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16           # 80 GB; measured peak was 32 GB
#SBATCH --time=24:00:00

# ── Paths ──────────────────────────────────────────────────────────────
BAM_DIR=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/9040-0001-human-analysis/02.BamFiles
GTF=/xdisk/haining/maarowosegbe/refs/ensembl/Homo_sapiens.GRCh38.110.gtf
OUTDIR=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/9040-0001-human-analysis/05.SPLAdder
LOGDIR=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/logs
THREADS=8
READ_LENGTH=151
EVENTS=exon_skip,intron_retention,alt_3prime,alt_5prime,mult_exon_skip   # build AND test must match

CONTROL="${BAM_DIR}/S9040-0001-0001_Aligned.sortedByCoord.out.bam"
DISEASE="${BAM_DIR}/S9040-0001-0002_Aligned.sortedByCoord.out.bam"
REMEDY="${BAM_DIR}/S9040-0001-0003_Aligned.sortedByCoord.out.bam"

# ── Environment ────────────────────────────────────────────────────────
export MAMBA_ROOT_PREFIX=/groups/haining/maarowosegbe/micromamba
eval "$(/opt/ohpc/pub/apps/micromamba/2.0.2-2/bin/micromamba shell hook --shell bash)"
micromamba activate spladder_env
set -eo pipefail

mkdir -p "${OUTDIR}/graph" "${LOGDIR}"

echo "========== SPLAdder started at: $(date) =========="
echo "Node: $(hostname) | CPUs: ${SLURM_CPUS_PER_TASK} | workers: ${THREADS}"

# ── Step 1: build — COMPLETED in job 24021154, kept for the record ─────
# spladder build \
#     -b "${CONTROL},${DISEASE},${REMEDY}" \
#     -a "${GTF}" \
#     -o "${OUTDIR}/graph" \
#     --readlen ${READ_LENGTH} \
#     --parallel ${THREADS} \
#     --event-types ${EVENTS} \
#     2>"${LOGDIR}/spladder_build.log"

# ── Step 2: differential splicing tests ────────────────────────────────
echo "--- Step 2: differential splicing tests ---"

run_test() {
    local LABEL=$1
    local COND_A=$2
    local COND_B=$3
    local LA="${LABEL%%_vs_*}"
    local LB="${LABEL##*_vs_}"
    echo "  Testing ${LABEL} at $(date)"
    spladder test \
        -o "${OUTDIR}/graph" \
        --conditionA "${COND_A}" \
        --conditionB "${COND_B}" \
        --labelA "${LA}" \
        --labelB "${LB}" \
        --event-types ${EVENTS} \
        2>"${LOGDIR}/spladder_test_${LABEL}.log"
    echo "  Done: ${OUTDIR}/graph/testing_${LA}_vs_${LB}/"
}

# run_test "Disease_vs_Control" "${DISEASE}" "${CONTROL}"   # completed in 24021154
run_test "Remedy_vs_Control"  "${REMEDY}"  "${CONTROL}"
run_test "Remedy_vs_Disease"  "${REMEDY}"  "${DISEASE}"

echo "========== SPLAdder finished at: $(date) =========="
