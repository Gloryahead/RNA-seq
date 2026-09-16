#!/usr/bin/env bash
# 05_quantify.sh
# Gene-level quantification (featureCounts) + transcript-level (Salmon).
# Both are needed: featureCounts for DEG, Salmon for SUPPA2 splicing.
#
# Usage:
#   bash scripts/05_quantify.sh
#   THREADS=16 bash scripts/05_quantify.sh
#
# Output:
#   results/counts/counts_raw.tsv      <- featureCounts (all samples, one file)
#   results/salmon/<sample>/quant.sf   <- Salmon per-sample transcript TPM

set -euo pipefail

SCRIPT_BASE="${SCRIPT_BASE:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
DATA_BASE="${DATA_BASE:-/xdisk/haining/maarowosegbe/RNA-seq}"

SAMPLES_TSV="${SCRIPT_BASE}/config/samples.tsv"
BAM_DIR="${DATA_BASE}/results/bam"
GTF="${DATA_BASE}/ref/annotation.gtf"
SALMON_INDEX="${DATA_BASE}/ref/salmon_index"
TRANSCRIPTS="${DATA_BASE}/ref/transcriptome.fa"
OUTDIR_COUNTS="${DATA_BASE}/results/counts"
OUTDIR_SALMON="${DATA_BASE}/results/salmon"
THREADS="${THREADS:-16}"

mkdir -p "${OUTDIR_COUNTS}" "${OUTDIR_SALMON}" "${DATA_BASE}/logs/quantify"

[[ -f "${GTF}" ]]             || { echo "ERROR: ${GTF} not found"; exit 1; }
[[ -d "${SALMON_INDEX}" ]]    || { echo "ERROR: Salmon index not found at ${SALMON_INDEX}. Run scripts/01_build_indices.sh"; exit 1; }

# ── Collect all BAMs ──────────────────────────────────────────────────
BAMS=()
while IFS=$'\t' read -r sample _rest; do
  [[ -z "${sample}" || "${sample}" =~ ^# || "${sample}" == "sample" ]] && continue
  bam="${BAM_DIR}/${sample}_Aligned.sortedByCoord.out.bam"
  [[ -f "${bam}" ]] || { echo "WARNING: BAM not found: ${bam}"; continue; }
  BAMS+=("${bam}")
done < "${SAMPLES_TSV}"

[[ ${#BAMS[@]} -eq 0 ]] && {
  echo "ERROR: No BAMs found in ${BAM_DIR}. Run scripts/04_align_star.sh first."
  exit 1
}

echo "Found ${#BAMS[@]} BAMs"

# ── 1. featureCounts (gene-level, all samples at once) ────────────────
COUNTS_FILE="${OUTDIR_COUNTS}/counts_raw.tsv"
if [[ ! -f "${COUNTS_FILE}" ]]; then
  echo ""
  echo "=== featureCounts: gene-level quantification ==="
  mamba-haining run -n rnaseq_env \
    featureCounts \
      -T "${THREADS}" \
      -a "${GTF}" \
      -o "${COUNTS_FILE}" \
      -p \
      --countReadPairs \
      -B \
      -C \
      -Q 10 \
      --extraAttributes gene_name,gene_biotype \
      "${BAMS[@]}" \
      2>"logs/quantify/featurecounts.log"
  echo "  -> ${COUNTS_FILE}"
  echo "  -> ${COUNTS_FILE}.summary"
else
  echo "[featureCounts] Already exists, skipping"
fi

# ── 2. Salmon per-sample (transcript-level TPM for SUPPA2/tximport) ──
echo ""
echo "=== Salmon: transcript-level quantification ==="

while IFS=$'\t' read -r sample group r1 r2 batch || [[ -n "$sample" ]]; do
  [[ -z "${sample}" || "${sample}" =~ ^# || "${sample}" == "sample" ]] && continue

  SALMON_OUT="${OUTDIR_SALMON}/${sample}"
  if [[ -f "${SALMON_OUT}/quant.sf" ]]; then
    echo "  [${sample}] Already quantified, skipping"
    continue
  fi

  echo "  Quantifying: ${sample}"
  mkdir -p "${SALMON_OUT}"

  # Determine read arguments
  if [[ -n "${r2:-}" ]]; then
    # Paired-end: look for trimmed FASTQs first
    r1_trim="results/trimmed/${sample}_val_1.fq.gz"
    r2_trim="results/trimmed/${sample}_val_2.fq.gz"
    [[ -f "${r1_trim}" ]] && r1="${r1_trim}"
    [[ -f "${r2_trim}" ]] && r2="${r2_trim}"
    read_args="-1 ${r1} -2 ${r2}"
  else
    r1_trim="results/trimmed/${sample}_trimmed.fq.gz"
    [[ -f "${r1_trim}" ]] && r1="${r1_trim}"
    read_args="-r ${r1}"
  fi

  mamba-haining run -n rnaseq_env \
    salmon quant \
      -i "${SALMON_INDEX}" \
      -l A \
      ${read_args} \
      --validateMappings \
      --gcBias \
      --numBootstraps 50 \
      --seqBias \
      -p "${THREADS}" \
      -o "${SALMON_OUT}" \
      2>"logs/quantify/salmon_${sample}.log"

  echo "  Done: ${SALMON_OUT}/quant.sf"
done < "${SAMPLES_TSV}"

echo ""
echo "=== Quantification complete ==="
echo "  featureCounts : ${COUNTS_FILE}"
echo "  Salmon TPM    : ${OUTDIR_SALMON}/<sample>/quant.sf"
echo ""
echo "Next step: bash scripts/06_multiqc.sh"
