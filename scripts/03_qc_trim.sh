#!/usr/bin/env bash
# 03_qc_trim.sh
# Per-sample QC (FastQC) + adapter trimming (Trim Galore).
# Reads config/samples.tsv; outputs to results/fastqc/ and results/trimmed/.
#
# Usage:
#   bash scripts/03_qc_trim.sh                    # all samples from samples.tsv
#   bash scripts/03_qc_trim.sh sample_id           # single sample (for testing)
#   THREADS=16 bash scripts/03_qc_trim.sh
#
# samples.tsv columns required: sample, fastq_r1, fastq_r2
#   (fastq_r2 can be empty for single-end)

set -euo pipefail

SAMPLES_TSV="config/samples.tsv"
OUTDIR_QC="results/fastqc"
OUTDIR_TRIM="results/trimmed"
THREADS="${THREADS:-8}"
SINGLE_SAMPLE="${1:-}"

mkdir -p "${OUTDIR_QC}" "${OUTDIR_TRIM}" logs/qc_trim

[[ -f "${SAMPLES_TSV}" ]] || { echo "ERROR: ${SAMPLES_TSV} not found"; exit 1; }

run_sample() {
  local sample="$1"
  local r1="$2"
  local r2="${3:-}"

  echo ""
  echo "═══ ${sample} ═══"

  # ── FastQC on raw reads ───────────────────────────────────────────
  echo "  [1/2] FastQC on raw reads..."
  if [[ -n "${r2}" ]]; then
    conda run -n rnaseq_env \
      fastqc -t "${THREADS}" -o "${OUTDIR_QC}" "${r1}" "${r2}" \
      2>"logs/qc_trim/${sample}_fastqc_raw.log"
  else
    conda run -n rnaseq_env \
      fastqc -t "${THREADS}" -o "${OUTDIR_QC}" "${r1}" \
      2>"logs/qc_trim/${sample}_fastqc_raw.log"
  fi

  # ── Trim Galore ───────────────────────────────────────────────────
  echo "  [2/2] Trim Galore..."
  local trim_args=(
    --cores "${THREADS}"
    --quality 20
    --length 20
    --fastqc
    --output_dir "${OUTDIR_TRIM}"
  )
  if [[ -n "${r2}" ]]; then
    trim_args+=(--paired)
    conda run -n rnaseq_env \
      trim_galore "${trim_args[@]}" "${r1}" "${r2}" \
      2>"logs/qc_trim/${sample}_trimgalore.log"
  else
    conda run -n rnaseq_env \
      trim_galore "${trim_args[@]}" "${r1}" \
      2>"logs/qc_trim/${sample}_trimgalore.log"
  fi
  echo "  Done: ${sample}"
}

# ── Read samples.tsv (skip header line) ──────────────────────────────
{
  read -r header   # skip header
  while IFS=$'\t' read -r sample group r1 r2 batch || [[ -n "$sample" ]]; do
    [[ -z "${sample}" || "${sample}" =~ ^# ]] && continue

    if [[ -n "${SINGLE_SAMPLE}" && "${sample}" != "${SINGLE_SAMPLE}" ]]; then
      continue
    fi

    [[ -f "${r1}" ]] || { echo "WARNING: R1 not found for ${sample}: ${r1}"; continue; }

    run_sample "${sample}" "${r1}" "${r2}"

  done
} < "${SAMPLES_TSV}"

echo ""
echo "=== QC + trimming complete ==="
echo "  FastQC reports : ${OUTDIR_QC}"
echo "  Trimmed FASTQs : ${OUTDIR_TRIM}"
echo ""
echo "Next step: bash scripts/04_align_star.sh"
