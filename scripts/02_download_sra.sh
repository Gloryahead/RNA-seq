#!/usr/bin/env bash
# 02_download_sra.sh
# Download FASTQ files from NCBI SRA using fasterq-dump.
# Reads a list of SRA accessions from a text file (one per line).
#
# Usage:
#   bash scripts/02_download_sra.sh data/sra_accessions.txt
#   bash scripts/02_download_sra.sh data/sra_accessions.txt --outdir /scratch/fastq
#
# sra_accessions.txt format (one SRA run accession per line):
#   SRR1234567
#   SRR7654321
#   ...
#
# Output: <outdir>/<accession>_1.fastq.gz and _2.fastq.gz (paired)
#         or <outdir>/<accession>.fastq.gz (single-end)

set -euo pipefail

DATA_BASE="${DATA_BASE:-/xdisk/haining/maarowosegbe/RNA-seq}"
ACCESSION_FILE="${1:-}"
OUTDIR="${2:-${DATA_BASE}/data/fastq}"
THREADS="${THREADS:-8}"
TMP_DIR="${DATA_BASE}/tmp/sra_tmp"

if [[ -z "${ACCESSION_FILE}" ]]; then
  echo "Usage: bash scripts/02_download_sra.sh <accession_list.txt> [outdir]"
  echo "       accession_list.txt: one SRA run accession per line (e.g. SRR1234567)"
  exit 1
fi

[[ -f "${ACCESSION_FILE}" ]] || { echo "ERROR: ${ACCESSION_FILE} not found"; exit 1; }

mkdir -p "${OUTDIR}" "${TMP_DIR}"

# Activate env with sra-tools
mamba-haining run -n rnaseq_env prefetch --version &>/dev/null \
  || { echo "ERROR: sra-tools not found. Activate rnaseq_env."; exit 1; }

echo "=== Downloading SRA accessions from ${ACCESSION_FILE} ==="

while IFS= read -r accession || [[ -n "$accession" ]]; do
  # Skip empty lines and comments
  [[ -z "${accession}" || "${accession}" =~ ^# ]] && continue
  accession=$(echo "${accession}" | tr -d '[:space:]')

  echo ""
  echo "--- ${accession} ---"

  # Check if already downloaded
  if ls "${OUTDIR}/${accession}"*.fastq.gz &>/dev/null; then
    echo "  Already exists, skipping"
    continue
  fi

  # Step 1: prefetch (downloads .sra file with resume support)
  echo "  [1/2] Prefetching ${accession}..."
  mamba-haining run -n rnaseq_env \
    prefetch \
      --output-directory "${TMP_DIR}" \
      --max-size 50G \
      "${accession}"

  # Step 2: fasterq-dump (converts .sra → .fastq, gzip in-place)
  echo "  [2/2] Converting to FASTQ..."
  mamba-haining run -n rnaseq_env \
    fasterq-dump \
      "${TMP_DIR}/${accession}/${accession}.sra" \
      --outdir "${OUTDIR}" \
      --threads "${THREADS}" \
      --split-files \
      --skip-technical \
      --progress

  # Compress with pigz (faster than gzip)
  echo "  Compressing..."
  for fq in "${OUTDIR}/${accession}"*.fastq; do
    [[ -f "$fq" ]] && pigz -p "${THREADS}" "$fq"
  done

  # Clean up prefetch cache
  rm -rf "${TMP_DIR}/${accession}"
  echo "  Done: $(ls "${OUTDIR}/${accession}"*.fastq.gz 2>/dev/null | xargs -I{} bash -c 'echo {} $(du -sh {} | cut -f1)')"

done < "${ACCESSION_FILE}"

echo ""
echo "=== SRA download complete ==="
echo "FASTQs in: ${OUTDIR}"
echo ""
echo "Next step: bash scripts/03_qc_trim.sh"
