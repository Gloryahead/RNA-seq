#!/usr/bin/env bash
# download_all_data.sh
# Downloads all NGS101 tutorial datasets from SRA and Zenodo.
# Each module's accessions live in data/accessions/<N>_<name>.txt.
#
# Usage:
#   bash scripts/download_all_data.sh            # all datasets
#   bash scripts/download_all_data.sh core       # Parts 1-5 (main RNA-seq)
#   bash scripts/download_all_data.sh mirna      # Part 15
#   bash scripts/download_all_data.sh circrna    # Parts 13/13.2
#   bash scripts/download_all_data.sh fusion     # Parts 16/16.2
#   bash scripts/download_all_data.sh wgcna      # WGCNA tutorial
#   bash scripts/download_all_data.sh crispr     # Part 20 (Zenodo)
#
# Outputs go to: ${DATA_BASE}/data/fastq/<accession>_1.fastq.gz

set -euo pipefail

SCRIPT_BASE="${SCRIPT_BASE:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
DATA_BASE="${DATA_BASE:-/xdisk/haining/maarowosegbe/ngs101-pipeline}"
FASTQ_DIR="${DATA_BASE}/data/fastq"
THREADS="${THREADS:-8}"
MODULE="${1:-all}"

mkdir -p "${FASTQ_DIR}" "${DATA_BASE}/tmp"
export TMPDIR="${DATA_BASE}/tmp"

# ── Helper: download one accession via prefetch + fasterq-dump ────────
download_sra() {
  local acc="$1"
  local outdir="$2"
  mkdir -p "${outdir}"

  # Check if already done
  if ls "${outdir}/${acc}_"*.fastq.gz &>/dev/null 2>&1; then
    echo "  [${acc}] Already downloaded — skipping"
    return
  fi

  echo "  [${acc}] Prefetching..."
  prefetch --output-directory "${DATA_BASE}/tmp" --max-size 100G "${acc}"

  echo "  [${acc}] Converting to FASTQ..."
  fasterq-dump \
    "${DATA_BASE}/tmp/${acc}/${acc}.sra" \
    --outdir "${outdir}" \
    --threads "${THREADS}" \
    --split-files \
    --skip-technical \
    --progress

  echo "  [${acc}] Compressing..."
  for fq in "${outdir}/${acc}"*.fastq; do
    [[ -f "$fq" ]] && pigz -p "${THREADS}" "$fq" && echo "    -> ${fq}.gz"
  done

  rm -rf "${DATA_BASE}/tmp/${acc}"
}

# ── Helper: download a whole accession file ───────────────────────────
download_accession_file() {
  local file="$1"
  local outdir="$2"
  echo ""
  echo "=== Downloading from: $(basename "${file}") ==="
  while IFS= read -r line || [[ -n "$line" ]]; do
    # Skip comments and empty lines
    [[ -z "$line" || "$line" =~ ^# ]] && continue
    local acc
    acc=$(echo "$line" | tr -d '[:space:]')
    download_sra "${acc}" "${outdir}"
  done < "${file}"
}

# ── CRISPR: Zenodo download (not SRA) ────────────────────────────────
download_crispr() {
  local crispr_dir="${DATA_BASE}/data/crispr"
  mkdir -p "${crispr_dir}"
  echo ""
  echo "=== CRISPR: Zenodo 5750854 ==="
  local base="https://zenodo.org/records/5750854/files"
  for fname in "brunello.tsv" "T0-Control.fastq.gz" "T8-APR-246.fastq.gz" "T8-Vehicle.fastq.gz"; do
    if [[ -f "${crispr_dir}/${fname}" ]]; then
      echo "  ${fname}: Already exists — skipping"
    else
      echo "  Downloading ${fname}..."
      wget -q --show-progress -O "${crispr_dir}/${fname}" "${base}/${fname}"
    fi
  done
  echo "  -> ${crispr_dir}/"
}

# ── WGCNA: fetch all SRR IDs from GEO accession (requires Entrez utils) ─
fetch_wgcna_accessions() {
  local acc_file="${SCRIPT_BASE}/data/accessions/05_wgcna.txt"
  local wgcna_dir="${FASTQ_DIR}/wgcna"
  echo ""
  echo "=== WGCNA: fetching SRR IDs for GSE261875 ==="
  if command -v esearch &>/dev/null; then
    esearch -db sra -query GSE261875 \
      | efetch -format runinfo \
      | cut -d',' -f1 \
      | grep '^SRR' \
      > /tmp/wgcna_srrs.txt
    echo "  Found $(wc -l < /tmp/wgcna_srrs.txt) SRR accessions"
    while IFS= read -r acc; do
      download_sra "${acc}" "${wgcna_dir}"
    done < /tmp/wgcna_srrs.txt
  else
    echo "  WARNING: esearch not found. Install entrez-direct:"
    echo "    conda install -n rnaseq_env -c bioconda entrez-direct"
    echo "  Then re-run: bash scripts/download_all_data.sh wgcna"
  fi
}

# ── Dispatch ─────────────────────────────────────────────────────────
case "${MODULE}" in
  core)
    download_accession_file \
      "${SCRIPT_BASE}/data/accessions/01_core_rnaseq.txt" \
      "${FASTQ_DIR}/core"
    ;;
  mirna)
    download_accession_file \
      "${SCRIPT_BASE}/data/accessions/02_mirna.txt" \
      "${FASTQ_DIR}/mirna"
    ;;
  circrna)
    download_accession_file \
      "${SCRIPT_BASE}/data/accessions/03_circrna.txt" \
      "${FASTQ_DIR}/circrna"
    ;;
  fusion)
    download_accession_file \
      "${SCRIPT_BASE}/data/accessions/04_fusion.txt" \
      "${FASTQ_DIR}/fusion"
    ;;
  wgcna)
    fetch_wgcna_accessions
    ;;
  crispr)
    download_crispr
    ;;
  all)
    echo "=== Downloading ALL NGS101 datasets ==="
    download_accession_file "${SCRIPT_BASE}/data/accessions/01_core_rnaseq.txt" "${FASTQ_DIR}/core"
    download_accession_file "${SCRIPT_BASE}/data/accessions/02_mirna.txt"       "${FASTQ_DIR}/mirna"
    download_accession_file "${SCRIPT_BASE}/data/accessions/03_circrna.txt"     "${FASTQ_DIR}/circrna"
    download_accession_file "${SCRIPT_BASE}/data/accessions/04_fusion.txt"      "${FASTQ_DIR}/fusion"
    download_crispr
    fetch_wgcna_accessions
    ;;
  *)
    echo "Unknown module: ${MODULE}"
    echo "Usage: bash $0 [core|mirna|circrna|fusion|wgcna|crispr|all]"
    exit 1
    ;;
esac

echo ""
echo "=== Download complete ==="
echo "  FASTQs in: ${FASTQ_DIR}/"
echo ""
echo "Next: bash scripts/00_download_refs.sh"
