#!/usr/bin/env bash
# 00_download_refs.sh
# Download genome FASTA, GTF, and miRBase files for human (hg38).
# Change GENOME_BUILD and GTF_RELEASE for mouse/rat.
# Run this ONCE before building indices.
#
# Usage:
#   bash scripts/00_download_refs.sh
#   bash scripts/00_download_refs.sh --organism mouse
#
# Output directory: ref/  (ignored by git — add to HPC /ref or /scratch)

set -euo pipefail

ORGANISM="${1:-human}"     # human | mouse | rat
THREADS="${2:-${THREADS:-8}}"

# Code stays in home; references go to xdisk
SCRIPT_BASE="${SCRIPT_BASE:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
DATA_BASE="${DATA_BASE:-/xdisk/haining/maarowosegbe/RNA-seq}"
REF_DIR="${DATA_BASE}/ref"
mkdir -p "${REF_DIR}"

# ── Organism-specific settings ────────────────────────────────────────
case "${ORGANISM}" in
  human)
    ENSEMBL_RELEASE="112"
    ENSEMBL_SPECIES="homo_sapiens"
    UCSC_GENOME="hg38"
    FASTA_URL="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/${ENSEMBL_SPECIES}/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    GTF_URL="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gtf/${ENSEMBL_SPECIES}/Homo_sapiens.GRCh38.${ENSEMBL_RELEASE}.gtf.gz"
    MIRBASE_SPECIES="hsa"
    ;;
  mouse)
    ENSEMBL_RELEASE="112"
    ENSEMBL_SPECIES="mus_musculus"
    UCSC_GENOME="mm39"
    FASTA_URL="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/${ENSEMBL_SPECIES}/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
    GTF_URL="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gtf/${ENSEMBL_SPECIES}/Mus_musculus.GRCm39.${ENSEMBL_RELEASE}.gtf.gz"
    MIRBASE_SPECIES="mmu"
    ;;
  rat)
    ENSEMBL_RELEASE="112"
    ENSEMBL_SPECIES="rattus_norvegicus"
    UCSC_GENOME="rn7"
    FASTA_URL="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/${ENSEMBL_SPECIES}/dna/Rattus_norvegicus.mRatBN7.2.dna.primary_assembly.fa.gz"
    GTF_URL="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gtf/${ENSEMBL_SPECIES}/Rattus_norvegicus.mRatBN7.2.${ENSEMBL_RELEASE}.gtf.gz"
    MIRBASE_SPECIES="rno"
    ;;
  *)
    echo "ERROR: Unknown organism '${ORGANISM}'. Use: human | mouse | rat"
    exit 1
    ;;
esac

FASTA_GZ="${REF_DIR}/genome.fa.gz"
FASTA="${REF_DIR}/genome.fa"
GTF_GZ="${REF_DIR}/annotation.gtf.gz"
GTF="${REF_DIR}/annotation.gtf"
REFFLAT="${REF_DIR}/refFlat.txt"                  # for CIRCexplorer2
MIRBASE_MATURE="${REF_DIR}/mirbase_mature.fa"
MIRBASE_HAIRPIN="${REF_DIR}/mirbase_hairpin.fa"

echo "=== Downloading ${ORGANISM} reference (Ensembl ${ENSEMBL_RELEASE}) ==="

# ── 1. Genome FASTA ───────────────────────────────────────────────────
if [[ ! -f "${FASTA}" ]]; then
  echo "[1/5] Downloading genome FASTA..."
  wget -q --show-progress -O "${FASTA_GZ}" "${FASTA_URL}"
  pigz -d -p "${THREADS}" "${FASTA_GZ}"
  echo "  -> ${FASTA}"
else
  echo "[1/5] ${FASTA} already exists, skipping"
fi

# ── 2. GTF annotation ─────────────────────────────────────────────────
if [[ ! -f "${GTF}" ]]; then
  echo "[2/5] Downloading GTF..."
  wget -q --show-progress -O "${GTF_GZ}" "${GTF_URL}"
  pigz -d -p "${THREADS}" "${GTF_GZ}"
  echo "  -> ${GTF}"
else
  echo "[2/5] ${GTF} already exists, skipping"
fi

# ── 3. RefFlat (for CIRCexplorer2) ───────────────────────────────────
if [[ ! -f "${REFFLAT}" ]]; then
  echo "[3/5] Downloading refFlat for CIRCexplorer2..."
  wget -q -O - \
    "https://hgdownload.soe.ucsc.edu/goldenPath/${UCSC_GENOME}/database/refFlat.txt.gz" \
    | gunzip -c > "${REFFLAT}"
  echo "  -> ${REFFLAT}"
else
  echo "[3/5] ${REFFLAT} already exists, skipping"
fi

# ── 4. miRBase (for miRDeep2) ─────────────────────────────────────────
MIRBASE_VER="22"
if [[ ! -f "${MIRBASE_MATURE}" ]]; then
  echo "[4/5] Downloading miRBase v${MIRBASE_VER} mature/hairpin sequences..."
  wget -q -O - \
    "https://mirbase.org/ftp/${MIRBASE_VER}/mature.fa.gz" \
    | gunzip -c | grep -A1 "^>${MIRBASE_SPECIES}-" > "${MIRBASE_MATURE}"
  wget -q -O - \
    "https://mirbase.org/ftp/${MIRBASE_VER}/hairpin.fa.gz" \
    | gunzip -c | grep -A1 "^>${MIRBASE_SPECIES}-" > "${MIRBASE_HAIRPIN}"
  echo "  -> ${MIRBASE_MATURE}"
  echo "  -> ${MIRBASE_HAIRPIN}"
else
  echo "[4/5] miRBase files already exist, skipping"
fi

# ── 5. Transcript FASTA for Salmon ────────────────────────────────────
TRANSCRIPTS="${REF_DIR}/transcriptome.fa"
if [[ ! -f "${TRANSCRIPTS}" ]]; then
  echo "[5/5] Extracting transcript sequences from genome + GTF..."
  mamba-haining run -n rnaseq_env \
    gffread "${GTF}" -g "${FASTA}" -w "${TRANSCRIPTS}"
  echo "  -> ${TRANSCRIPTS}"
else
  echo "[5/5] ${TRANSCRIPTS} already exists, skipping"
fi

echo ""
echo "=== Reference download complete ==="
echo "  Genome FASTA : ${FASTA}"
echo "  GTF          : ${GTF}"
echo "  Transcriptome: ${TRANSCRIPTS}"
echo ""
echo "Next step: bash scripts/01_build_indices.sh"
