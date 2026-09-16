#!/usr/bin/env bash
# 01_build_indices.sh
# Build STAR, Salmon, Bowtie v1, and BWA genome indices.
# Run ONCE after downloading references (scripts/00_download_refs.sh).
# Each index type lives in its own subdirectory under ref/.
#
# Usage:
#   bash scripts/01_build_indices.sh            # all indices
#   bash scripts/01_build_indices.sh --star     # STAR only
#   bash scripts/01_build_indices.sh --salmon   # Salmon only
#
# Resource requirements (human genome):
#   STAR:   ~32 GB RAM, ~2 hours
#   Salmon: ~8  GB RAM, ~10 minutes
#   Bowtie: ~4  GB RAM, ~10 minutes
#   BWA:    ~8  GB RAM, ~2 hours

set -euo pipefail

DATA_BASE="${DATA_BASE:-/xdisk/haining/maarowosegbe/ngs101-pipeline}"
REF_DIR="${DATA_BASE}/ref"
FASTA="${REF_DIR}/genome.fa"
GTF="${REF_DIR}/annotation.gtf"
TRANSCRIPTS="${REF_DIR}/transcriptome.fa"
THREADS="${THREADS:-16}"

BUILD_STAR=true
BUILD_SALMON=true
BUILD_BOWTIE=true
BUILD_BWA=true

for arg in "$@"; do
  case "$arg" in
    --star)   BUILD_SALMON=false; BUILD_BOWTIE=false; BUILD_BWA=false ;;
    --salmon) BUILD_STAR=false; BUILD_BOWTIE=false; BUILD_BWA=false ;;
    --bowtie) BUILD_STAR=false; BUILD_SALMON=false; BUILD_BWA=false ;;
    --bwa)    BUILD_STAR=false; BUILD_SALMON=false; BUILD_BOWTIE=false ;;
  esac
done

module load apptainer 2>/dev/null || true  # HPC module

# ── Verify references exist ───────────────────────────────────────────
for f in "${FASTA}" "${GTF}"; do
  [[ -f "$f" ]] || { echo "ERROR: $f not found. Run scripts/00_download_refs.sh first."; exit 1; }
done

# ── 1. STAR genome index ──────────────────────────────────────────────
STAR_INDEX="${REF_DIR}/star_index"
if $BUILD_STAR && [[ ! -f "${STAR_INDEX}/Genome" ]]; then
  echo "=== Building STAR index (this takes ~2 hours for hg38) ==="
  mkdir -p "${STAR_INDEX}"
  conda run -n rnaseq_env \
    STAR \
      --runMode          genomeGenerate \
      --genomeDir        "${STAR_INDEX}" \
      --genomeFastaFiles "${FASTA}" \
      --sjdbGTFfile      "${GTF}" \
      --sjdbOverhang     149 \
      --runThreadN       "${THREADS}" \
      --genomeSAindexNbases 14
  echo "  -> ${STAR_INDEX}"
elif $BUILD_STAR; then
  echo "[STAR] Index already exists at ${STAR_INDEX}, skipping"
fi

# ── 2. Salmon transcript index ────────────────────────────────────────
SALMON_INDEX="${REF_DIR}/salmon_index"
if $BUILD_SALMON && [[ ! -d "${SALMON_INDEX}" ]]; then
  [[ -f "${TRANSCRIPTS}" ]] || { echo "ERROR: ${TRANSCRIPTS} missing. Re-run 00_download_refs.sh"; exit 1; }
  echo "=== Building Salmon index ==="
  conda run -n rnaseq_env \
    salmon index \
      -t "${TRANSCRIPTS}" \
      -d <(grep "^>" "${FASTA}" | cut -d ' ' -f 1 | tr -d '>') \
      -i "${SALMON_INDEX}" \
      --gencode \
      -p "${THREADS}"
  echo "  -> ${SALMON_INDEX}"
elif $BUILD_SALMON; then
  echo "[Salmon] Index already exists, skipping"
fi

# ── 3. Bowtie v1 index (for miRDeep2) ────────────────────────────────
BOWTIE_INDEX="${REF_DIR}/bowtie_index/genome"
if $BUILD_BOWTIE && [[ ! -f "${BOWTIE_INDEX}.1.ebwt" ]]; then
  echo "=== Building Bowtie v1 index (for miRDeep2) ==="
  mkdir -p "$(dirname "${BOWTIE_INDEX}")"
  conda run -n mirna_env \
    bowtie-build --threads "${THREADS}" "${FASTA}" "${BOWTIE_INDEX}"
  echo "  -> ${BOWTIE_INDEX}*"
elif $BUILD_BOWTIE; then
  echo "[Bowtie] Index already exists, skipping"
fi

# ── 4. BWA index (for CIRI3 circRNA detection) ────────────────────────
if $BUILD_BWA && [[ ! -f "${FASTA}.bwt" ]]; then
  echo "=== Building BWA index (for CIRI3, ~2 hours for hg38) ==="
  conda run -n ciri3_env \
    bwa index "${FASTA}"
  echo "  -> ${FASTA}.bwt"
elif $BUILD_BWA; then
  echo "[BWA] Index already exists, skipping"
fi

echo ""
echo "=== All indices built ==="
echo "  STAR   : ${STAR_INDEX}"
echo "  Salmon : ${SALMON_INDEX}"
echo "  Bowtie : ${BOWTIE_INDEX}*"
echo "  BWA    : ${FASTA}.bwt"
echo ""
echo "Next step: bash scripts/02_download_sra.sh  (if using SRA data)"
echo "  OR       bash scripts/03_qc_trim.sh       (if you have FASTQs)"
