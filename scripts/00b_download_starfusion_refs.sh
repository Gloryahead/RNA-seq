#!/usr/bin/env bash
# 00b_download_starfusion_refs.sh
# Download the STAR-Fusion CTAT plug-n-play genome library for GRCh38.
# Required before running STAR-Fusion via Apptainer (fusiongene_env step).
# Run ONCE; the library is ~30 GB extracted.
#
# Prerequisites:
#   - Apptainer SIF already pulled:
#       apptainer pull --dir /xdisk/haining/maarowosegbe/apptainer_images \
#           docker://trinityctat/starfusion:latest
#
# Usage:
#   bash scripts/00b_download_starfusion_refs.sh
#
# Output: /xdisk/haining/maarowosegbe/refs/GRCh38_gencode_v44_CTAT_lib_Apr012024.plug-n-play/

set -euo pipefail

# ── Paths ──────────────────────────────────────────────────────────────────────
REFS_DIR="/xdisk/haining/maarowosegbe/refs"
CTAT_TARBALL="GRCh38_gencode_v44_CTAT_lib_Apr012024.plug-n-play.tar.gz"
CTAT_URL="https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/${CTAT_TARBALL}"
CTAT_DIR="${REFS_DIR}/GRCh38_gencode_v44_CTAT_lib_Apr012024.plug-n-play"

mkdir -p "${REFS_DIR}"

# ── Download ───────────────────────────────────────────────────────────────────
if [[ -d "${CTAT_DIR}" ]]; then
    echo "CTAT library already exists at ${CTAT_DIR} — skipping download."
    exit 0
fi

echo "=== Downloading STAR-Fusion CTAT genome library (~8 GB compressed) ==="
wget -q --show-progress -O "${REFS_DIR}/${CTAT_TARBALL}" "${CTAT_URL}"

# ── Extract ────────────────────────────────────────────────────────────────────
echo "=== Extracting (~30 GB) ==="
tar -xzf "${REFS_DIR}/${CTAT_TARBALL}" -C "${REFS_DIR}"

# ── Clean up tarball ───────────────────────────────────────────────────────────
rm "${REFS_DIR}/${CTAT_TARBALL}"
echo "=== Done — tarball removed to save space ==="

echo ""
echo "CTAT library: ${CTAT_DIR}"
echo ""
echo "Test STAR-Fusion with:"
echo "  SIF=/xdisk/haining/maarowosegbe/apptainer_images/starfusion_latest.sif"
echo "  apptainer exec \\"
echo "      --bind /xdisk/haining/maarowosegbe:/xdisk/haining/maarowosegbe \\"
echo "      \$SIF STAR-Fusion \\"
echo "      --left_fq R1.fastq.gz \\"
echo "      --right_fq R2.fastq.gz \\"
echo "      --genome_lib_dir ${CTAT_DIR} \\"
echo "      --CPU 8 \\"
echo "      --output_dir star_fusion_output"
