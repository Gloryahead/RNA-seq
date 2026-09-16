#!/usr/bin/env bash
# 07_convert_formats.sh
# Common BAM/FASTQ format conversions needed by specialized analysis tools.
#
#   BAM → bigWig     (for genome browser tracks; deepTools)
#   BAM → BED        (for BEDTools-based analyses)
#   BAM → FASTQ      (recover reads from BAM, e.g. for re-alignment)
#   GTF → BED12      (for RSeQC, various splicing tools)
#   GTF → genePred   (for CIRCexplorer2 refFlat)
#   SAM → sorted BAM (if you have raw SAM from a manual step)
#
# Usage (on HPC: wrap in srun or sbatch):
#   bash scripts/07_convert_formats.sh bam2bw
#   bash scripts/07_convert_formats.sh bam2bed
#   bash scripts/07_convert_formats.sh bam2fastq
#   bash scripts/07_convert_formats.sh gtf2bed
#   bash scripts/07_convert_formats.sh sam2bam

set -euo pipefail

SCRIPT_BASE="${SCRIPT_BASE:-/home/u11/maarowosegbe/RNA-seq}"
DATA_BASE="${DATA_BASE:-/xdisk/haining/maarowosegbe/RNA-seq}"
COMMAND="${1:-help}"
BAM_DIR="${DATA_BASE}/results/bam"
THREADS="${THREADS:-8}"
CHROM_SIZES="${DATA_BASE}/ref/chrom.sizes"
GTF="${DATA_BASE}/ref/annotation.gtf"
TMP_DIR="${DATA_BASE}/tmp"

mkdir -p "${DATA_BASE}/logs/convert" "${TMP_DIR}"

# ── Generate chrom.sizes if needed ────────────────────────────────────
make_chrom_sizes() {
  if [[ ! -f "${CHROM_SIZES}" ]]; then
    echo "  Generating chrom.sizes..."
    mamba-haining run -n rnaseq_env \
      samtools view -H "${BAM_DIR}/$(ls "${BAM_DIR}"/*.bam | head -1 | xargs basename)" \
      | grep "^@SQ" \
      | awk '{gsub("SN:|LN:","",$2" "$3); print $2"\t"$3}' \
      > "${CHROM_SIZES}"
  fi
}

# ── BAM → bigWig (normalized coverage track for genome browser) ───────
bam2bw() {
  mkdir -p "${DATA_BASE}/results/bigwig"
  for bam in "${BAM_DIR}"/*.bam; do
    sample=$(basename "${bam}" _Aligned.sortedByCoord.out.bam)
    out="${DATA_BASE}/results/bigwig/${sample}.bw"
    [[ -f "${out}" ]] && { echo "  ${sample}: bigWig exists, skipping"; continue; }
    echo "  ${sample}: BAM → bigWig (RPKM normalized)..."
    mamba-haining run -n rnaseq_env \
      bamCoverage \
        -b "${bam}" \
        -o "${out}" \
        --normalizeUsing RPKM \
        --binSize 10 \
        --smoothLength 50 \
        -p "${THREADS}" \
        2>"${DATA_BASE}/logs/convert/${sample}_bam2bw.log"
    echo "    -> ${out}"
  done
}

# ── BAM → BED (for BEDTools analyses) ────────────────────────────────
bam2bed() {
  mkdir -p "${DATA_BASE}/results/bed"
  for bam in "${BAM_DIR}"/*.bam; do
    sample=$(basename "${bam}" _Aligned.sortedByCoord.out.bam)
    out="${DATA_BASE}/results/bed/${sample}.bed"
    [[ -f "${out}" ]] && { echo "  ${sample}: BED exists, skipping"; continue; }
    echo "  ${sample}: BAM → BED..."
    mamba-haining run -n rnaseq_env \
      bedtools bamtobed -i "${bam}" | sort -k1,1 -k2,2n > "${out}"
    echo "    -> ${out}"
  done
}

# ── BAM → FASTQ (recover reads) ───────────────────────────────────────
bam2fastq() {
  mkdir -p "${DATA_BASE}/results/recovered_fastq"
  for bam in "${BAM_DIR}"/*.bam; do
    sample=$(basename "${bam}" _Aligned.sortedByCoord.out.bam)
    out_r1="${DATA_BASE}/results/recovered_fastq/${sample}_R1.fastq.gz"
    [[ -f "${out_r1}" ]] && { echo "  ${sample}: FASTQ exists, skipping"; continue; }
    echo "  ${sample}: BAM → FASTQ..."
    mamba-haining run -n rnaseq_env \
      samtools sort -n -@ "${THREADS}" -o "${TMP_DIR}/${sample}_namesorted.bam" "${bam}"
    mamba-haining run -n rnaseq_env \
      samtools fastq -@ "${THREADS}" \
        -1 "${DATA_BASE}/results/recovered_fastq/${sample}_R1.fastq.gz" \
        -2 "${DATA_BASE}/results/recovered_fastq/${sample}_R2.fastq.gz" \
        -0 /dev/null \
        -s /dev/null \
        "${TMP_DIR}/${sample}_namesorted.bam" \
        2>"${DATA_BASE}/logs/convert/${sample}_bam2fq.log"
    rm -f "${TMP_DIR}/${sample}_namesorted.bam"
    echo "    -> ${DATA_BASE}/results/recovered_fastq/${sample}_R{1,2}.fastq.gz"
  done
}

# ── GTF → BED12 (for RSeQC / infer_experiment.py) ────────────────────
gtf2bed() {
  OUT="${DATA_BASE}/ref/annotation.bed12"
  [[ -f "${OUT}" ]] && { echo "  ${OUT} exists, skipping"; return; }
  echo "  GTF → BED12..."
  mamba-haining run -n rnaseq_env \
    gtfToGenePred "${GTF}" "${TMP_DIR}/annotation.genePred"
  mamba-haining run -n rnaseq_env \
    genePredToBed "${TMP_DIR}/annotation.genePred" "${OUT}"
  echo "  -> ${OUT}"
}

# ── SAM → sorted BAM + index ──────────────────────────────────────────
sam2bam() {
  for sam in "${BAM_DIR}"/*.sam; do
    [[ -f "${sam}" ]] || continue
    sample=$(basename "${sam}" .sam)
    out="${BAM_DIR}/${sample}.sorted.bam"
    echo "  ${sample}: SAM → sorted BAM..."
    mamba-haining run -n rnaseq_env \
      samtools sort -@ "${THREADS}" -o "${out}" "${sam}"
    mamba-haining run -n rnaseq_env \
      samtools index "${out}"
    rm "${sam}"
    echo "    -> ${out}"
  done
}

# ── Dispatch ──────────────────────────────────────────────────────────
case "${COMMAND}" in
  bam2bw)    make_chrom_sizes; bam2bw   ;;
  bam2bed)   bam2bed   ;;
  bam2fastq) bam2fastq ;;
  gtf2bed)   gtf2bed   ;;
  sam2bam)   sam2bam   ;;
  help|*)
    echo "Usage: bash scripts/07_convert_formats.sh <command>"
    echo ""
    echo "Commands:"
    echo "  bam2bw      BAM → normalized bigWig (for IGV / UCSC browser)"
    echo "  bam2bed     BAM → BED (for BEDTools operations)"
    echo "  bam2fastq   BAM → FASTQ (recover reads for re-alignment)"
    echo "  gtf2bed     GTF → BED12 (for RSeQC infer_experiment.py)"
    echo "  sam2bam     SAM → sorted, indexed BAM"
    ;;
esac
