#!/usr/bin/env bash
# 07_convert_formats.sh
# Common BAM/FASTQ format conversions needed by specialized analysis tools.
#
#   BAM → bigWig     (for genome browser tracks; deepTools)
#   BAM → BED        (for BEDTools-based analyses)
#   BAM → FASTQ      (recover reads from BAM, e.g. for re-alignment)
#   FASTQ → FASTA    (some tools require FASTA input)
#   GTF → BED12      (for RSeQC, various splicing tools)
#   GTF → genePred   (for CIRCexplorer2 refFlat)
#   SAM → sorted BAM (if you have raw SAM from a manual step)
#
# Usage:
#   bash scripts/07_convert_formats.sh bam2bw         # BAM → bigWig
#   bash scripts/07_convert_formats.sh bam2bed        # BAM → BED
#   bash scripts/07_convert_formats.sh bam2fastq      # BAM → FASTQ
#   bash scripts/07_convert_formats.sh gtf2bed        # GTF → BED12
#   bash scripts/07_convert_formats.sh sam2bam        # sort + index SAM

set -euo pipefail

COMMAND="${1:-help}"
BAM_DIR="results/bam"
THREADS="${THREADS:-8}"
CHROM_SIZES="ref/chrom.sizes"
GTF="ref/annotation.gtf"

# ── Generate chrom.sizes if needed ────────────────────────────────────
make_chrom_sizes() {
  if [[ ! -f "${CHROM_SIZES}" ]]; then
    echo "  Generating chrom.sizes..."
    conda run -n rnaseq_env \
      samtools view -H "${BAM_DIR}/$(ls "${BAM_DIR}"/*.bam | head -1 | xargs basename)" \
      | grep "^@SQ" \
      | awk '{gsub("SN:|LN:","",$2" "$3); print $2"\t"$3}' \
      > "${CHROM_SIZES}"
  fi
}

# ── BAM → bigWig (normalized coverage track for genome browser) ───────
bam2bw() {
  mkdir -p results/bigwig
  for bam in "${BAM_DIR}"/*.bam; do
    sample=$(basename "${bam}" _Aligned.sortedByCoord.out.bam)
    out="results/bigwig/${sample}.bw"
    [[ -f "${out}" ]] && { echo "  ${sample}: bigWig exists, skipping"; continue; }
    echo "  ${sample}: BAM → bigWig (RPKM normalized)..."
    conda run -n rnaseq_env \
      bamCoverage \
        -b "${bam}" \
        -o "${out}" \
        --normalizeUsing RPKM \
        --binSize 10 \
        --smoothLength 50 \
        -p "${THREADS}" \
        2>"logs/convert/${sample}_bam2bw.log"
    echo "    -> ${out}"
  done
}

# ── BAM → BED (for BEDTools analyses) ────────────────────────────────
bam2bed() {
  mkdir -p results/bed
  for bam in "${BAM_DIR}"/*.bam; do
    sample=$(basename "${bam}" _Aligned.sortedByCoord.out.bam)
    out="results/bed/${sample}.bed"
    [[ -f "${out}" ]] && { echo "  ${sample}: BED exists, skipping"; continue; }
    echo "  ${sample}: BAM → BED..."
    conda run -n rnaseq_env \
      bedtools bamtobed -i "${bam}" | sort -k1,1 -k2,2n > "${out}"
    echo "    -> ${out}"
  done
}

# ── BAM → FASTQ (recover reads) ───────────────────────────────────────
bam2fastq() {
  mkdir -p results/recovered_fastq
  for bam in "${BAM_DIR}"/*.bam; do
    sample=$(basename "${bam}" _Aligned.sortedByCoord.out.bam)
    out_r1="results/recovered_fastq/${sample}_R1.fastq.gz"
    [[ -f "${out_r1}" ]] && { echo "  ${sample}: FASTQ exists, skipping"; continue; }
    echo "  ${sample}: BAM → FASTQ..."
    conda run -n rnaseq_env \
      samtools sort -n -@ "${THREADS}" -o /tmp/${sample}_namesorted.bam "${bam}"
    conda run -n rnaseq_env \
      samtools fastq -@ "${THREADS}" \
        -1 results/recovered_fastq/${sample}_R1.fastq.gz \
        -2 results/recovered_fastq/${sample}_R2.fastq.gz \
        -0 /dev/null \
        -s /dev/null \
        /tmp/${sample}_namesorted.bam \
        2>"logs/convert/${sample}_bam2fq.log"
    rm -f /tmp/${sample}_namesorted.bam
    echo "    -> results/recovered_fastq/${sample}_R{1,2}.fastq.gz"
  done
}

# ── GTF → BED12 (for RSeQC / infer_experiment.py) ────────────────────
gtf2bed() {
  OUT="ref/annotation.bed12"
  [[ -f "${OUT}" ]] && { echo "  ${OUT} exists, skipping"; return; }
  echo "  GTF → BED12..."
  conda run -n rnaseq_env \
    gtfToGenePred "${GTF}" /tmp/annotation.genePred
  conda run -n rnaseq_env \
    genePredToBed /tmp/annotation.genePred "${OUT}"
  echo "  -> ${OUT}"
}

# ── SAM → sorted BAM + index ──────────────────────────────────────────
sam2bam() {
  mkdir -p logs/convert
  for sam in results/bam/*.sam; do
    [[ -f "${sam}" ]] || continue
    sample=$(basename "${sam}" .sam)
    out="${BAM_DIR}/${sample}.sorted.bam"
    echo "  ${sample}: SAM → sorted BAM..."
    conda run -n rnaseq_env \
      samtools sort -@ "${THREADS}" -o "${out}" "${sam}"
    conda run -n rnaseq_env \
      samtools index "${out}"
    rm "${sam}"
    echo "    -> ${out}"
  done
}

# ── Dispatch ──────────────────────────────────────────────────────────
mkdir -p logs/convert

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
