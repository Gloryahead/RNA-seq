#!/usr/bin/env bash
# 04_align_star.sh
# STAR 2-pass alignment for all samples.
# Produces: sorted BAM + BAM index + chimeric junctions (for circRNA/fusion).
# Reads trimmed FASTQs from results/trimmed/ (output of 03_qc_trim.sh).
#
# Usage:
#   bash scripts/04_align_star.sh
#   bash scripts/04_align_star.sh sample_id    # single sample
#   THREADS=16 bash scripts/04_align_star.sh
#
# Key STAR flags:
#   --chimOutType Junctions   needed for STAR-Fusion and CIRCexplorer2
#   --outSAMtype BAM SortedByCoordinate   sorted BAM directly

set -euo pipefail

SAMPLES_TSV="config/samples.tsv"
STAR_INDEX="ref/star_index"
GTF="ref/annotation.gtf"
INDIR="results/trimmed"
OUTDIR="results/bam"
THREADS="${THREADS:-16}"
SINGLE_SAMPLE="${1:-}"

mkdir -p "${OUTDIR}" logs/star

[[ -f "${STAR_INDEX}/Genome" ]] || {
  echo "ERROR: STAR index not found at ${STAR_INDEX}. Run scripts/01_build_indices.sh first."
  exit 1
}

align_sample() {
  local sample="$1"
  local r1="$2"
  local r2="${3:-}"

  local bam="${OUTDIR}/${sample}_Aligned.sortedByCoord.out.bam"
  if [[ -f "${bam}" ]]; then
    echo "  ${sample}: BAM already exists, skipping"
    return
  fi

  echo ""
  echo "═══ STAR: ${sample} ═══"

  # Detect trimmed FASTQ names (Trim Galore adds _val_1/_val_2 suffix)
  local star_r1 star_r2=""
  # Paired-end
  if [[ -n "${r2}" ]]; then
    star_r1="${INDIR}/${sample}_val_1.fq.gz"
    star_r2="${INDIR}/${sample}_val_2.fq.gz"
    # Fallback to original naming if Trim Galore wasn't used
    [[ -f "${star_r1}" ]] || star_r1="${r1}"
    [[ -f "${star_r2}" ]] || star_r2="${r2}"
    local reads_arg="${star_r1} ${star_r2}"
  else
    star_r1="${INDIR}/${sample}_trimmed.fq.gz"
    [[ -f "${star_r1}" ]] || star_r1="${r1}"
    local reads_arg="${star_r1}"
  fi

  echo "  R1: ${star_r1}"
  [[ -n "${star_r2}" ]] && echo "  R2: ${star_r2}"

  conda run -n rnaseq_env \
    STAR \
      --runMode            alignReads \
      --runThreadN         "${THREADS}" \
      --genomeDir          "${STAR_INDEX}" \
      --sjdbGTFfile        "${GTF}" \
      --readFilesIn        ${reads_arg} \
      --readFilesCommand   zcat \
      --outSAMtype         BAM SortedByCoordinate \
      --outSAMattributes   NH HI NM MD AS XS \
      --outSAMstrandField  intronMotif \
      --outFilterType      BySJout \
      --outFilterMultimapNmax 20 \
      --alignSJoverhangMin   8 \
      --alignSJDBoverhangMin 1 \
      --outFilterMismatchNmax 999 \
      --outFilterMismatchNoverReadLmax 0.04 \
      --alignIntronMin     20 \
      --alignIntronMax     1000000 \
      --alignMatesGapMax   1000000 \
      --chimSegmentMin     12 \
      --chimJunctionOverhangMin 12 \
      --chimOutType        Junctions \
      --outWigType         None \
      --outFileNamePrefix  "${OUTDIR}/${sample}_" \
      2>"logs/star/${sample}.log"

  echo "  Indexing BAM..."
  conda run -n rnaseq_env \
    samtools index -@ 4 "${bam}"

  echo "  Flagstat:"
  conda run -n rnaseq_env \
    samtools flagstat "${bam}" | sed 's/^/    /'

  echo "  Done: ${bam}"
}

# ── Read samples.tsv ─────────────────────────────────────────────────
{
  read -r _header
  while IFS=$'\t' read -r sample group r1 r2 batch || [[ -n "$sample" ]]; do
    [[ -z "${sample}" || "${sample}" =~ ^# ]] && continue
    [[ -n "${SINGLE_SAMPLE}" && "${sample}" != "${SINGLE_SAMPLE}" ]] && continue
    align_sample "${sample}" "${r1}" "${r2}"
  done
} < "${SAMPLES_TSV}"

echo ""
echo "=== STAR alignment complete ==="
echo "  BAMs in: ${OUTDIR}"
echo "  Chimeric junctions: ${OUTDIR}/*_Chimeric.out.junction"
echo ""
echo "Next step: bash scripts/05_quantify.sh"
