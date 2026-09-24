#!/bin/bash
# ──────────────────────────────────────────────────────────────────────
# StringTie transcript quantification — Chunlong Ma RNA-seq (9040-0001)
# 3 human samples, GRCh38/GENCODE v44, paired-end 151 bp
#
# Workflow:
#   Step 1: guided assembly per sample       → assembly/{sample}.gtf
#   Step 2: merge all assemblies             → merged.gtf
#   Step 3: re-estimate expression (-e -B)   → quant/{sample}/
#            gene_abund.tab  — gene-level TPM
#            t_data.ctab     — Ballgown transcript table
# ──────────────────────────────────────────────────────────────────────

#SBATCH --job-name=stringtie_chunlong
#SBATCH --output=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/logs/stringtie_%j.out
#SBATCH --error=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/logs/stringtie_%j.err
#SBATCH --account=haining
#SBATCH --partition=standard
#SBATCH --ntasks=8
#SBATCH --time=04:00:00

# ── Paths ──────────────────────────────────────────────────────────────
BAM_DIR=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/9040-0001-human-analysis/02.BamFiles
GTF=/xdisk/haining/maarowosegbe/refs/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf
OUTDIR=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/9040-0001-human-analysis/03.StringTie
LOGDIR=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/logs
THREADS=${SLURM_NTASKS}

SAMPLES=(
    S9040-0001-0001
    S9040-0001-0002
    S9040-0001-0003
)

# ── Environment ────────────────────────────────────────────────────────
export MAMBA_ROOT_PREFIX=/groups/haining/maarowosegbe/micromamba
eval "$(/opt/ohpc/pub/apps/micromamba/2.0.2-2/bin/micromamba shell hook --shell bash)"
micromamba activate transcript_splicing_env

mkdir -p "${OUTDIR}/assembly" "${OUTDIR}/quant" "${LOGDIR}"

echo "========== StringTie started at: $(date) =========="

# ── Step 1: guided assembly per sample ────────────────────────────────
echo "--- Step 1: per-sample assembly ---"
for SAMPLE in "${SAMPLES[@]}"; do
    BAM="${BAM_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam"
    OUT_GTF="${OUTDIR}/assembly/${SAMPLE}.gtf"
    echo "  Assembling ${SAMPLE}..."
    stringtie "${BAM}" \
        -G "${GTF}" \
        -o "${OUT_GTF}" \
        -p "${THREADS}" \
        2>"${LOGDIR}/${SAMPLE}_assemble.log"
    echo "  Done: ${OUT_GTF}"
done

# ── Step 2: merge all sample GTFs ─────────────────────────────────────
echo "--- Step 2: merge assemblies ---"
ASSEMBLY_GTFS=$(printf '%s\n' "${SAMPLES[@]}" | \
    sed "s|.*|${OUTDIR}/assembly/&.gtf|")

stringtie --merge \
    -G "${GTF}" \
    -o "${OUTDIR}/merged.gtf" \
    ${ASSEMBLY_GTFS} \
    2>"${LOGDIR}/stringtie_merge.log"
echo "  Merged GTF: ${OUTDIR}/merged.gtf"

# ── Step 3: re-estimate expression against merged GTF ─────────────────
echo "--- Step 3: expression quantification ---"
for SAMPLE in "${SAMPLES[@]}"; do
    BAM="${BAM_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam"
    SAMPLE_OUTDIR="${OUTDIR}/quant/${SAMPLE}"
    mkdir -p "${SAMPLE_OUTDIR}"
    echo "  Quantifying ${SAMPLE}..."
    stringtie "${BAM}" \
        -G "${OUTDIR}/merged.gtf" \
        -e -B \
        -o "${SAMPLE_OUTDIR}/quant.gtf" \
        -A "${SAMPLE_OUTDIR}/gene_abund.tab" \
        -p "${THREADS}" \
        2>"${LOGDIR}/${SAMPLE}_quant.log"
    echo "  Done: ${SAMPLE_OUTDIR}/"
done

echo "========== StringTie finished at: $(date) =========="
echo ""
echo "Outputs:"
echo "  Merged transcriptome : ${OUTDIR}/merged.gtf"
echo "  Gene TPMs per sample : ${OUTDIR}/quant/{sample}/gene_abund.tab"
echo "  Ballgown tables      : ${OUTDIR}/quant/{sample}/t_data.ctab"
