#!/bin/bash
# --------------------------------------------------------------
### Chunlong Ma RNA-seq — 9040-0001 — Human hg38
### n=1 per group; use results for exploration only
# --------------------------------------------------------------

#SBATCH --job-name=chunlong_rnaseq
#SBATCH --output=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/logs/slurm/%x-%j.out
#SBATCH --error=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/logs/slurm/%x-%j.err
#SBATCH --account=haining
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=240:00:00

echo "========== Pipeline started at: $(date) =========="
echo "  Node    : $(hostname)"
echo "  Job ID  : ${SLURM_JOB_ID}"

SCRIPT_BASE=/home/u11/maarowosegbe/RNA-seq
LOG_BASE=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq

mkdir -p "${LOG_BASE}/logs/slurm"

export MAMBA_ROOT_PREFIX=/groups/haining/maarowosegbe/micromamba
export MAMBA_EXE=/opt/ohpc/pub/apps/micromamba/2.0.2-2/bin/micromamba
MC="${MAMBA_EXE}"

echo ""
echo "--- Snakemake dry run ---"
"${MC}" run -n snakemake_env \
  snakemake \
    --snakefile "${SCRIPT_BASE}/workflow/Snakefile" \
    --configfile "${SCRIPT_BASE}/config/config_chunlong.yaml" \
    --profile "${SCRIPT_BASE}/workflow/profiles/slurm" \
    --dry-run \
    --quiet

echo ""
echo "--- Full pipeline run ---"
"${MC}" run -n snakemake_env \
  snakemake \
    --snakefile "${SCRIPT_BASE}/workflow/Snakefile" \
    --configfile "${SCRIPT_BASE}/config/config_chunlong.yaml" \
    --profile "${SCRIPT_BASE}/workflow/profiles/slurm" \
    --rerun-incomplete \
    --keep-going \
    --latency-wait 60 \
    --jobs 100 \
    --printshellcmds \
  2>&1 | tee "${LOG_BASE}/logs/slurm/snakemake_${SLURM_JOB_ID}.log"

echo "========== Pipeline finished at: $(date) =========="
