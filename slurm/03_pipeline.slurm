#!/bin/bash
# --------------------------------------------------------------
### PART 1: Request resources
# --------------------------------------------------------------

#SBATCH --job-name=ngs101_pipeline
#SBATCH --output=/xdisk/haining/maarowosegbe/RNA-seq/logs/slurm/%x-%j.out
#SBATCH --error=/xdisk/haining/maarowosegbe/RNA-seq/logs/slurm/%x-%j.err
#SBATCH --account=haining
#SBATCH --partition=standard

### The SLURM job itself only runs the Snakemake controller (1 CPU, low RAM).
### Snakemake then submits individual rules as child SLURM jobs automatically.
### Each child job requests its own CPUs/memory (defined in the SLURM profile).
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=240:00:00       # controller must stay alive for the whole run

# --------------------------------------------------------------
### PART 2: Run Snakemake pipeline
# --------------------------------------------------------------

echo "========== Pipeline started at: $(date) =========="
echo "  Node    : $(hostname)"
echo "  Job ID  : ${SLURM_JOB_ID}"

SCRIPT_BASE=/home/u11/maarowosegbe/RNA-seq
DATA_BASE=/xdisk/haining/maarowosegbe/RNA-seq

export DATA_BASE

set +eu; source ~/.bashrc; set -eu
export MAMBA_ROOT_PREFIX=/groups/haining/maarowosegbe/micromamba
export MAMBA_EXE=/opt/ohpc/pub/apps/micromamba/2.0.2-2/bin/micromamba
eval "$(micromamba shell hook --shell=bash)"
micromamba activate snakemake_env   # see NOTE below

echo ""
echo "--- Snakemake dry run (DAG check) ---"
snakemake \
  --snakefile "${SCRIPT_BASE}/workflow/Snakefile" \
  --configfile "${SCRIPT_BASE}/config/config.yaml" \
  --profile "${SCRIPT_BASE}/workflow/profiles/slurm" \
  --dry-run \
  --quiet

echo ""
echo "--- Full pipeline run ---"
snakemake \
  --snakefile "${SCRIPT_BASE}/workflow/Snakefile" \
  --configfile "${SCRIPT_BASE}/config/config.yaml" \
  --profile "${SCRIPT_BASE}/workflow/profiles/slurm" \
  --use-conda \
  --conda-frontend mamba \
  --conda-prefix "${DATA_BASE}/.snakemake/conda" \
  --rerun-incomplete \
  --keep-going \
  --latency-wait 60 \
  --jobs 100 \
  --printshellcmds \
  2>&1 | tee "${DATA_BASE}/logs/slurm/snakemake_${SLURM_JOB_ID}.log"

echo ""
echo "--- MultiQC summary ---"
micromamba activate rnaseq_env
bash "${SCRIPT_BASE}/scripts/06_multiqc.sh"

echo "========== Pipeline finished at: $(date) =========="

# NOTE: create a minimal snakemake_env if you don't have one:
#   mamba-haining create -n snakemake_env -c conda-forge -c bioconda \
#     snakemake=8.20 python=3.11 -y
