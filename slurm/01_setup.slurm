#!/bin/bash
# --------------------------------------------------------------
### PART 1: Request resources
# --------------------------------------------------------------

### Job name
#SBATCH --job-name=ngs101_setup

### Output log — goes to xdisk so it survives after the job
#SBATCH --output=/xdisk/haining/maarowosegbe/RNA-seq/logs/slurm/%x-%j.out
#SBATCH --error=/xdisk/haining/maarowosegbe/RNA-seq/logs/slurm/%x-%j.err

### PI group
#SBATCH --account=haining

### Partition (standard is fine for setup; use high_priority if available)
#SBATCH --partition=standard

### Resources — index building is memory-heavy (STAR needs ~32 GB for hg38)
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=4gb      # 16 × 4gb = 64 GB total

### Wall time — downloading + indexing can take 4–6 hours
#SBATCH --time=12:00:00

# --------------------------------------------------------------
### PART 2: Setup
# --------------------------------------------------------------

echo "========== Setup job started at: $(date) =========="
echo "  Node    : $(hostname)"
echo "  Job ID  : ${SLURM_JOB_ID}"

# Two-path pattern:
#   SCRIPT_BASE → home (quota-safe, version-controlled code)
#   DATA_BASE   → xdisk (large storage: refs, data, results, logs)
SCRIPT_BASE=/home/u11/maarowosegbe/RNA-seq
DATA_BASE=/xdisk/haining/maarowosegbe/RNA-seq

export THREADS=${SLURM_NTASKS}
export DATA_BASE

# Ensure xdisk tree exists
bash "${SCRIPT_BASE}/scripts/setup_xdisk.sh"

# Activate the core conda environment
set +eu
source ~/.bashrc
set -eu
eval "$(mamba-haining shell hook --shell=bash)"
mamba-haining activate rnaseq_env

echo ""
echo "--- Step 1: Download references ---"
ORGANISM="${ORGANISM:-human}"
bash "${SCRIPT_BASE}/scripts/00_download_refs.sh" "${ORGANISM}" "${THREADS}"

echo ""
echo "--- Step 2: Build indices ---"
bash "${SCRIPT_BASE}/scripts/01_build_indices.sh"

echo ""
echo "========== Setup job finished at: $(date) =========="
echo "Next: sbatch ${SCRIPT_BASE}/slurm/03_pipeline.slurm"
