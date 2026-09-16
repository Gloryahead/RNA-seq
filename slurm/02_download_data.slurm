#!/bin/bash
# --------------------------------------------------------------
### PART 1: Request resources
# --------------------------------------------------------------

#SBATCH --job-name=ngs101_download
#SBATCH --output=/xdisk/haining/maarowosegbe/RNA-seq/logs/slurm/%x-%j.out
#SBATCH --error=/xdisk/haining/maarowosegbe/RNA-seq/logs/slurm/%x-%j.err
#SBATCH --account=haining
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=24:00:00        # SRA downloads can be slow — give it a full day

# --------------------------------------------------------------
### PART 2: Download FASTQs from SRA
# --------------------------------------------------------------

echo "========== Download job started at: $(date) =========="

SCRIPT_BASE=/home/u11/maarowosegbe/RNA-seq
DATA_BASE=/xdisk/haining/maarowosegbe/RNA-seq

export THREADS=${SLURM_NTASKS}
export TMPDIR="${DATA_BASE}/tmp"
mkdir -p "${TMPDIR}"

set +eu; source ~/.bashrc; set -eu
export MAMBA_ROOT_PREFIX=/groups/haining/maarowosegbe/micromamba
export MAMBA_EXE=/opt/ohpc/pub/apps/micromamba/2.0.2-2/bin/micromamba
eval "$(micromamba shell hook --shell=bash)"
micromamba activate rnaseq_env

ACCESSION_FILE="${DATA_BASE}/data/sra_accessions.txt"
FASTQ_DIR="${DATA_BASE}/data/fastq"

if [[ ! -s "${ACCESSION_FILE}" ]]; then
  echo "ERROR: ${ACCESSION_FILE} is empty."
  echo "  Add your SRA accessions (one per line) and re-run this job."
  exit 1
fi

echo "Downloading SRA accessions from: ${ACCESSION_FILE}"
bash "${SCRIPT_BASE}/scripts/02_download_sra.sh" \
  "${ACCESSION_FILE}" \
  "${FASTQ_DIR}"

echo "========== Download job finished at: $(date) =========="
