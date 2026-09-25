#!/bin/bash
#SBATCH --job-name=rmats_summary
#SBATCH --output=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/logs/rmats_summary_%j.out
#SBATCH --error=/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/logs/rmats_summary_%j.err
#SBATCH --account=haining
#SBATCH --partition=standard
#SBATCH --ntasks=1
#SBATCH --time=00:30:00

export MAMBA_ROOT_PREFIX=/groups/haining/maarowosegbe/micromamba
eval "$(/opt/ohpc/pub/apps/micromamba/2.0.2-2/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq_r_env

cd ~/RNA-seq
Rscript workflow/scripts/rmats_summary.R
