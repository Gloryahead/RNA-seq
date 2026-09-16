#!/usr/bin/env bash
# setup_xdisk.sh
# Create the full xdisk directory tree for the pipeline.
# Run ONCE on the HPC before any other script.
#
# Usage:
#   bash scripts/setup_xdisk.sh
#
# After this, xdisk layout will be:
#   /xdisk/haining/maarowosegbe/RNA-seq/
#     ref/          ← genome FASTA, GTF, indices
#     data/fastq/   ← raw FASTQs (SRA or your own)
#     results/      ← all pipeline outputs
#     logs/         ← all SLURM + tool logs
#     containers/   ← Apptainer SIF files (if using --singularity)

set -euo pipefail

DATA_BASE="${DATA_BASE:-/xdisk/haining/maarowosegbe/RNA-seq}"

echo "Creating xdisk directory structure under: ${DATA_BASE}"

mkdir -p \
  "${DATA_BASE}/ref" \
  "${DATA_BASE}/data/fastq" \
  "${DATA_BASE}/data/sra_accessions" \
  "${DATA_BASE}/results/fastqc" \
  "${DATA_BASE}/results/trimmed" \
  "${DATA_BASE}/results/bam" \
  "${DATA_BASE}/results/salmon" \
  "${DATA_BASE}/results/counts" \
  "${DATA_BASE}/results/multiqc" \
  "${DATA_BASE}/results/deg" \
  "${DATA_BASE}/results/pathways" \
  "${DATA_BASE}/results/splicing/spladder" \
  "${DATA_BASE}/results/splicing/suppa2" \
  "${DATA_BASE}/results/circrna/circexplorer2" \
  "${DATA_BASE}/results/circrna/ciri3" \
  "${DATA_BASE}/results/mirna" \
  "${DATA_BASE}/results/fusion/starfusion" \
  "${DATA_BASE}/results/fusion/fusioncatcher" \
  "${DATA_BASE}/results/viral/esviritu" \
  "${DATA_BASE}/results/crispr/counts" \
  "${DATA_BASE}/results/crispr/results" \
  "${DATA_BASE}/results/crispr/figures" \
  "${DATA_BASE}/results/networks/wgcna" \
  "${DATA_BASE}/results/networks/genie3" \
  "${DATA_BASE}/results/networks/masterreg" \
  "${DATA_BASE}/results/cancer_subtype" \
  "${DATA_BASE}/results/bigwig" \
  "${DATA_BASE}/results/normalization" \
  "${DATA_BASE}/logs/slurm" \
  "${DATA_BASE}/logs/star" \
  "${DATA_BASE}/logs/trim" \
  "${DATA_BASE}/logs/salmon" \
  "${DATA_BASE}/logs/featurecounts" \
  "${DATA_BASE}/logs/mirdeep2" \
  "${DATA_BASE}/logs/circexplorer2" \
  "${DATA_BASE}/logs/ciri3" \
  "${DATA_BASE}/logs/starfusion" \
  "${DATA_BASE}/logs/mageck" \
  "${DATA_BASE}/logs/networks" \
  "${DATA_BASE}/containers"

# Drop a placeholder accessions file so the user knows where to put their list
cat > "${DATA_BASE}/data/sra_accessions.txt" <<'EOF'
# SRA accession list — one accession per line
# Example:
# SRR1234567
# SRR7654321
EOF

echo ""
echo "Done. Directory tree created at: ${DATA_BASE}"
echo ""
echo "Next steps:"
echo "  1. Edit config/config.yaml  (genome paths → ${DATA_BASE}/ref/...)"
echo "  2. Edit config/samples.tsv  (fastq paths → ${DATA_BASE}/data/fastq/...)"
echo "  3. sbatch slurm/01_setup.slurm    (download refs + build indices)"
echo "  4. sbatch slurm/02_pipeline.slurm (full pipeline)"
