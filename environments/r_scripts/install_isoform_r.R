#!/usr/bin/env Rscript
# install_isoform_r.R
# Run after activating isoform_env to install Bioconductor packages required
# for tximport + DRIMSeq DTU analysis.
# Usage: micromamba activate isoform_env && Rscript environments/r_scripts/install_isoform_r.R

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

bioc_pkgs <- c(
  "tximport",         # Import Salmon/Kallisto transcript-level quantification
  "tximeta",          # tximport with provenance tracking (optional but recommended)
  "DRIMSeq",          # Differential Transcript Usage via Dirichlet-multinomial model
  "GenomicFeatures",  # Build TxDb from GTF for transcript → gene mapping
  "BiocParallel",     # Parallel back-ends (MulticoreParam) for DRIMSeq
  "AnnotationDbi"     # select() to extract tx→gene mappings from TxDb
)
BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE)

message("Isoform DTU Bioconductor packages installed successfully.")
