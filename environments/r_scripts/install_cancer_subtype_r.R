#!/usr/bin/env Rscript
# install_cancer_subtype_r.R
# Run after activating cancer_subtype_env.
# Usage: conda activate cancer_subtype_env && Rscript environments/r_scripts/install_cancer_subtype_r.R

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

bioc_pkgs <- c(
  "edgeR",             # Differential expression (shared with rnaseq_r_env)
  "limma",             # Linear models for RNA-seq
  "GEOquery",          # GEO data retrieval (GSE209998)
  "GSVA",              # Gene set variation analysis
  "genefu",            # PAM50 and breast cancer gene signatures
  "EnsDb.Hsapiens.v86" # Ensembl annotation (hg38 / GRCh38)
)
BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE)

# rmeta: meta-analysis plots (CRAN only)
if (!require("rmeta", quietly = TRUE))
  install.packages("rmeta", repos = "https://cloud.r-project.org")

# PAM50 classifier — GitHub only (no CRAN/Bioconductor release)
if (!require("PAM50", quietly = TRUE))
  devtools::install_github("ccchang0111/PAM50")

message("Cancer subtype R packages installed successfully.")
message("Example data: GSE209998 (retrieved via GEOquery in R)")
