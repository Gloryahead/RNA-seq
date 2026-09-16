#!/usr/bin/env Rscript
# install_cancer_subtype_r.R
# Run after activating cancer_subtype_env to install packages not in bioconda.
# Usage: conda activate cancer_subtype_env && Rscript r_scripts/install_cancer_subtype_r.R

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

# rmeta: meta-analysis plots (CRAN only, not in conda-forge)
if (!require("rmeta", quietly = TRUE))
  install.packages("rmeta", repos = "https://cloud.r-project.org")

# PAM50 classifier — GitHub only (no CRAN/Bioconductor release)
if (!require("devtools", quietly = TRUE))
  install.packages("devtools", repos = "https://cloud.r-project.org")

if (!require("PAM50", quietly = TRUE))
  devtools::install_github("ccchang0111/PAM50")

message("Cancer subtype R extras installed successfully.")
message("Example data: GSE209998 (retrieved via GEOquery in R)")
