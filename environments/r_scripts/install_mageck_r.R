#!/usr/bin/env Rscript
# install_mageck_r.R
# Run after activating mageck_env to install R packages not in bioconda.
# Usage: conda activate mageck_env && Rscript r_scripts/install_mageck_r.R

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

# ggprism: CRAN only
if (!require("ggprism", quietly = TRUE))
  install.packages("ggprism", repos = "https://cloud.r-project.org")

# msigdbr: MSigDB gene sets in R (CRAN / Bioconductor)
if (!require("msigdbr", quietly = TRUE))
  install.packages("msigdbr", repos = "https://cloud.r-project.org")

# ReactomePA: Reactome pathway enrichment
if (!require("ReactomePA", quietly = TRUE))
  BiocManager::install("ReactomePA")

message("MAGeCK R extras installed successfully.")
