#!/usr/bin/env Rscript
# install_core_r_extras.R
# Run after activating rnaseq_r_env to install CRAN packages not in bioconda.
# Usage: conda activate rnaseq_r_env && Rscript r_scripts/install_core_r_extras.R

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

# ggprism: ggplot2 extension for prism-style publication figures (Part 4)
if (!require("ggprism", quietly = TRUE))
  install.packages("ggprism", repos = "https://cloud.r-project.org")

message("Core R extras installed successfully.")
message("MSigDB gene sets for GSEA (Part 5) must be downloaded manually from:")
message("  https://www.gsea-msigdb.org/gsea/msigdb")
message("  File used in tutorial: msigdb.v2024.1.Mm.symbols.gmt (mouse)")
