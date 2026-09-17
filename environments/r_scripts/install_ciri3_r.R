#!/usr/bin/env Rscript
# install_ciri3_r.R
# Run after activating ciri3_env to install Bioconductor packages for DE.
# Usage: micromamba activate ciri3_env && Rscript environments/r_scripts/install_ciri3_r.R

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

BiocManager::install(c("edgeR", "limma"), ask = FALSE, update = FALSE)

message("CIRI3 R packages installed successfully.")
message("Note: rMATS is a command-line tool. Install via conda or pip separately:")
message("  conda install -c bioconda rmats")
message("  OR: pip install rmats")
message("  OR: https://github.com/Xinglab/rmats-turbo")
