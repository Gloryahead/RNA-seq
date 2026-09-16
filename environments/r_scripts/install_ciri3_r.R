#!/usr/bin/env Rscript
# install_ciri3_r.R
# Run after activating ciri3_env to install rMATS for differential splicing.
# Usage: conda activate ciri3_env && Rscript r_scripts/install_ciri3_r.R

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

# rMATS is primarily a Python/command-line tool but has an R companion package.
# If the conda bioconda package is unavailable, install via pip in the shell:
#   pip install rmats
# Or use rMATS-turbo from: https://github.com/Xinglab/rmats-turbo

message("Note: rMATS is a command-line tool. Install via conda or pip separately:")
message("  conda install -c bioconda rmats")
message("  OR: pip install rmats")
message("  OR: https://github.com/Xinglab/rmats-turbo")
