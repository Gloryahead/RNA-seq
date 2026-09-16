#!/usr/bin/env Rscript
# install_genie3_masterreg_r.R
# Run after activating genie3_masterreg_env to install Bioconductor packages
# not available in bioconda (RegEnrich, RedeR).
# Usage: conda activate genie3_masterreg_env && Rscript r_scripts/install_genie3_masterreg_r.R

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

# RegEnrich: master regulator scoring (Bioconductor, not in bioconda)
if (!require("RegEnrich", quietly = TRUE))
  BiocManager::install("RegEnrich")

# RedeR: interactive network visualization (Bioconductor, not in bioconda)
if (!require("RedeR", quietly = TRUE))
  BiocManager::install("RedeR")

message("GENIE3 / Master Regulator R extras installed successfully.")
message("Example datasets:")
message("  GENIE3:        GSE261875 (TDP-43, human motor neurons)")
message("  Master Reg:    GSE53239  (6 treatment + 6 control samples)")
