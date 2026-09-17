#!/usr/bin/env Rscript
# install_genie3_masterreg_r.R
# Run after activating genie3_masterreg_env to install ALL Bioconductor packages.
# Usage: conda activate genie3_masterreg_env && Rscript environments/r_scripts/install_genie3_masterreg_r.R

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

bioc_pkgs <- c(
  "GENIE3",            # Gene regulatory network inference
  "dorothea",          # TF regulon database (DoRothEA)
  "RTN",               # Regulatory network and enrichment
  "DESeq2",            # Differential expression (input preparation)
  "edgeR",             # Differential expression
  "limma",             # Linear models
  "clusterProfiler",   # GO/KEGG pathway enrichment
  "org.Hs.eg.db",      # Human gene annotation
  "RegEnrich",         # Master regulator scoring
  "RedeR"              # Interactive network visualization
)
BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE)

message("GENIE3 / Master Regulator Bioconductor packages installed successfully.")
message("Example datasets:")
message("  GENIE3:     GSE261875 (TDP-43, human motor neurons)")
message("  Master Reg: GSE53239  (6 treatment + 6 control samples)")
