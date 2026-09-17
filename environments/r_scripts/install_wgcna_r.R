#!/usr/bin/env Rscript
# install_wgcna_r.R
# Run after activating wgcna_env to install WGCNA and Bioconductor packages.
# r-wgcna is NOT in the conda YAML (Bioconductor transitive deps conflict on UA HPC).
# Usage: micromamba activate wgcna_env && Rscript environments/r_scripts/install_wgcna_r.R

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

bioc_pkgs <- c(
  "WGCNA",             # Gene co-expression network analysis (primary package)
  "impute",            # WGCNA dependency: missing-value imputation
  "preprocessCore",    # WGCNA dependency: quantile normalization
  "GO.db",             # Gene Ontology database
  "limma",             # Linear models (input preprocessing)
  "DESeq2",            # Differential expression (input preparation)
  "clusterProfiler",   # GO/KEGG pathway enrichment
  "org.Hs.eg.db",      # Human gene annotation
  "STRINGdb"           # PPI network from STRING database
)
BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE)

message("WGCNA and Bioconductor packages installed successfully.")
message("STRING database: species 9606 (Homo sapiens), version 11.5")
