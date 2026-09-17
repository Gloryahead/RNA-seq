#!/usr/bin/env Rscript
# install_wgcna_r.R
# Run after activating wgcna_env to install Bioconductor packages.
# Usage: conda activate wgcna_env && Rscript environments/r_scripts/install_wgcna_r.R

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

bioc_pkgs <- c(
  "limma",             # Linear models (input preprocessing)
  "DESeq2",            # Differential expression (input preparation)
  "clusterProfiler",   # GO/KEGG pathway enrichment
  "org.Hs.eg.db",      # Human gene annotation
  "GO.db",             # Gene Ontology database
  "STRINGdb"           # PPI network from STRING database
)
BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE)

message("WGCNA Bioconductor packages installed successfully.")
message("STRING database: species 9606 (Homo sapiens), version 11.5")
