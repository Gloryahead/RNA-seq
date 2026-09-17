#!/usr/bin/env Rscript
# install_mageck_r.R
# Run after activating mageck_env to install Bioconductor + CRAN extras.
# Usage: conda activate mageck_env && Rscript environments/r_scripts/install_mageck_r.R

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

bioc_pkgs <- c(
  "clusterProfiler",    # GO/KEGG enrichment
  "org.Hs.eg.db",       # Human gene annotation
  "pathview",           # KEGG pathway visualization
  "enrichplot",         # Enrichment visualization
  "DOSE",               # Disease ontology enrichment
  "ReactomePA"          # Reactome pathway enrichment
)
BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE)

# ggprism: CRAN only
if (!require("ggprism", quietly = TRUE))
  install.packages("ggprism", repos = "https://cloud.r-project.org")

# msigdbr: MSigDB gene sets
if (!require("msigdbr", quietly = TRUE))
  install.packages("msigdbr", repos = "https://cloud.r-project.org")

message("MAGeCK R packages installed successfully.")
