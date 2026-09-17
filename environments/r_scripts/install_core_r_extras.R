#!/usr/bin/env Rscript
# install_core_r_extras.R
# Run after activating rnaseq_r_env to install Bioconductor + CRAN extras.
# Usage: micromamba activate rnaseq_r_env && Rscript environments/r_scripts/install_core_r_extras.R

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

bioc_pkgs <- c(
  "DESeq2",             # Differential expression
  "edgeR",              # Differential expression (TMM normalization)
  "limma",              # Linear models for microarray/RNA-seq
  "tximport",           # Import Salmon/Kallisto quantification
  "biomaRt",            # Bioconductor interface to Ensembl BioMart
  "clusterProfiler",    # GO/KEGG pathway enrichment
  "EnhancedVolcano",    # Publication-quality volcano plots
  "enrichplot",         # Visualization for enrichment results
  "rtracklayer",        # GTF/GFF import/export
  "sva"                 # Surrogate variable analysis (batch correction)
)
BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE)

# ggprism: prism-style publication figures (Part 4), CRAN only
if (!require("ggprism", quietly = TRUE))
  install.packages("ggprism", repos = "https://cloud.r-project.org")

message("Core R + Bioconductor packages installed successfully.")
message("MSigDB gene sets for GSEA (Part 5) must be downloaded manually from:")
message("  https://www.gsea-msigdb.org/gsea/msigdb")
message("  File used in tutorial: msigdb.v2024.1.Mm.symbols.gmt (mouse)")
