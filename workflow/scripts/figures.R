#!/usr/bin/env Rscript
# figures.R — Publication-ready figures (Part 4)
# Volcano (EnhancedVolcano), top-50 DEG heatmap
# Input:  DEG results TSV, DGE/DESeq2 RDS object
# Usage:  Called by Snakemake rule publication_figures

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(EnhancedVolcano)
})

opt_list <- list(
  make_option("--degs",        type="character"),
  make_option("--rds",         type="character"),
  make_option("--treatment",   type="character"),
  make_option("--control",     type="character"),
  make_option("--alpha",       type="double",  default=0.05),
  make_option("--lfc",         type="double",  default=1.0),
  make_option("--out_volcano", type="character"),
  make_option("--out_heatmap", type="character")
)
opt <- parse_args(OptionParser(option_list=opt_list))

degs <- fread(opt$degs, data.table=FALSE)
lfc_col <- intersect(c("log2FC","log2FoldChange","logFC"), colnames(degs))[1]

# ── EnhancedVolcano ───────────────────────────────────────────────────
pdf(opt$out_volcano, width=9, height=7)
EnhancedVolcano(degs,
  lab           = degs$gene,
  x             = lfc_col,
  y             = "FDR",
  pCutoff       = opt$alpha,
  FCcutoff      = opt$lfc,
  title         = paste(opt$treatment, "vs", opt$control),
  subtitle      = paste0("FDR < ", opt$alpha, "  |log2FC| > ", opt$lfc),
  legendPosition= "right",
  pointSize     = 2,
  labSize       = 3.5,
  col           = c("grey70","grey70","steelblue","firebrick"))
dev.off()

# ── Top-50 DEG heatmap ────────────────────────────────────────────────
obj <- readRDS(opt$rds)
if (inherits(obj, "list") && "vsd" %in% names(obj)) {
  norm_mat <- assay(obj$vsd)
} else if (inherits(obj, "DESeqDataSet")) {
  norm_mat <- assay(SummarizedExperiment::assay(obj))
} else if (inherits(obj, "DGEList")) {
  norm_mat <- edgeR::cpm(obj, log=TRUE)
} else {
  # limma-voom list
  norm_mat <- obj$v$E
}

sig_degs <- degs[!is.na(degs$FDR) & degs$FDR < opt$alpha & abs(degs[[lfc_col]]) >= opt$lfc, ]
top50    <- head(sig_degs$gene[order(sig_degs$FDR)], 50)
mat      <- norm_mat[rownames(norm_mat) %in% top50, , drop=FALSE]
mat      <- mat[match(intersect(top50, rownames(mat)), rownames(mat)), ]

if (nrow(mat) >= 2) {
  mat_scaled <- t(scale(t(mat)))
  pdf(opt$out_heatmap, width=9, height=max(6, nrow(mat)*0.22 + 2))
  pheatmap(mat_scaled,
           color           = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),
           breaks          = seq(-3, 3, length.out=101),
           cluster_rows    = TRUE,
           cluster_cols    = TRUE,
           show_rownames   = TRUE,
           fontsize_row    = 8,
           main            = paste0("Top ", nrow(mat), " DEGs — z-score"))
  dev.off()
} else {
  message("Fewer than 2 significant DEGs; skipping heatmap.")
  pdf(opt$out_heatmap); plot.new(); dev.off()
}

message("Figures written.")
