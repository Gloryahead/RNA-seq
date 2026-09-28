#!/usr/bin/env Rscript
# figures.R — Publication-ready figures (Part 4)
# Volcano (EnhancedVolcano, NGS101 style), DEG heatmap
# Input:  DEG results TSV, DGE/DESeq2 RDS object
# Usage:  Called by Snakemake rule publication_figures

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(EnhancedVolcano)
  library(SummarizedExperiment)
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

degs    <- fread(opt$degs, data.table=FALSE)
lfc_col <- intersect(c("log2FC","log2FoldChange","logFC"), colnames(degs))[1]
p_col   <- intersect(c("FDR","adj.P.Val","padj"), colnames(degs))[1]

# ── EnhancedVolcano (NGS101 style) ────────────────────────────────────
n_up   <- sum(degs[[lfc_col]] >  opt$lfc & degs[[p_col]] < opt$alpha, na.rm=TRUE)
n_down <- sum(degs[[lfc_col]] < -opt$lfc & degs[[p_col]] < opt$alpha, na.rm=TRUE)

keyvals <- ifelse(
  degs[[lfc_col]] < -opt$lfc & degs[[p_col]] < opt$alpha, "blue",
  ifelse(degs[[lfc_col]] >  opt$lfc & degs[[p_col]] < opt$alpha, "red", "grey30")
)
names(keyvals)[keyvals == "blue"]   <- paste0("Downregulated (n=", n_down, ")")
names(keyvals)[keyvals == "red"]    <- paste0("Upregulated (n=", n_up, ")")
names(keyvals)[keyvals == "grey30"] <- "Non-significant"

pdf(opt$out_volcano, width=16/2.54, height=18/2.54)
print(EnhancedVolcano(degs,
  lab              = NA,
  x                = lfc_col,
  y                = p_col,
  title            = "",
  subtitle         = "",
  pCutoff          = opt$alpha,
  FCcutoff         = opt$lfc,
  gridlines.major  = FALSE,
  gridlines.minor  = FALSE,
  ylim             = c(0, max(-log10(degs[[p_col]]), na.rm=TRUE) + 0.5),
  colCustom        = keyvals,
  axisLabSize      = 20,
  labSize          = 3,
  legendLabSize    = 16,
  legendIconSize   = 7,
  captionLabSize   = 16,
  colAlpha         = 1,
  pointSize        = 0.3
))
dev.off()

# ── Heatmap (NGS101 style) ────────────────────────────────────────────
obj <- readRDS(opt$rds)
if (inherits(obj, "list") && "vsd" %in% names(obj)) {
  norm_mat  <- assay(obj$vsd)
  col_group <- obj$dds$group
  sample_names <- colnames(norm_mat)
} else if (inherits(obj, "DESeqDataSet")) {
  norm_mat  <- assay(obj)
  col_group <- obj$group
  sample_names <- colnames(norm_mat)
} else if (inherits(obj, "DGEList")) {
  norm_mat  <- edgeR::cpm(obj, log=TRUE)
  col_group <- obj$samples$group
  sample_names <- colnames(norm_mat)
} else {
  # limma-voom list
  norm_mat  <- obj$v$E
  col_group <- obj$v$targets$group
  sample_names <- colnames(norm_mat)
}

# All significant DEGs (not capped at 50)
sig_genes <- degs$gene[!is.na(degs[[p_col]]) &
                        degs[[p_col]] < opt$alpha &
                        abs(degs[[lfc_col]]) >= opt$lfc]
mat <- norm_mat[rownames(norm_mat) %in% sig_genes, , drop=FALSE]

# Rename rows to gene symbols if available
if ("gene_symbol" %in% colnames(degs)) {
  sym_map <- setNames(degs$gene_symbol, degs$gene)
  rownames(mat) <- ifelse(is.na(sym_map[rownames(mat)]), rownames(mat), sym_map[rownames(mat)])
}

# Order columns by group (control first)
if (!is.null(col_group)) {
  col_order <- order(col_group)
  mat        <- mat[, col_order, drop=FALSE]
  col_group  <- col_group[col_order]
  col_annot  <- data.frame(Treatment=as.character(col_group),
                            row.names=colnames(mat))
} else {
  col_annot <- NULL
}

if (nrow(mat) >= 2) {
  pdf(opt$out_heatmap, width=4, height=4)
  pheatmap(mat,
           scale           = "row",
           cluster_rows    = TRUE,
           cluster_cols    = FALSE,
           show_rownames   = FALSE,
           show_colnames   = TRUE,
           color           = colorRampPalette(c("navy","white","red"))(100),
           annotation_col  = col_annot,
           main            = "Differential Expression Heatmap")
  dev.off()
} else {
  message("Fewer than 2 significant DEGs at lfc >= ", opt$lfc, "; skipping heatmap.")
  pdf(opt$out_heatmap); plot.new(); dev.off()
}

message("Figures written.")
