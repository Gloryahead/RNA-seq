#!/usr/bin/env Rscript
# normalization.R — library size / batch correction diagnostics
# Called by normalization_report rule in 02_deg_r.smk
# Outputs: per-method PCA/box plots; RLE plots; corrected matrix (optional)

suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(edgeR)
  library(sva)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(limma)
})

option_list <- list(
  make_option("--counts",  type="character", help="Raw counts matrix (featureCounts TSV)"),
  make_option("--samples", type="character", help="Sample metadata TSV (columns: sample,group,batch)"),
  make_option("--outdir",  type="character", default="results/normalization"),
  make_option("--methods", type="character", default="tmm,deseq2vst,rpkm",
              help="Comma-separated: tmm,deseq2vst,rpkm,rle"),
  make_option("--batch_correct", type="logical", default=FALSE,
              help="Run ComBat batch correction (requires batch column in samples)")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

# ── Load data ─────────────────────────────────────────────────────────
meta <- read.table(opt$samples, sep="\t", header=TRUE, stringsAsFactors=FALSE)
rownames(meta) <- meta$sample

raw <- read.table(opt$counts, sep="\t", header=TRUE, skip=1, check.names=FALSE)
gene_ids <- raw[[1]]
count_cols <- grep("\\.bam$", colnames(raw), value=TRUE)
mat <- as.matrix(raw[, count_cols])
colnames(mat) <- sub(".*/(.*)\\.bam$", "\\1", colnames(mat))
rownames(mat) <- gene_ids
mat <- mat[, meta$sample]  # align columns to sample order

# ── Helper: PCA plot ──────────────────────────────────────────────────
pca_plot <- function(norm_mat, meta, title, outfile) {
  pca <- prcomp(t(norm_mat), scale.=FALSE)
  var_pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  df <- data.frame(
    PC1 = pca$x[,1], PC2 = pca$x[,2],
    group = meta$group,
    sample = rownames(meta)
  )
  p <- ggplot(df, aes(PC1, PC2, color=group, label=sample)) +
    geom_point(size=3) +
    geom_text_repel(size=3, max.overlaps=20) +
    labs(
      title = title,
      x = paste0("PC1 (", var_pct[1], "%)"),
      y = paste0("PC2 (", var_pct[2], "%)")
    ) +
    theme_bw(base_size=13)
  ggsave(outfile, p, width=7, height=5)
  invisible(p)
}

# ── Helper: RLE plot ──────────────────────────────────────────────────
rle_plot <- function(norm_mat, meta, title, outfile) {
  rle <- norm_mat - apply(norm_mat, 1, median)
  df  <- reshape2::melt(rle)
  colnames(df) <- c("gene","sample","rle")
  df$group <- meta[as.character(df$sample), "group"]
  p <- ggplot(df, aes(sample, rle, fill=group)) +
    geom_boxplot(outlier.size=0.3, lwd=0.3) +
    geom_hline(yintercept=0, linetype="dashed", color="red") +
    labs(title=title, x=NULL, y="RLE") +
    theme_bw(base_size=11) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=7))
  ggsave(outfile, p, width=max(8, ncol(norm_mat)*0.4), height=5)
  invisible(p)
}

methods <- strsplit(opt$methods, ",")[[1]]

# ── TMM normalization (edgeR) ─────────────────────────────────────────
if ("tmm" %in% methods) {
  dge   <- DGEList(counts=mat, group=meta$group)
  dge   <- calcNormFactors(dge, method="TMM")
  tmm   <- cpm(dge, log=TRUE, prior.count=1)
  write.table(tmm, file.path(opt$outdir, "tmm_log2cpm.tsv"),
              sep="\t", quote=FALSE)
  pca_plot(tmm, meta, "PCA — TMM log2CPM",
           file.path(opt$outdir, "pca_tmm.pdf"))
  rle_plot(tmm, meta, "RLE — TMM log2CPM",
           file.path(opt$outdir, "rle_tmm.pdf"))
}

# ── DESeq2 VST ────────────────────────────────────────────────────────
if ("deseq2vst" %in% methods) {
  dds <- DESeqDataSetFromMatrix(mat, colData=meta, design=~group)
  dds <- estimateSizeFactors(dds)
  vst <- assay(varianceStabilizingTransformation(dds, blind=TRUE))
  write.table(vst, file.path(opt$outdir, "vst.tsv"),
              sep="\t", quote=FALSE)
  pca_plot(vst, meta, "PCA — DESeq2 VST",
           file.path(opt$outdir, "pca_vst.pdf"))
  rle_plot(vst, meta, "RLE — DESeq2 VST",
           file.path(opt$outdir, "rle_vst.pdf"))
}

# ── ComBat batch correction (on VST, optional) ────────────────────────
if (opt$batch_correct && "deseq2vst" %in% methods && "batch" %in% colnames(meta)) {
  vst <- read.table(file.path(opt$outdir, "vst.tsv"), sep="\t", header=TRUE)
  corrected <- ComBat(dat=as.matrix(vst), batch=meta$batch)
  write.table(corrected, file.path(opt$outdir, "vst_batch_corrected.tsv"),
              sep="\t", quote=FALSE)
  pca_plot(corrected, meta, "PCA — VST + ComBat",
           file.path(opt$outdir, "pca_vst_batch_corrected.pdf"))
}

message("normalization.R complete. Outputs in: ", opt$outdir)
