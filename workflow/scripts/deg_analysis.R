#!/usr/bin/env Rscript
# deg_analysis.R — Differential expression analysis (Parts 3, 20)
# Methods: DESeq2, edgeR, or limma-voom (set via --method flag)
# Input:   count matrix (featureCounts format), sample sheet (TSV)
# Output:  DEG results TSV, PCA PDF, MA PDF, sample-distance heatmap PDF, RDS object
# Usage:   Called by Snakemake rule 02_deg_r; run standalone with:
#   Rscript workflow/scripts/deg_analysis.R \
#     --counts results/counts/counts_raw.tsv \
#     --samples config/samples.tsv \
#     --treatment KO --control WT \
#     --method DESeq2 --out_degs results/deg/KO_vs_WT_DEG_results.tsv [...]

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

# ── CLI arguments ─────────────────────────────────────────────────────
opt_list <- list(
  make_option("--counts",      type="character", help="featureCounts output TSV"),
  make_option("--samples",     type="character", help="Sample sheet TSV (sample, group, ...)"),
  make_option("--treatment",   type="character", help="Treatment group name"),
  make_option("--control",     type="character", help="Control group name"),
  make_option("--method",      type="character", default="DESeq2",
              help="DESeq2 | edgeR | limma-voom [default %default]"),
  make_option("--min_count",   type="integer",   default=10,
              help="CPM filter: min count in at least --min_samples samples [default %default]"),
  make_option("--min_samples", type="integer",   default=2,
              help="Minimum samples passing CPM filter [default %default]"),
  make_option("--lfc",         type="double",    default=0,
              help="log2FC threshold for DESeq2::results() lfcThreshold [default %default]"),
  make_option("--alpha",       type="double",    default=0.05,
              help="FDR cutoff [default %default]"),
  make_option("--organism",    type="character", default="mouse",
              help="human | mouse | rat [default %default]"),
  make_option("--out_degs",    type="character", help="Output DEG TSV"),
  make_option("--out_pca",     type="character", help="Output PCA PDF"),
  make_option("--out_ma",      type="character", help="Output MA plot PDF"),
  make_option("--out_heatmap", type="character", help="Output sample-distance heatmap PDF"),
  make_option("--out_rds",     type="character", help="Output RDS (DGE/DESeq2 object)")
)
opt <- parse_args(OptionParser(option_list=opt_list))

for (dir in unique(dirname(c(opt$out_degs, opt$out_pca, opt$out_ma, opt$out_heatmap, opt$out_rds))))
  dir.create(dir, showWarnings=FALSE, recursive=TRUE)

message("── DEG analysis: ", opt$treatment, " vs ", opt$control, " (", opt$method, ") ──")

# ── Load and format count matrix ──────────────────────────────────────
# featureCounts header: Geneid, Chr, Start, End, Strand, Length, sample1.bam, ...
raw <- fread(opt$counts, skip=1, data.table=FALSE)
count_cols <- grep("\\.bam$", colnames(raw), value=TRUE)
counts <- as.matrix(raw[, count_cols])
rownames(counts) <- raw$Geneid
colnames(counts) <- gsub(".*/|_Aligned.*|\\.bam$", "", colnames(counts))

# ── Load sample sheet and subset to this comparison ───────────────────
meta <- read.delim(opt$samples, stringsAsFactors=FALSE, comment.char="#")
meta <- meta[meta$group %in% c(opt$treatment, opt$control), ]
meta$group <- factor(meta$group, levels=c(opt$control, opt$treatment))
meta <- meta[order(meta$group), ]
counts <- counts[, meta$sample]

message("Samples: ", paste(colnames(counts), collapse=", "))
message("Dimensions before filter: ", nrow(counts), " genes x ", ncol(counts), " samples")

# ── CPM filter ────────────────────────────────────────────────────────
suppressPackageStartupMessages(library(edgeR))
dge_raw <- DGEList(counts=counts, group=meta$group)
keep     <- rowSums(cpm(dge_raw) >= opt$min_count) >= opt$min_samples
counts_f <- counts[keep, ]
message("After CPM filter: ", nrow(counts_f), " genes retained")

# ── Run chosen DEG method ─────────────────────────────────────────────
if (opt$method == "DESeq2") {
  suppressPackageStartupMessages(library(DESeq2))
  dds <- DESeqDataSetFromMatrix(
    countData = counts_f,
    colData   = meta,
    design    = ~ group
  )
  dds    <- DESeq(dds, quiet=TRUE)
  res    <- results(dds, contrast=c("group", opt$treatment, opt$control),
                    alpha=opt$alpha, lfcThreshold=opt$lfc)
  res    <- lfcShrink(dds, contrast=c("group", opt$treatment, opt$control),
                      type="ashr", res=res, quiet=TRUE)
  res_df <- as.data.frame(res)
  res_df <- data.frame(gene=rownames(res_df), res_df, stringsAsFactors=FALSE)
  setnames_map <- c(log2FoldChange="log2FC", pvalue="pvalue", padj="FDR")
  colnames(res_df) <- ifelse(colnames(res_df) %in% names(setnames_map),
                             setnames_map[colnames(res_df)], colnames(res_df))
  vsd    <- vst(dds, blind=FALSE)
  norm_counts <- assay(vsd)
  dge_obj <- dds
  saveRDS(list(dds=dds, vsd=vsd, res=res), opt$out_rds)

} else if (opt$method == "edgeR") {
  suppressPackageStartupMessages(library(edgeR))
  dge   <- DGEList(counts=counts_f, group=meta$group)
  dge   <- calcNormFactors(dge, method="TMM")
  design <- model.matrix(~ group, data=meta)
  dge   <- estimateDisp(dge, design)
  fit   <- glmQLFit(dge, design)
  qlf   <- glmQLFTest(fit, coef=ncol(design))
  res_df <- topTags(qlf, n=Inf, adjust.method="BH")$table
  res_df <- data.frame(gene=rownames(res_df), res_df, stringsAsFactors=FALSE)
  norm_counts <- cpm(dge, log=TRUE)
  saveRDS(dge, opt$out_rds)

} else if (opt$method == "limma-voom") {
  suppressPackageStartupMessages({library(edgeR); library(limma)})
  dge    <- DGEList(counts=counts_f, group=meta$group)
  dge    <- calcNormFactors(dge, method="TMM")
  design <- model.matrix(~ group, data=meta)
  v      <- voom(dge, design, plot=FALSE)
  fit    <- lmFit(v, design)
  fit    <- eBayes(fit)
  res_df <- topTable(fit, coef=ncol(design), n=Inf, adjust.method="BH")
  res_df <- data.frame(gene=rownames(res_df), res_df, stringsAsFactors=FALSE)
  norm_counts <- v$E
  saveRDS(list(dge=dge, v=v, fit=fit), opt$out_rds)

} else {
  stop("Unknown method: ", opt$method, ". Use DESeq2, edgeR, or limma-voom.")
}

# ── Write DEG table ───────────────────────────────────────────────────
res_df <- res_df[order(res_df$FDR, na.last=TRUE), ]
fwrite(res_df, opt$out_degs, sep="\t", na="NA")
n_sig <- sum(!is.na(res_df$FDR) & res_df$FDR < opt$alpha, na.rm=TRUE)
message("Significant DEGs (FDR < ", opt$alpha, "): ", n_sig)

# ── PCA plot ──────────────────────────────────────────────────────────
pca     <- prcomp(t(norm_counts), scale.=FALSE)
pca_df  <- as.data.frame(pca$x[, 1:2])
pca_df$sample <- rownames(pca_df)
pca_df$group  <- meta$group[match(pca_df$sample, meta$sample)]
var_pct <- round(100 * summary(pca)$importance[2, 1:2], 1)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color=group, label=sample)) +
  geom_point(size=4, alpha=.9) +
  ggrepel::geom_text_repel(size=3, max.overlaps=20) +
  labs(x=paste0("PC1 (", var_pct[1], "%)"),
       y=paste0("PC2 (", var_pct[2], "%)"),
       title=paste(opt$method, "PCA:", opt$treatment, "vs", opt$control),
       color="Group") +
  theme_classic(base_size=12) +
  theme(legend.position="right")
ggsave(opt$out_pca, p_pca, width=7, height=5)

# ── MA plot ───────────────────────────────────────────────────────────
ma_df <- res_df
ma_col <- ifelse(opt$method == "DESeq2", "baseMean",
          ifelse(opt$method == "edgeR",  "logCPM", "AveExpr"))
if (!ma_col %in% colnames(ma_df)) ma_col <- colnames(ma_df)[2]
ma_df$sig <- !is.na(ma_df$FDR) & ma_df$FDR < opt$alpha

pdf(opt$out_ma, width=7, height=5)
plot(log2(ma_df[[ma_col]]+1), ma_df$log2FC,
     pch=20, cex=.4, col=ifelse(ma_df$sig, "firebrick", "grey60"),
     xlab="log2(Mean expression)", ylab="log2 Fold Change",
     main=paste("MA plot:", opt$treatment, "vs", opt$control))
abline(h=0, lty=2, col="navy")
legend("topright", legend=c(paste0("FDR < ", opt$alpha), "n.s."),
       col=c("firebrick","grey60"), pch=20)
dev.off()

# ── Sample-distance heatmap ───────────────────────────────────────────
dist_mat <- as.matrix(dist(t(norm_counts)))
ann      <- data.frame(Group=meta$group, row.names=meta$sample)
pdf(opt$out_heatmap, width=6, height=5)
pheatmap(dist_mat,
         annotation_col=ann,
         clustering_distance_rows="euclidean",
         clustering_distance_cols="euclidean",
         color=colorRampPalette(rev(brewer.pal(9,"Blues")))(100),
         main="Sample-to-sample Euclidean distances (VST/log-CPM)")
dev.off()

message("Done. DEG results: ", opt$out_degs)
