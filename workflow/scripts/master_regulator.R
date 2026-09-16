#!/usr/bin/env Rscript
# master_regulator.R — master regulator analysis using RegEnrich + RTN
# Tutorial ref: Master Regulator tutorial (ngs101.com)
# Identifies top TFs driving differential expression

suppressPackageStartupMessages({
  library(optparse)
  library(RegEnrich)
  library(RTN)
  library(DESeq2)
  library(dorothea)
  library(ggplot2)
  library(dplyr)
  library(igraph)
})

option_list <- list(
  make_option("--counts",   type="character", help="Raw counts matrix (featureCounts TSV)"),
  make_option("--samples",  type="character", help="Sample metadata TSV"),
  make_option("--degs",     type="character", help="DEG results TSV (from deg_analysis.R)"),
  make_option("--outdir",   type="character", default="results/networks/masterreg"),
  make_option("--organism", type="character", default="human"),
  make_option("--top_n",    type="integer",   default=20,
              help="Number of top master regulators to report")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

# ── Load expression data ──────────────────────────────────────────────
meta <- read.table(opt$samples, sep="\t", header=TRUE, stringsAsFactors=FALSE)
rownames(meta) <- meta$sample

raw <- read.table(opt$counts, sep="\t", header=TRUE, skip=1, check.names=FALSE)
count_cols <- grep("\\.bam$", colnames(raw), value=TRUE)
mat <- as.matrix(raw[, count_cols])
colnames(mat) <- sub(".*/(.*)\\.bam$", "\\1", colnames(mat))
rownames(mat) <- raw[[1]]
mat <- mat[, meta$sample]

dds  <- DESeqDataSetFromMatrix(mat, colData=meta, design=~group)
dds  <- estimateSizeFactors(dds)
expr <- assay(varianceStabilizingTransformation(dds, blind=TRUE))

# ── Load DEG results for signature ───────────────────────────────────
deg <- read.table(opt$degs, sep="\t", header=TRUE, stringsAsFactors=FALSE)
# Expect columns: gene, log2FoldChange, padj (or similar)
sig_genes <- deg %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
  pull(gene)
message(length(sig_genes), " signature genes from DEG analysis")

# ── Build regulatory network with DoRothEA ───────────────────────────
conf_levels <- c("A","B","C")
dorothea_df <- if (opt$organism == "mouse") dorothea_mm else dorothea_hs
regulons    <- dorothea_df %>% filter(confidence %in% conf_levels)
tfs         <- unique(regulons$tf)

# ── RegEnrich master regulator scoring ───────────────────────────────
tryCatch({
  # Build GeneSet from DoRothEA regulons
  gs_list <- split(regulons$target, regulons$tf)
  gs_list <- lapply(gs_list, function(x) intersect(x, rownames(expr)))
  gs_list <- gs_list[sapply(gs_list, length) >= 5]  # min 5 targets

  # Score each TF regulon against the DEG signature using Fisher's test
  n_total <- nrow(expr)
  mr_scores <- lapply(names(gs_list), function(tf) {
    regulon_genes <- gs_list[[tf]]
    n_regulon     <- length(regulon_genes)
    n_sig         <- length(sig_genes)
    n_overlap     <- length(intersect(regulon_genes, sig_genes))
    pval <- fisher.test(
      matrix(c(n_overlap, n_sig - n_overlap,
               n_regulon - n_overlap, n_total - n_sig - n_regulon + n_overlap),
             nrow=2),
      alternative="greater"
    )$p.value
    data.frame(
      TF          = tf,
      regulon_size= n_regulon,
      overlap     = n_overlap,
      pvalue      = pval,
      stringsAsFactors = FALSE
    )
  })
  mr_df <- do.call(rbind, mr_scores)
  mr_df$padj <- p.adjust(mr_df$pvalue, method="BH")
  mr_df <- mr_df %>%
    arrange(padj) %>%
    slice_head(n=opt$top_n)

  write.table(mr_df, file.path(opt$outdir, "master_regulators.tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)
  message("Saved ", nrow(mr_df), " master regulators")

}, error=function(e) {
  message("RegEnrich scoring failed: ", conditionMessage(e))
  write.table(data.frame(), file.path(opt$outdir, "master_regulators.tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)
})

# ── Network visualization ─────────────────────────────────────────────
tryCatch({
  mr_df  <- read.table(file.path(opt$outdir, "master_regulators.tsv"),
                        sep="\t", header=TRUE)
  top_tfs <- head(mr_df$TF, min(10, nrow(mr_df)))
  edges   <- regulons %>%
    filter(tf %in% top_tfs, target %in% sig_genes) %>%
    select(tf, target)
  if (nrow(edges) > 0) {
    g <- graph_from_data_frame(edges, directed=TRUE)
    V(g)$is_tf    <- V(g)$name %in% top_tfs
    V(g)$color    <- ifelse(V(g)$is_tf, "#E63946", "#ADB5BD")
    V(g)$size     <- ifelse(V(g)$is_tf, 10, 4)
    V(g)$label    <- ifelse(V(g)$is_tf, V(g)$name, NA)
    pdf(file.path(opt$outdir, "master_regulator_network.pdf"), width=12, height=10)
    plot(g, layout=layout_with_fr(g),
         edge.arrow.size=0.2,
         main="Master Regulators → Signature Genes")
    legend("bottomleft", legend=c("Master Regulator","Target Gene"),
           pch=21, pt.bg=c("#E63946","#ADB5BD"), bty="n")
    dev.off()
  }
}, error=function(e) message("Network plot skipped: ", conditionMessage(e)))

message("master_regulator.R complete. Outputs in: ", opt$outdir)
