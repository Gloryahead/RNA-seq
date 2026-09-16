#!/usr/bin/env Rscript
# genie3.R — Gene Regulatory Network inference with GENIE3 + DoRothEA TF enrichment
# Tutorial ref: GENIE3 tutorial (ngs101.com)

suppressPackageStartupMessages({
  library(optparse)
  library(GENIE3)
  library(DESeq2)
  library(dorothea)
  library(dplyr)
  library(ggplot2)
  library(igraph)
  library(reshape2)
})

option_list <- list(
  make_option("--counts",        type="character", help="Raw counts matrix (featureCounts TSV)"),
  make_option("--samples",       type="character", help="Sample metadata TSV"),
  make_option("--outdir",        type="character", default="results/networks/genie3"),
  make_option("--organism",      type="character", default="human"),
  make_option("--top_links",     type="integer",   default=5000,
              help="Number of top regulatory links to keep"),
  make_option("--tf_confidence", type="character", default="A,B",
              help="DoRothEA confidence levels (A,B,C,D)")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

# ── Load and normalize ────────────────────────────────────────────────
meta <- read.table(opt$samples, sep="\t", header=TRUE, stringsAsFactors=FALSE)
rownames(meta) <- meta$sample

raw <- read.table(opt$counts, sep="\t", header=TRUE, skip=1, check.names=FALSE)
count_cols <- grep("\\.bam$", colnames(raw), value=TRUE)
mat <- as.matrix(raw[, count_cols])
colnames(mat) <- sub(".*/(.*)\\.bam$", "\\1", colnames(mat))
rownames(mat) <- raw[[1]]
mat <- mat[, meta$sample]

dds     <- DESeqDataSetFromMatrix(mat, colData=meta, design=~group)
dds     <- estimateSizeFactors(dds)
exprMat <- assay(varianceStabilizingTransformation(dds, blind=TRUE))

# ── Identify TF list from DoRothEA ───────────────────────────────────
conf_levels <- strsplit(opt$tf_confidence, ",")[[1]]
dorothea_df <- if (opt$organism == "mouse") {
  dorothea_mm
} else {
  dorothea_hs
}
regulons <- dorothea_df %>%
  filter(confidence %in% conf_levels) %>%
  pull(tf) %>%
  unique()
tf_in_data <- intersect(regulons, rownames(exprMat))
message(length(tf_in_data), " TFs from DoRothEA found in expression matrix")

# ── GENIE3 regulatory network inference ──────────────────────────────
message("Running GENIE3 (this may take several minutes)...")
weight_mat <- GENIE3(exprMat, regulators=tf_in_data, nCores=1)

# Extract top links
links <- getLinkList(weight_mat, reportMax=opt$top_links)
colnames(links) <- c("regulator", "target", "importance")
write.table(links, file.path(opt$outdir, "regulatory_links.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)
message("Saved top ", nrow(links), " regulatory links")

# ── Network visualization with igraph ────────────────────────────────
# Keep only the top 200 links for plotting
top_plot <- head(links, 200)
g <- graph_from_data_frame(top_plot, directed=TRUE)
V(g)$is_tf <- V(g)$name %in% tf_in_data

pdf(file.path(opt$outdir, "network_plot.pdf"), width=12, height=10)
plot(
  g,
  vertex.size    = ifelse(V(g)$is_tf, 8, 4),
  vertex.color   = ifelse(V(g)$is_tf, "#E63946", "#457B9D"),
  vertex.label   = ifelse(degree(g, mode="out") > 5, V(g)$name, NA),
  vertex.label.cex = 0.6,
  edge.arrow.size  = 0.2,
  edge.width       = top_plot$importance * 5,
  layout           = layout_with_fr(g),
  main             = "GENIE3 Gene Regulatory Network (top 200 links)"
)
legend("bottomleft", legend=c("TF","Target"),
       pch=21, pt.bg=c("#E63946","#457B9D"), bty="n")
dev.off()

# ── TF hub summary: regulators sorted by out-degree ──────────────────
hub_summary <- links %>%
  group_by(regulator) %>%
  summarise(n_targets=n(), mean_importance=mean(importance)) %>%
  arrange(desc(n_targets))
write.table(hub_summary, file.path(opt$outdir, "tf_hub_summary.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

message("genie3.R complete. Outputs in: ", opt$outdir)
