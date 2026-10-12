#!/usr/bin/env Rscript
# isoform_analysis.R
# Differential Transcript Usage (DTU) analysis using tximport + DRIMSeq.
#
# Pipeline:
#   1. Build transcript-to-gene mapping from GTF (GenomicFeatures)
#   2. Import Salmon quant.sf files at transcript level (tximport)
#   3. Subset to treatment + control samples
#   4. Create DRIMSeq dmDSdata object → dmFilter → dmPrecision → dmFit → dmTest
#   5. Extract gene- and transcript-level results
#   6. Write TSV + proportion plots PDF
#
# Usage (called from Snakemake rule isoform_analysis):
#   Rscript isoform_analysis.R \
#       --quant_files "s1/quant.sf s2/quant.sf ..." \
#       --gtf /ref/annotation.gtf \
#       --samples config/samples.tsv \
#       --treatment KRAS_SPIB \
#       --control KRAS \
#       --organism mouse \
#       --min_feature_expr 10 \
#       --min_gene_expr 10 \
#       --min_feature_prop 0.1 \
#       --threads 8 \
#       --out_dtu results/isoform/comp_dtu_results.tsv \
#       --out_plots results/isoform/comp_dtu_plots.pdf

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(BiocParallel)
  library(GenomicFeatures)
  library(AnnotationDbi)
  library(tximport)
  library(DRIMSeq)
})

# ── CLI arguments ─────────────────────────────────────────────────────
opt_list <- list(
  make_option("--quant_files",      type = "character",
              help = "Space-separated paths to quant.sf files"),
  make_option("--gtf",              type = "character",
              help = "Path to GTF annotation file"),
  make_option("--samples",          type = "character",
              help = "Sample sheet TSV (columns: sample, group, ...)"),
  make_option("--treatment",        type = "character",
              help = "Treatment group name"),
  make_option("--control",          type = "character",
              help = "Control group name"),
  make_option("--min_feature_expr", type = "integer", default = 10L,
              help = "DRIMSeq min transcript count [default %default]"),
  make_option("--min_gene_expr",    type = "integer", default = 10L,
              help = "DRIMSeq min gene count [default %default]"),
  make_option("--min_feature_prop", type = "double",  default = 0.1,
              help = "DRIMSeq min transcript proportion [default %default]"),
  make_option("--organism",         type = "character", default = "mouse",
              help = "Organism: human | mouse | rat [default %default]"),
  make_option("--threads",          type = "integer", default = 4L,
              help = "Number of threads for BiocParallel [default %default]"),
  make_option("--out_dtu",          type = "character",
              help = "Output DTU results TSV"),
  make_option("--out_plots",        type = "character",
              help = "Output plots PDF")
)
opt <- parse_args(OptionParser(option_list = opt_list))

# ── Ensure output directories exist ──────────────────────────────────
dir.create(dirname(opt$out_dtu),   recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$out_plots), recursive = TRUE, showWarnings = FALSE)

# ── Configure parallelism ─────────────────────────────────────────────
bpparam <- MulticoreParam(workers = opt$threads)

# ── Load sample sheet ─────────────────────────────────────────────────
samples_df <- fread(opt$samples, sep = "\t", header = TRUE,
                    colClasses = "character")
samples_df <- samples_df[group %in% c(opt$treatment, opt$control)]
samples_df[, group := factor(group, levels = c(opt$control, opt$treatment))]
setkey(samples_df, sample)

# ── Parse quant file paths ────────────────────────────────────────────
quant_paths <- strsplit(trimws(opt$quant_files), "\\s+")[[1]]
# Derive sample names from parent directory names  (e.g., .../salmon/SAMPLE/quant.sf)
path_samples <- basename(dirname(quant_paths))
names(quant_paths) <- path_samples

# Subset to samples in this comparison, in the order they appear in samples_df
comp_samples <- samples_df$sample
missing <- setdiff(comp_samples, path_samples)
if (length(missing) > 0)
  stop("quant.sf paths missing for samples: ", paste(missing, collapse = ", "))
quant_paths <- quant_paths[comp_samples]

# ── Build transcript → gene map from GTF ─────────────────────────────
message("Building transcript-to-gene map from GTF...")
txdb    <- makeTxDbFromGFF(opt$gtf, format = "GTF")
tx2gene <- AnnotationDbi::select(txdb,
                                  keys    = keys(txdb, keytype = "TXNAME"),
                                  columns = c("TXNAME","GENEID"),
                                  keytype = "TXNAME")
tx2gene <- as.data.table(tx2gene)[!is.na(GENEID)]
setnames(tx2gene, c("TXNAME","GENEID"), c("tx_id","gene_id"))

# ── Import Salmon quant files at transcript level ─────────────────────
message("Importing Salmon quant files with tximport...")
txi <- tximport(
  quant_paths,
  type    = "salmon",
  tx2gene = tx2gene[, .(tx_id, gene_id)],
  txOut   = TRUE          # transcript-level counts
)

# Round scaled counts for DRIMSeq (expects integer-like input)
counts_mat <- round(txi$counts)
mode(counts_mat) <- "integer"

# ── Build DRIMSeq sample data ─────────────────────────────────────────
sample_data <- data.frame(
  sample_id = comp_samples,
  group     = samples_df[comp_samples, group],
  stringsAsFactors = FALSE
)

# ── Create DRIMSeq data object ────────────────────────────────────────
message("Building DRIMSeq data object...")
d <- dmDSdata(counts = counts_mat, samples = sample_data)

# ── Filter ────────────────────────────────────────────────────────────
message("Filtering low-count transcripts/genes...")
d <- dmFilter(d,
              min_samps_feature_expr = 1L,
              min_samps_feature_prop = 1L,
              min_samps_gene_expr    = 1L,
              min_feature_expr       = opt$min_feature_expr,
              min_feature_prop       = opt$min_feature_prop,
              min_gene_expr          = opt$min_gene_expr)

message("Genes after filter: ", length(d))

if (length(d) == 0) {
  message("No genes passed DRIMSeq filter. Writing empty results.")
  empty_dt <- data.table(
    gene_id      = character(),
    transcript_id = character(),
    lr           = numeric(),
    df           = integer(),
    pvalue       = numeric(),
    adj_pvalue   = numeric()
  )
  fwrite(empty_dt, opt$out_dtu, sep = "\t")
  pdf(opt$out_plots, width = 8, height = 6)
  plot.new(); text(0.5, 0.5, "No genes passed DRIMSeq filter.", cex = 1.2)
  dev.off()
  quit(save = "no", status = 0)
}

# ── Precision, Fit, Test ──────────────────────────────────────────────
message("Estimating precision...")
design_full <- model.matrix(~group, data = DRIMSeq::samples(d))
d <- dmPrecision(d, design = design_full, BPPARAM = bpparam)

message("Fitting model...")
d <- dmFit(d, design = design_full, BPPARAM = bpparam)

message("Testing DTU...")
d <- dmTest(d, coef = paste0("group", opt$treatment), BPPARAM = bpparam)

# ── Extract results ───────────────────────────────────────────────────
res_feature <- as.data.table(results(d, level = "feature"))
res_gene    <- as.data.table(results(d, level = "gene"))

# Merge gene-level adj_pvalue into feature table as gene_adj_pvalue
setnames(res_gene, "adj_pvalue", "gene_adj_pvalue")
res_merged <- merge(res_feature,
                    res_gene[, .(gene_id, gene_adj_pvalue)],
                    by = "gene_id", all.x = TRUE)

setorder(res_merged, adj_pvalue, na.last = TRUE)
fwrite(res_merged, opt$out_dtu, sep = "\t")
message("DTU results written: ", nrow(res_merged), " transcripts.")

# ── Significant DTU genes ─────────────────────────────────────────────
sig_genes <- res_gene[!is.na(gene_adj_pvalue) & gene_adj_pvalue < 0.05,
                      gene_id]
n_sig_genes        <- length(sig_genes)
n_sig_transcripts  <- nrow(res_merged[!is.na(adj_pvalue) & adj_pvalue < 0.05])
message("Significant DTU genes (FDR<0.05): ", n_sig_genes)
message("Significant DTU transcripts (FDR<0.05): ", n_sig_transcripts)

# ── Plots ─────────────────────────────────────────────────────────────
pdf(opt$out_plots, width = 10, height = 8)

# 1. Proportion plots for top 9 significant DTU genes
if (n_sig_genes > 0) {
  top_genes <- sig_genes[seq_len(min(9, n_sig_genes))]
  for (gid in top_genes) {
    tryCatch(
      plotProportions(d, gene_id = gid,
                      group_variable = "group",
                      plot_type = "ribbonplot"),
      error = function(e) {
        message("plotProportions failed for gene ", gid, ": ", e$message)
        plot.new()
        text(0.5, 0.5, paste("plotProportions failed for:", gid), cex = 0.9)
      }
    )
  }
} else {
  plot.new()
  text(0.5, 0.5, "No significant DTU genes at FDR < 0.05 for proportion plots.",
       cex = 1.0)
}

# 2. Volcano-style: log2(proportion change) vs -log10(pvalue)
if (nrow(res_merged) > 0 && "lr" %in% names(res_merged)) {
  plot_dt <- res_merged[!is.na(pvalue) & !is.na(lr)]
  plot_dt[, neg_log10p := -log10(pmax(pvalue, 1e-300))]
  plot_dt[, significant := !is.na(adj_pvalue) & adj_pvalue < 0.05]
  p_vol <- ggplot(plot_dt, aes(x = lr, y = neg_log10p, colour = significant)) +
    geom_point(alpha = 0.4, size = 1.2) +
    scale_colour_manual(values = c("FALSE" = "grey60", "TRUE" = "#d62728"),
                        labels = c("FALSE" = "FDR ≥ 0.05", "TRUE" = "FDR < 0.05"),
                        name = "") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "black") +
    labs(title = "DTU: likelihood-ratio statistic vs significance",
         x = "Likelihood-ratio (lr)", y = expression(-log[10](p-value))) +
    theme_bw(base_size = 12)
  print(p_vol)
} else {
  plot.new()
  text(0.5, 0.5, "Insufficient data for volcano plot.", cex = 1.0)
}

# 3. Summary bar chart: significant genes / transcripts at FDR < 0.05
summary_bar <- data.table(
  Category = c(
    paste0("DTU genes\n(", opt$treatment, " vs ", opt$control, ")"),
    paste0("DTU transcripts\n(", opt$treatment, " vs ", opt$control, ")")
  ),
  Count = c(n_sig_genes, n_sig_transcripts)
)
p_sum <- ggplot(summary_bar, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 0.5, show.legend = FALSE) +
  scale_fill_manual(values = c("#1f77b4","#ff7f0e")) +
  geom_text(aes(label = Count), vjust = -0.4, size = 4) +
  labs(title = "Significant DTU at FDR < 0.05",
       x = NULL, y = "Count") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(size = 10))
print(p_sum)

dev.off()

message("Done. DTU table: ", opt$out_dtu)
message("      Plots PDF: ", opt$out_plots)
