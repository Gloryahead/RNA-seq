#!/usr/bin/env Rscript
# rna_editing.R
# Downstream analysis of JACUSA2 RNA editing output.
# Filters A-to-I (A>G / T>C) editing sites from JACUSA2 call-2 differential
# output, produces diagnostic plots, and writes a per-comparison summary.
#
# Usage (called from Snakemake rule rna_edit_analysis):
#   Rscript rna_editing.R \
#       --call2 /path/to/comp_differential.tsv \
#       --call1_dir /path/to/call1_dir \
#       --samples config/samples.tsv \
#       --treatment KRAS_SPIB \
#       --control KRAS \
#       --min_coverage 10 \
#       --min_frequency 0.05 \
#       --min_edited_reads 3 \
#       --known_snps "" \
#       --out_table results/rna_editing/comp_editing_filtered.tsv \
#       --out_plots results/rna_editing/comp_editing_plots.pdf \
#       --out_summary results/rna_editing/comp_editing_summary.tsv

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
})

# ── CLI arguments ─────────────────────────────────────────────────────
opt_list <- list(
  make_option("--call2",            type = "character", help = "JACUSA2 call-2 TSV"),
  make_option("--call1_dir",        type = "character", help = "Directory with per-sample call-1 TSVs"),
  make_option("--samples",          type = "character", help = "Sample sheet TSV"),
  make_option("--treatment",        type = "character", help = "Treatment group name"),
  make_option("--control",          type = "character", help = "Control group name"),
  make_option("--min_coverage",     type = "integer",   default = 10L,
              help = "Min read depth at site in both conditions [default %default]"),
  make_option("--min_frequency",    type = "double",    default = 0.05,
              help = "Min editing frequency [default %default]"),
  make_option("--min_edited_reads", type = "integer",   default = 3L,
              help = "Min edited-read count [default %default]"),
  make_option("--known_snps",       type = "character", default = "",
              help = "Optional VCF of known SNPs to exclude (path or empty)"),
  make_option("--out_table",        type = "character", help = "Output filtered TSV"),
  make_option("--out_plots",        type = "character", help = "Output PDF of plots"),
  make_option("--out_summary",      type = "character", help = "Output summary TSV")
)
opt <- parse_args(OptionParser(option_list = opt_list))

# ── Ensure output directories exist ──────────────────────────────────
dir.create(dirname(opt$out_table),   recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$out_plots),   recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$out_summary), recursive = TRUE, showWarnings = FALSE)

# ── Load sample sheet ─────────────────────────────────────────────────
samples_df <- fread(opt$samples, sep = "\t", header = TRUE, colClasses = "character")
trt_samples <- samples_df[group == opt$treatment, sample]
ctl_samples <- samples_df[group == opt$control,   sample]
all_samples  <- c(trt_samples, ctl_samples)

# ── Parse JACUSA2 call-2 TSV ──────────────────────────────────────────
# Skip comment lines starting with '#'; first non-comment line is the header.
raw_lines <- readLines(opt$call2)
comment_idx <- grep("^#", raw_lines)
if (length(comment_idx) > 0) {
  header_line <- raw_lines[max(comment_idx)]
  # Strip leading '#' from header
  header_line <- sub("^#+\\s*", "", header_line)
  data_lines  <- raw_lines[-seq_len(max(comment_idx))]
} else {
  header_line <- raw_lines[1]
  data_lines  <- raw_lines[-1]
}

if (length(data_lines) == 0) {
  message("No data rows in call-2 output. Writing empty outputs.")
  empty_cols <- c("CHROM","POS","NAME","SCORE","STRAND","INFO","FILTER",
                  "trt_total","trt_A","trt_C","trt_G","trt_T",
                  "ctl_total","ctl_A","ctl_C","ctl_G","ctl_T",
                  "editing_freq_trt","editing_freq_ctl","substitution",
                  "stat","pValue","fdr_adj")
  empty_dt <- data.table(matrix(nrow = 0, ncol = length(empty_cols),
                                dimnames = list(NULL, empty_cols)))
  fwrite(empty_dt, opt$out_table, sep = "\t")
  fwrite(data.table(metric = character(), value = character()),
         opt$out_summary, sep = "\t")
  pdf(opt$out_plots, width = 8, height = 6)
  plot.new(); text(0.5, 0.5, "No editing sites detected.", cex = 1.2)
  dev.off()
  quit(save = "no", status = 0)
}

col_names <- strsplit(header_line, "\t")[[1]]
con <- textConnection(paste(data_lines, collapse = "\n"))
dt  <- fread(con, sep = "\t", header = FALSE, col.names = col_names,
             fill = TRUE, na.strings = c("NA","N/A","."))
close(con)

n_raw <- nrow(dt)
message("Loaded ", n_raw, " sites from call-2 output.")

# ── Identify BASES columns (after FILTER) ────────────────────────────
# JACUSA2 call-2 format:
#   CHROM POS NAME SCORE STRAND INFO FILTER BASES[trt_s1] ... BASES[trt_sN]
#   BASES[ctl_s1] ... BASES[ctl_sN]  [stat  pValue  fdr_adj]
fixed_cols <- c("CHROM","POS","NAME","SCORE","STRAND","INFO","FILTER")
stat_cols  <- intersect(col_names, c("stat","pValue","fdr_adj","p.value","P.value","FDR"))
bases_cols <- setdiff(col_names, c(fixed_cols, stat_cols))

n_trt <- length(trt_samples)
n_ctl <- length(ctl_samples)

# Assign treatment / control BASES columns by position
trt_bases_cols <- bases_cols[seq_len(n_trt)]
ctl_bases_cols <- bases_cols[seq(n_trt + 1, n_trt + n_ctl)]

# ── Parse BASES format: total:A:C:G:T ────────────────────────────────
parse_bases_cols <- function(dt, cols, prefix) {
  # Average across replicate columns → integer counts
  counts <- lapply(cols, function(col) {
    mat <- do.call(rbind, strsplit(as.character(dt[[col]]), ":"))
    mode(mat) <- "numeric"
    mat
  })
  # Element-wise mean across replicates
  avg_mat <- Reduce("+", counts) / length(counts)
  colnames(avg_mat) <- c("total","A","C","G","T")
  for (cn in colnames(avg_mat)) {
    dt[, paste0(prefix, "_", cn) := avg_mat[, cn]]
  }
  dt
}

dt <- parse_bases_cols(dt, trt_bases_cols, "trt")
dt <- parse_bases_cols(dt, ctl_bases_cols, "ctl")

# ── Determine substitution type from dominant edited base ─────────────
# For each site: ref base inferred from strand
# + strand A>G editing: G_count elevated, A is ref
# - strand T>C editing  (same biological event): C_count elevated on - strand
dt[, substitution := {
  # Dominant non-reference base on + strand
  trt_counts <- cbind(A = trt_A, C = trt_C, G = trt_G, T = trt_T)
  ctl_counts <- cbind(A = ctl_A, C = ctl_C, G = ctl_G, T = ctl_T)
  combined   <- trt_counts + ctl_counts
  # Determine ref base: highest-count base overall
  bases       <- c("A","C","G","T")
  ref_base    <- bases[apply(combined, 1, which.max)]
  # Dominant alt base (exclude ref)
  alt_base    <- vapply(seq_len(nrow(combined)), function(i) {
    vec <- combined[i, ]
    vec[bases == ref_base[i]] <- -Inf
    bases[which.max(vec)]
  }, character(1))
  paste0(ref_base, ">", alt_base)
}]

# ── Compute editing frequencies ───────────────────────────────────────
dt[, editing_freq_trt := ifelse(
  STRAND == "+",
  trt_G / pmax(trt_total, 1),
  trt_C / pmax(trt_total, 1)
)]
dt[, editing_freq_ctl := ifelse(
  STRAND == "+",
  ctl_G / pmax(ctl_total, 1),
  ctl_C / pmax(ctl_total, 1)
)]
dt[, edited_reads_trt := ifelse(STRAND == "+", trt_G, trt_C)]
dt[, edited_reads_ctl := ifelse(STRAND == "+", ctl_G, ctl_C)]

# ── Filter sites ──────────────────────────────────────────────────────
dt_filtered <- dt[
  trt_total >= opt$min_coverage &
  ctl_total >= opt$min_coverage &
  (editing_freq_trt >= opt$min_frequency | editing_freq_ctl >= opt$min_frequency) &
  (edited_reads_trt >= opt$min_edited_reads | edited_reads_ctl >= opt$min_edited_reads)
]

n_filtered <- nrow(dt_filtered)
message("Sites after coverage/frequency filter: ", n_filtered)

# Keep only A>G (+ strand) and T>C (- strand) for canonical A-to-I editing
dt_ai <- dt_filtered[substitution %in% c("A>G","T>C")]
n_ai  <- nrow(dt_ai)
message("A-to-I (A>G / T>C) sites: ", n_ai)

# ── Select output columns ─────────────────────────────────────────────
keep_cols <- c(fixed_cols,
               "trt_total","trt_A","trt_C","trt_G","trt_T",
               "ctl_total","ctl_A","ctl_C","ctl_G","ctl_T",
               "editing_freq_trt","editing_freq_ctl",
               "edited_reads_trt","edited_reads_ctl",
               "substitution")
keep_cols <- c(keep_cols, intersect(stat_cols, col_names))
fwrite(dt_filtered[, ..keep_cols], opt$out_table, sep = "\t")

# ── Summary table ─────────────────────────────────────────────────────
ai_frac <- if (n_filtered > 0) round(n_ai / n_filtered, 4) else NA_real_
summary_dt <- data.table(
  metric = c("total_sites_detected",
             "sites_after_filter",
             "A_to_I_sites",
             "A_to_I_fraction"),
  value  = as.character(c(n_raw, n_filtered, n_ai, ai_frac))
)
fwrite(summary_dt, opt$out_summary, sep = "\t")

# ── Plots ─────────────────────────────────────────────────────────────
pdf(opt$out_plots, width = 9, height = 6)

# 1. Bar chart: editing event counts by substitution type
if (nrow(dt_filtered) > 0) {
  sub_counts <- dt_filtered[, .N, by = substitution][order(-N)]
  sub_counts[, is_AI := substitution %in% c("A>G","T>C")]
  p1 <- ggplot(sub_counts, aes(x = reorder(substitution, -N), y = N,
                               fill = is_AI)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = c("TRUE" = "#d62728", "FALSE" = "#aec7e8"),
                      labels = c("TRUE" = "A>I (A>G/T>C)", "FALSE" = "Other"),
                      name   = "Type") +
    labs(title = "Editing site counts by substitution type",
         x = "Substitution", y = "Number of sites") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p1)
} else {
  plot.new(); text(0.5, 0.5, "No sites passed filter.", cex = 1.2)
}

# 2. Scatter plot: editing frequency treatment vs control
if (nrow(dt_filtered) > 0) {
  plot_dt <- copy(dt_filtered)
  has_pval <- "pValue" %in% names(plot_dt)
  if (has_pval) {
    plot_dt[, sig := !is.na(pValue) & pValue < 0.05]
  } else {
    plot_dt[, sig := FALSE]
  }
  p2 <- ggplot(plot_dt, aes(x = editing_freq_ctl, y = editing_freq_trt,
                             colour = sig)) +
    geom_point(alpha = 0.5, size = 1.5) +
    scale_colour_manual(values = c("FALSE" = "grey60", "TRUE" = "#d62728"),
                        labels = c("FALSE" = "Not significant",
                                   "TRUE"  = "p < 0.05"),
                        name   = "") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "black") +
    labs(title = paste("Editing frequency:", opt$treatment, "vs", opt$control),
         x = paste("Editing frequency —", opt$control),
         y = paste("Editing frequency —", opt$treatment)) +
    theme_bw(base_size = 12)
  print(p2)
} else {
  plot.new(); text(0.5, 0.5, "No sites to plot.", cex = 1.2)
}

# 3. Histograms of editing frequency per condition
if (nrow(dt_filtered) > 0) {
  freq_long <- rbind(
    data.table(group = opt$treatment, freq = dt_filtered$editing_freq_trt),
    data.table(group = opt$control,   freq = dt_filtered$editing_freq_ctl)
  )
  p3 <- ggplot(freq_long[is.finite(freq)], aes(x = freq, fill = group)) +
    geom_histogram(bins = 40, alpha = 0.65, position = "identity") +
    scale_fill_manual(values = c("#1f77b4","#ff7f0e"), name = "Group") +
    labs(title = "Distribution of editing frequencies",
         x = "Editing frequency", y = "Number of sites") +
    theme_bw(base_size = 12)
  print(p3)
} else {
  plot.new(); text(0.5, 0.5, "No sites to plot.", cex = 1.2)
}

dev.off()

message("Done. Filtered table:  ", opt$out_table)
message("       Plots PDF:       ", opt$out_plots)
message("       Summary table:   ", opt$out_summary)
