#!/usr/bin/env Rscript
# rmats_summary.R — Summarize rMATS alternative splicing results
# Chunlong Ma RNA-seq (9040-0001): 3 comparisons, n=1 per group
#
# b1/b2 assignments in 06_rmats_chunlong.slurm:
#   Disease_vs_Control : b1 = Disease  , b2 = Control  → IncLevel1=Disease,  IncLevel2=Control
#   Remedy_vs_Control  : b1 = Remedy   , b2 = Control  → IncLevel1=Remedy,   IncLevel2=Control
#   Remedy_vs_Disease  : b1 = Remedy   , b2 = Disease  → IncLevel1=Remedy,   IncLevel2=Disease
#   IncLevelDifference = IncLevel1 - IncLevel2
#
# Usage (on HPC, after activating rnaseq_r_env):
#   Rscript workflow/scripts/rmats_summary.R
#
# Outputs (written to OUTDIR):
#   summary_counts.tsv          — event counts passing filters per comparison/type
#   all_filtered_events.tsv     — full filtered event table across all comparisons
#   rescue_events.tsv           — events changed in Disease and reversed by Remedy
#   MAPT_events.tsv             — all MAPT events (unfiltered)
#   MAPT_PSI_barplot.pdf        — PSI per condition per MAPT event type
#   rescue_scatter.pdf          — dPSI scatter: Disease_vs_Control vs Remedy_vs_Disease
#   top_rescue_barplot.pdf      — PSI bar chart for top 10 rescue events

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

# ── Config ─────────────────────────────────────────────────────────────────────
RMATS_DIR <- "/xdisk/haining/maarowosegbe/Chunlong_Ma_RNA_seq/9040-0001-human-analysis/04.rMATS"
OUTDIR    <- file.path(RMATS_DIR, "summary")
GENE_OF_INTEREST <- "MAPT"
MIN_DPSI  <- 0.10   # |IncLevelDifference| threshold
MIN_READS <- 10     # minimum IJC + SJC per sample

COMPARISONS  <- c("Disease_vs_Control", "Remedy_vs_Control", "Remedy_vs_Disease")
EVENT_TYPES  <- c("SE", "RI", "A3SS", "A5SS", "MXE")

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
message("Output dir: ", OUTDIR)

# ── Load all rMATS JC files ────────────────────────────────────────────────────
load_events <- function(comp, event) {
  f <- file.path(RMATS_DIR, comp, paste0(event, ".MATS.JC.txt"))
  if (!file.exists(f)) {
    warning("Missing: ", f); return(NULL)
  }
  dt <- fread(f, data.table = FALSE)
  dt$comparison <- comp
  dt$event_type <- event
  dt
}

message("Loading rMATS files...")
raw <- rbindlist(lapply(COMPARISONS, function(comp)
  lapply(EVENT_TYPES, function(ev) load_events(comp, ev))),
  fill = TRUE, use.names = TRUE)
raw <- as.data.frame(raw)

message("  Loaded ", nrow(raw), " total events across all comparisons and event types")

# ── Parse key columns ──────────────────────────────────────────────────────────
# For n=1, IJC/SJC are single integers; IncLevel1/2 are single floats
raw$IJC1 <- suppressWarnings(as.numeric(raw$IJC_SAMPLE_1))
raw$SJC1 <- suppressWarnings(as.numeric(raw$SJC_SAMPLE_1))
raw$IJC2 <- suppressWarnings(as.numeric(raw$IJC_SAMPLE_2))
raw$SJC2 <- suppressWarnings(as.numeric(raw$SJC_SAMPLE_2))
raw$PSI1 <- suppressWarnings(as.numeric(raw$IncLevel1))
raw$PSI2 <- suppressWarnings(as.numeric(raw$IncLevel2))
raw$dPSI <- raw$IncLevelDifference  # already numeric

# ── Filter ─────────────────────────────────────────────────────────────────────
valid_psi  <- !is.na(raw$PSI1) & !is.na(raw$PSI2) & raw$PSI1 != -1 & raw$PSI2 != -1
covered    <- (raw$IJC1 + raw$SJC1 >= MIN_READS) & (raw$IJC2 + raw$SJC2 >= MIN_READS)
big_change <- abs(raw$dPSI) >= MIN_DPSI

filtered <- raw[valid_psi & covered & big_change, ]
message("  After filtering (|dPSI| >= ", MIN_DPSI, ", reads >= ", MIN_READS, "): ",
        nrow(filtered), " events")

# ── Summary counts ─────────────────────────────────────────────────────────────
counts <- do.call(rbind, lapply(COMPARISONS, function(comp) {
  do.call(rbind, lapply(EVENT_TYPES, function(ev) {
    sub <- filtered[filtered$comparison == comp & filtered$event_type == ev, ]
    data.frame(comparison  = comp,
               event_type  = ev,
               n_changed   = nrow(sub),
               n_increased = sum(sub$dPSI > 0, na.rm = TRUE),
               n_decreased = sum(sub$dPSI < 0, na.rm = TRUE))
  }))
}))

message("\n── Events per comparison / event type ──")
print(counts, row.names = FALSE)
fwrite(counts, file.path(OUTDIR, "summary_counts.tsv"), sep = "\t")

# Save full filtered table
out_cols <- c("comparison", "event_type", "geneSymbol", "GeneID", "chr", "strand",
              "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2",
              "IncLevel1", "IncLevel2", "IncLevelDifference", "PValue", "FDR")
out_cols <- intersect(out_cols, colnames(filtered))
fwrite(filtered[order(abs(filtered$dPSI), decreasing = TRUE), out_cols],
       file.path(OUTDIR, "all_filtered_events.tsv"), sep = "\t")

# ── Rescue events ──────────────────────────────────────────────────────────────
# Rescue: sign(dPSI_Disease_vs_Control) != sign(dPSI_Remedy_vs_Disease)
# i.e., the change in Disease is reversed by Remedy
dc <- raw[raw$comparison == "Disease_vs_Control" & valid_psi & covered &
            abs(raw$dPSI) >= MIN_DPSI,
          c("geneSymbol", "event_type", "chr", "strand", "PSI1", "PSI2", "dPSI")]
colnames(dc)[5:7] <- c("PSI_Disease", "PSI_Control", "dPSI_Disease_vs_Control")

rd <- raw[raw$comparison == "Remedy_vs_Disease" & valid_psi,
          c("geneSymbol", "event_type", "chr", "strand", "PSI1", "PSI2", "dPSI")]
colnames(rd)[5:7] <- c("PSI_Remedy", "PSI_Disease_rd", "dPSI_Remedy_vs_Disease")

rescue <- merge(dc, rd, by = c("geneSymbol", "event_type", "chr", "strand"))

# Direction reversal: Disease went up → Remedy went down (or vice versa)
rescue <- rescue[sign(rescue$dPSI_Disease_vs_Control) != sign(rescue$dPSI_Remedy_vs_Disease), ]
rescue$PSI_Disease_rd <- NULL
rescue <- rescue[order(abs(rescue$dPSI_Disease_vs_Control), decreasing = TRUE), ]

message("\n── Rescue events: ", nrow(rescue),
        " (changed in Disease, reversed direction by Remedy) ──")
print(head(rescue, 15), row.names = FALSE)
fwrite(rescue, file.path(OUTDIR, "rescue_events.tsv"), sep = "\t")

# ── MAPT events ────────────────────────────────────────────────────────────────
mapt_cols <- c("comparison", "event_type", "chr", "strand",
               "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2",
               "IncLevel1", "IncLevel2", "IncLevelDifference", "PValue", "FDR",
               "PSI1", "PSI2", "dPSI")
mapt_cols <- intersect(mapt_cols, colnames(raw))

mapt <- raw[raw$geneSymbol == GENE_OF_INTEREST & valid_psi, mapt_cols]
mapt <- mapt[order(mapt$comparison, mapt$event_type), ]

message("\n── ", GENE_OF_INTEREST, " events (all comparisons, all event types) ──")
print(mapt[, c("comparison", "event_type", "PSI1", "PSI2", "dPSI",
                "IJC_SAMPLE_1", "SJC_SAMPLE_1")], row.names = FALSE)
fwrite(mapt, file.path(OUTDIR, paste0(GENE_OF_INTEREST, "_events.tsv")), sep = "\t")

# ── Plot: MAPT PSI per condition ────────────────────────────────────────────────
if (nrow(mapt) > 0) {
  # Build one row per (event_type × comparison × condition)
  mapt_long <- rbind(
    # b1 = first-named condition
    data.frame(event_type = mapt$event_type,
               comparison = mapt$comparison,
               condition  = sub("(.*)_vs_.*", "\\1", mapt$comparison),
               PSI_pct    = mapt$PSI1 * 100,
               reads      = mapt$IJC_SAMPLE_1 + mapt$SJC_SAMPLE_1),
    # b2 = second-named condition
    data.frame(event_type = mapt$event_type,
               comparison = mapt$comparison,
               condition  = sub(".*_vs_(.*)", "\\1", mapt$comparison),
               PSI_pct    = mapt$PSI2 * 100,
               reads      = mapt$IJC_SAMPLE_2 + mapt$SJC_SAMPLE_2)
  )
  mapt_long$condition <- factor(mapt_long$condition,
                                levels = c("Control", "Disease", "Remedy"))

  p_mapt <- ggplot(mapt_long, aes(x = condition, y = PSI_pct, fill = condition)) +
    geom_col(width = 0.65, color = "white") +
    geom_text(aes(label = paste0(round(PSI_pct, 1), "%")),
              vjust = -0.4, size = 3) +
    facet_grid(event_type ~ comparison, scales = "free_y") +
    scale_fill_manual(values = c(Control = "#4E79A7",
                                 Disease = "#E15759",
                                 Remedy  = "#59A14F")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(title  = paste(GENE_OF_INTEREST, "— PSI (%) per condition and event type"),
         subtitle = "n=1 per group; values are exploratory",
         y = "PSI (%)", x = NULL) +
    theme_classic(base_size = 11) +
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text       = element_text(size = 9, face = "bold"),
          axis.text.x      = element_text(angle = 30, hjust = 1))

  out_mapt <- file.path(OUTDIR, paste0(GENE_OF_INTEREST, "_PSI_barplot.pdf"))
  ggsave(out_mapt, p_mapt,
         width  = 3 * length(COMPARISONS) + 1,
         height = 2.5 * length(unique(mapt$event_type)) + 1.5)
  message("Saved: ", out_mapt)
}

# ── Plot: rescue scatter (dPSI Disease_vs_Control vs Remedy_vs_Disease) ─────────
dc_all <- raw[raw$comparison == "Disease_vs_Control" & valid_psi & covered, ]
rd_all <- raw[raw$comparison == "Remedy_vs_Disease"  & valid_psi, ]
scatter_df <- merge(
  dc_all[, c("geneSymbol", "event_type", "chr", "strand", "dPSI")],
  rd_all[, c("geneSymbol", "event_type", "chr", "strand", "dPSI")],
  by = c("geneSymbol", "event_type", "chr", "strand"),
  suffixes = c("_DvsC", "_RvsD")
)

# Label top rescue events and MAPT
scatter_df$is_rescue <- abs(scatter_df$dPSI_DvsC) >= MIN_DPSI &
                        sign(scatter_df$dPSI_DvsC) != sign(scatter_df$dPSI_RvsD)
scatter_df$label <- ifelse(scatter_df$geneSymbol == GENE_OF_INTEREST |
                             (scatter_df$is_rescue &
                              rank(-abs(scatter_df$dPSI_DvsC)) <= 10),
                           paste(scatter_df$geneSymbol, scatter_df$event_type), "")

p_scatter <- ggplot(scatter_df,
                    aes(x = dPSI_DvsC, y = dPSI_RvsD,
                        color = is_rescue, alpha = is_rescue)) +
  geom_point(size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept =  c(-MIN_DPSI, MIN_DPSI), linetype = "dotted", color = "grey70") +
  geom_hline(yintercept = c(-MIN_DPSI, MIN_DPSI),  linetype = "dotted", color = "grey70") +
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20, color = "black") +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#E15759"),
                     labels = c("Other", "Rescue candidate")) +
  scale_alpha_manual(values = c("FALSE" = 0.4, "TRUE" = 1)) +
  facet_wrap(~ event_type, ncol = 3) +
  labs(title    = "Rescue scatter: Disease vs Control  ×  Remedy vs Disease",
       subtitle = paste("Q2 (top-left) and Q4 (bottom-right) = Remedy reverses Disease change.",
                        " n=1; exploratory only."),
       x     = "dPSI  Disease vs Control",
       y     = "dPSI  Remedy vs Disease",
       color = NULL, alpha = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position  = "bottom",
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold"))

out_scatter <- file.path(OUTDIR, "rescue_scatter.pdf")
ggsave(out_scatter, p_scatter, width = 12, height = 8)
message("Saved: ", out_scatter)

# ── Plot: top rescue events PSI bar chart ──────────────────────────────────────
if (nrow(rescue) > 0) {
  top10 <- head(rescue, 10)
  top10$label <- paste0(top10$geneSymbol, " (", top10$event_type, ")")
  top10$label <- factor(top10$label, levels = rev(top10$label))

  bar_df <- rbind(
    data.frame(label = top10$label, condition = "Control",
               PSI   = top10$PSI_Control * 100),
    data.frame(label = top10$label, condition = "Disease",
               PSI   = top10$PSI_Disease * 100),
    data.frame(label = top10$label, condition = "Remedy",
               PSI   = top10$PSI_Remedy  * 100)
  )
  bar_df$condition <- factor(bar_df$condition,
                             levels = c("Control", "Disease", "Remedy"))

  p_top <- ggplot(bar_df, aes(x = condition, y = PSI, fill = condition)) +
    geom_col(width = 0.65, color = "white") +
    geom_text(aes(label = paste0(round(PSI, 1), "%")),
              vjust = -0.4, size = 2.8) +
    facet_wrap(~ label, ncol = 5) +
    scale_fill_manual(values = c(Control = "#4E79A7",
                                 Disease = "#E15759",
                                 Remedy  = "#59A14F")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(title    = "Top rescue events — PSI per condition",
         subtitle = "Events with largest dPSI in Disease that are reversed by Remedy (n=1, exploratory)",
         y = "PSI (%)", x = NULL) +
    theme_classic(base_size = 10) +
    theme(legend.position  = "bottom",
          strip.background = element_blank(),
          strip.text       = element_text(size = 8, face = "bold"),
          axis.text.x      = element_text(angle = 30, hjust = 1))

  out_top <- file.path(OUTDIR, "top_rescue_barplot.pdf")
  ggsave(out_top, p_top, width = 14, height = 6)
  message("Saved: ", out_top)
}

message("\n── Done. All outputs in: ", OUTDIR, " ──")
