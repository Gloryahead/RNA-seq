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
#   MAPT_RI_PSI_barplot.pdf     — intron retention PSI per condition (same style as A5SS)
#   MAPT_splicing_efficiency.pdf— intron retention flipped to splicing efficiency (%)
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
file_list <- unlist(lapply(COMPARISONS, function(comp)
  lapply(EVENT_TYPES, function(ev) load_events(comp, ev))),
  recursive = FALSE)
file_list <- Filter(Negate(is.null), file_list)
raw <- as.data.frame(rbindlist(file_list, fill = TRUE, use.names = TRUE))

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

mapt <- raw[!is.na(raw$geneSymbol) & raw$geneSymbol == GENE_OF_INTEREST &
              !is.na(raw$comparison) & valid_psi, mapt_cols]
mapt <- mapt[order(mapt$comparison, mapt$event_type), ]

message("\n── ", GENE_OF_INTEREST, " events (all comparisons, all event types) ──")
print(mapt[, c("comparison", "event_type", "PSI1", "PSI2", "dPSI",
                "IJC_SAMPLE_1", "SJC_SAMPLE_1")], row.names = FALSE)
fwrite(mapt, file.path(OUTDIR, paste0(GENE_OF_INTEREST, "_events.tsv")), sep = "\t")

# ── Plot: MAPT PSI — 3 conditions side by side per event ──────────────────────
# Extract PSI for each condition directly from each comparison:
#   Control = IncLevel2 from Disease_vs_Control
#   Disease = IncLevel1 from Disease_vs_Control
#   Remedy  = IncLevel1 from Remedy_vs_Control
#
# IMPORTANT: merge on event_type+chr+strand alone creates a cartesian product when
# there are multiple MAPT events of the same type on the same strand (e.g. multiple SE
# events on chr17:+). We build a coordinate-based event_key to match events 1-to-1.
if (nrow(mapt) > 0) {
  # Build event_key from the first non-NA exon-start coordinate for each event type:
  #   SE         → exonStart_0base
  #   RI         → riExonStart_0base
  #   A3SS/A5SS  → longExonStart_0base
  #   MXE        → 1stExonStart_0base (R may rename to X1stExonStart_0base)
  coord_val <- rep("", nrow(mapt))
  if ("exonStart_0base"     %in% colnames(mapt))
    coord_val <- ifelse(!is.na(mapt$exonStart_0base),
                        as.character(mapt$exonStart_0base), coord_val)
  if ("riExonStart_0base"   %in% colnames(mapt))
    coord_val <- ifelse(!is.na(mapt$riExonStart_0base),
                        as.character(mapt$riExonStart_0base), coord_val)
  if ("longExonStart_0base" %in% colnames(mapt))
    coord_val <- ifelse(!is.na(mapt$longExonStart_0base),
                        as.character(mapt$longExonStart_0base), coord_val)
  mxe_col <- intersect(c("X1stExonStart_0base", "1stExonStart_0base"), colnames(mapt))
  if (length(mxe_col) > 0)
    coord_val <- ifelse(!is.na(mapt[[mxe_col[1]]]),
                        as.character(mapt[[mxe_col[1]]]), coord_val)

  mapt$event_key <- paste(mapt$event_type, mapt$chr, mapt$strand, coord_val, sep = "_")

  dc_m <- mapt[mapt$comparison == "Disease_vs_Control",
               c("event_key", "event_type", "chr", "strand", "PSI1", "PSI2")]
  rc_m <- mapt[mapt$comparison == "Remedy_vs_Control",
               c("event_key", "PSI1")]

  # 1-to-1 merge; keep all Disease_vs_Control events, NA Remedy if not detected
  mapt3 <- merge(dc_m, rc_m, by = "event_key", all.x = TRUE)
  names(mapt3)[names(mapt3) == "PSI1.x"] <- "PSI_Disease"
  names(mapt3)[names(mapt3) == "PSI2"]   <- "PSI_Control"
  names(mapt3)[names(mapt3) == "PSI1.y"] <- "PSI_Remedy"

  if (nrow(mapt3) > 0) {
    # Build readable per-event label (event type + genomic position of the key exon)
    mapt3$pos_label <- mapply(function(key, etype, chr, strand) {
      prefix <- paste0(etype, "_", chr, "_", strand, "_")
      sub(prefix, "", key, fixed = TRUE)
    }, mapt3$event_key, mapt3$event_type, mapt3$chr, mapt3$strand)
    mapt3$label <- paste0(mapt3$event_type, "\n", mapt3$chr, ":", mapt3$strand,
                          "\n@", mapt3$pos_label)
    mapt3$label <- factor(mapt3$label, levels = unique(mapt3$label))

    # Reshape to long format
    mapt_long <- rbind(
      data.frame(label = mapt3$label, event_type = mapt3$event_type,
                 condition = "Control", PSI_pct = mapt3$PSI_Control * 100),
      data.frame(label = mapt3$label, event_type = mapt3$event_type,
                 condition = "Disease", PSI_pct = mapt3$PSI_Disease * 100),
      data.frame(label = mapt3$label, event_type = mapt3$event_type,
                 condition = "Remedy",  PSI_pct = mapt3$PSI_Remedy  * 100)
    )
    mapt_long$condition <- factor(mapt_long$condition,
                                  levels = c("Control", "Disease", "Remedy"))

    p_mapt <- ggplot(mapt_long, aes(x = condition, y = PSI_pct, fill = condition)) +
      geom_col(width = 0.65, color = "white", na.rm = TRUE) +
      geom_text(aes(label = ifelse(is.na(PSI_pct), "",
                                   paste0(round(PSI_pct, 1), "%"))),
                vjust = -0.4, size = 2.8, na.rm = TRUE) +
      facet_wrap(~ label, scales = "free_y") +
      scale_fill_manual(values = c(Control = "#4E79A7",
                                   Disease = "#E15759",
                                   Remedy  = "#59A14F")) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
      labs(title    = paste(GENE_OF_INTEREST, "- PSI (%) per condition"),
           subtitle = "Each panel = one splicing event (matched by genomic position). n=1, exploratory.",
           y = "PSI (%)", x = NULL) +
      theme_classic(base_size = 11) +
      theme(legend.position  = "bottom",
            strip.background = element_blank(),
            strip.text       = element_text(size = 7, face = "bold"),
            axis.text.x      = element_text(angle = 30, hjust = 1))

    out_mapt <- file.path(OUTDIR, paste0(GENE_OF_INTEREST, "_PSI_barplot.pdf"))
    ggsave(out_mapt, p_mapt, width = 14, height = 8)
    message("Saved: ", out_mapt)
  } else {
    message("No matching Disease_vs_Control + Remedy_vs_Control events for ", GENE_OF_INTEREST)
  }
}

# ── Plot: MAPT splicing efficiency from intron retention (RI) events ──────────
# Splicing efficiency = (1 - RI_PSI) × 100%
# High bar = intron efficiently spliced out; low bar = intron retained / splicing stalled
# Uses ALL MAPT RI events (no dPSI threshold) so every intron is shown.
mapt_ri <- raw[!is.na(raw$geneSymbol) & raw$geneSymbol == GENE_OF_INTEREST &
               !is.na(raw$event_type)  & raw$event_type == "RI" & valid_psi, ]

if (nrow(mapt_ri) > 0) {
  ri_coord <- if ("riExonStart_0base" %in% colnames(mapt_ri))
    ifelse(!is.na(mapt_ri$riExonStart_0base),
           as.character(mapt_ri$riExonStart_0base), "")
  else rep("", nrow(mapt_ri))

  mapt_ri$event_key <- paste("RI", mapt_ri$chr, mapt_ri$strand, ri_coord, sep = "_")

  dc_ri <- mapt_ri[mapt_ri$comparison == "Disease_vs_Control",
                   c("event_key", "chr", "strand", "PSI1", "PSI2")]
  rc_ri <- mapt_ri[mapt_ri$comparison == "Remedy_vs_Control",
                   c("event_key", "PSI1")]

  ri3 <- merge(dc_ri, rc_ri, by = "event_key", all.x = TRUE)
  names(ri3)[names(ri3) == "PSI1.x"] <- "RET_Disease"
  names(ri3)[names(ri3) == "PSI2"]   <- "RET_Control"
  names(ri3)[names(ri3) == "PSI1.y"] <- "RET_Remedy"

  if (nrow(ri3) > 0) {
    ri3$SE_Control <- (1 - ri3$RET_Control) * 100
    ri3$SE_Disease <- (1 - ri3$RET_Disease) * 100
    ri3$SE_Remedy  <- (1 - ri3$RET_Remedy)  * 100

    ri3$pos <- mapply(function(key, chr, strand) {
      sub(paste0("RI_", chr, "_", strand, "_"), "", key, fixed = TRUE)
    }, ri3$event_key, ri3$chr, ri3$strand)
    ri3$label <- paste0("Intron\n", ri3$chr, ":", ri3$strand, "\n@", ri3$pos)
    ri3$label <- factor(ri3$label, levels = unique(ri3$label))

    ri_long <- rbind(
      data.frame(label = ri3$label, condition = "Control", efficiency = ri3$SE_Control),
      data.frame(label = ri3$label, condition = "Disease", efficiency = ri3$SE_Disease),
      data.frame(label = ri3$label, condition = "Remedy",
                 efficiency = ifelse(is.na(ri3$SE_Remedy), NA, ri3$SE_Remedy))
    )
    ri_long$condition <- factor(ri_long$condition, levels = c("Control", "Disease", "Remedy"))

    # ── Plot 1: Intron Retention PSI — same style as A5SS barplot ─────────
    ret_long <- rbind(
      data.frame(label = ri3$label, condition = "Control",
                 PSI_pct = ri3$RET_Control * 100),
      data.frame(label = ri3$label, condition = "Disease",
                 PSI_pct = ri3$RET_Disease * 100),
      data.frame(label = ri3$label, condition = "Remedy",
                 PSI_pct = ifelse(is.na(ri3$RET_Remedy), NA, ri3$RET_Remedy * 100))
    )
    ret_long$condition <- factor(ret_long$condition, levels = c("Control", "Disease", "Remedy"))

    p_ret <- ggplot(ret_long, aes(x = condition, y = PSI_pct, fill = condition)) +
      geom_col(width = 0.65, color = "white", na.rm = TRUE) +
      geom_text(aes(label = ifelse(is.na(PSI_pct), "",
                                   paste0(round(PSI_pct, 1), "%"))),
                vjust = -0.4, size = 2.8, na.rm = TRUE) +
      facet_wrap(~ label, scales = "free_y") +
      scale_fill_manual(values = c(Control = "#4E79A7",
                                   Disease = "#E15759",
                                   Remedy  = "#59A14F")) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
      labs(title    = paste(GENE_OF_INTEREST, "- Intron Retention PSI (%) per condition"),
           subtitle = paste("Higher bar = more intron retained = less efficient splicing.",
                            "All MAPT introns with valid read counts. n=1, exploratory."),
           y = "Intron Retention PSI (%)", x = NULL) +
      theme_classic(base_size = 11) +
      theme(legend.position  = "bottom",
            strip.background = element_blank(),
            strip.text       = element_text(size = 7, face = "bold"),
            axis.text.x      = element_text(angle = 30, hjust = 1))

    out_ret <- file.path(OUTDIR, paste0(GENE_OF_INTEREST, "_RI_PSI_barplot.pdf"))
    ggsave(out_ret, p_ret, width = 14, height = 8)
    message("Saved: ", out_ret)

    # ── Plot 2: Splicing Efficiency — inverted view ────────────────────────
    p_ri <- ggplot(ri_long, aes(x = condition, y = efficiency, fill = condition)) +
      geom_col(width = 0.65, color = "white", na.rm = TRUE) +
      geom_text(aes(label = ifelse(is.na(efficiency), "",
                                   paste0(round(efficiency, 1), "%"))),
                vjust = -0.4, size = 2.8, na.rm = TRUE) +
      facet_wrap(~ label, scales = "free_y") +
      scale_fill_manual(values = c(Control = "#4E79A7",
                                   Disease = "#E15759",
                                   Remedy  = "#59A14F")) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
      labs(title    = paste(GENE_OF_INTEREST, "- Splicing Efficiency per intron (%)"),
           subtitle = paste("Splicing efficiency = (1 - intron retention PSI) x 100%.",
                            "Higher bar = intron efficiently spliced out. n=1, exploratory."),
           y = "Splicing Efficiency (%)", x = NULL) +
      theme_classic(base_size = 11) +
      theme(legend.position  = "bottom",
            strip.background = element_blank(),
            strip.text       = element_text(size = 7, face = "bold"),
            axis.text.x      = element_text(angle = 30, hjust = 1))

    out_ri <- file.path(OUTDIR, paste0(GENE_OF_INTEREST, "_splicing_efficiency.pdf"))
    ggsave(out_ri, p_ri, width = 14, height = 8)
    message("Saved: ", out_ri)
  }
} else {
  message("No MAPT RI events with valid PSI found.")
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
  # One best event per gene (largest |dPSI_Disease_vs_Control|) to avoid duplicate labels
  rescue_dedup <- rescue[order(abs(rescue$dPSI_Disease_vs_Control), decreasing = TRUE), ]
  rescue_dedup <- rescue_dedup[!duplicated(rescue_dedup$geneSymbol), ]
  top10 <- head(rescue_dedup, 10)
  top10$label <- paste0(top10$geneSymbol, " (", top10$event_type, ")")
  top10$label <- factor(top10$label, levels = rev(unique(top10$label)))

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
