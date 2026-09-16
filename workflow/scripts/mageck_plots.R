#!/usr/bin/env Rscript
# mageck_plots.R — MAGeCK CRISPR screen visualization
# Called by mageck_plot rule (10_crispr.smk)

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(RColorBrewer)
  library(clusterProfiler)
})

option_list <- list(
  make_option("--indir",    type="character", help="Directory with .gene_summary.txt files"),
  make_option("--outdir",   type="character", default="results/crispr/figures"),
  make_option("--organism", type="character", default="human")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

org_db_map <- list(
  human="org.Hs.eg.db", mouse="org.Mm.eg.db", rat="org.Rn.eg.db"
)
org_db <- org_db_map[[opt$organism]]
if (!requireNamespace(org_db, quietly=TRUE)) {
  warning("OrgDb package not available: ", org_db, "; skipping GSEA")
  org_db <- NULL
}

gene_summaries <- list.files(opt$indir, pattern="\\.gene_summary\\.txt$", full.names=TRUE)

for (gs_file in gene_summaries) {
  comp <- sub("\\.gene_summary\\.txt$", "", basename(gs_file))
  gs   <- read.table(gs_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)

  # MAGeCK column names: id, num, neg.score, neg.p.value, neg.fdr, ...
  # pos/neg refer to enrichment/depletion
  if (!"id" %in% colnames(gs)) next
  gs$neg_log10_fdr <- -log10(gs$neg.fdr + 1e-300)
  gs$pos_log10_fdr <- -log10(gs$pos.fdr + 1e-300)
  gs$label         <- ifelse(gs$neg.fdr < 0.1 | gs$pos.fdr < 0.1, gs$id, "")

  # Volcano: pos score (y) vs neg score (x)
  p_vol <- ggplot(gs, aes(neg_log10_fdr, pos_log10_fdr, label=label)) +
    geom_point(
      aes(color=case_when(
        neg.fdr < 0.1 ~ "Depleted",
        pos.fdr < 0.1 ~ "Enriched",
        TRUE ~ "ns"
      )),
      size=1.5, alpha=0.7
    ) +
    geom_text_repel(size=2.5, max.overlaps=20) +
    scale_color_manual(values=c(Depleted="#2166AC", Enriched="#D6604D", ns="grey60")) +
    labs(
      title = paste0("CRISPR Screen: ", comp),
      x = "-log10 FDR (negative selection)",
      y = "-log10 FDR (positive selection)",
      color = NULL
    ) +
    theme_bw(base_size=12)
  ggsave(file.path(opt$outdir, paste0(comp, "_volcano.pdf")), p_vol, width=7, height=6)

  # Rank plot: top enriched + top depleted
  top_dep <- gs %>% arrange(neg.fdr)  %>% slice_head(n=20)
  top_enr <- gs %>% arrange(pos.fdr)  %>% slice_head(n=20)
  plot_df  <- bind_rows(
    mutate(top_dep, direction="Depleted", score=neg.score),
    mutate(top_enr, direction="Enriched", score=pos.score)
  )
  p_rank <- ggplot(plot_df, aes(reorder(id, score), score, fill=direction)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values=c(Depleted="#2166AC", Enriched="#D6604D")) +
    labs(title=paste0(comp, " — Top hits"), x=NULL, y="MAGeCK score") +
    theme_bw(base_size=11) +
    facet_wrap(~direction, scales="free")
  ggsave(file.path(opt$outdir, paste0(comp, "_top_hits.pdf")), p_rank, width=10, height=6)

  # GSEA on depleted genes (using FDR rank)
  if (!is.null(org_db)) {
    tryCatch({
      library(org_db, character.only=TRUE)
      gene_list <- setNames(-log10(gs$neg.fdr + 1e-300), gs$id)
      gene_list <- sort(gene_list, decreasing=TRUE)
      gsea_res  <- gseGO(
        geneList   = gene_list,
        OrgDb      = get(org_db),
        ont        = "BP",
        keyType    = "SYMBOL",
        minGSSize  = 15,
        maxGSSize  = 500,
        pvalueCutoff = 0.1,
        verbose    = FALSE
      )
      if (nrow(gsea_res@result) > 0) {
        p_gsea <- dotplot(gsea_res, showCategory=15,
                          title=paste0(comp, " — GSEA GO BP (depleted)"))
        ggsave(file.path(opt$outdir, paste0(comp, "_gsea_bp.pdf")), p_gsea, width=10, height=8)
      }
    }, error=function(e) message("GSEA skipped for ", comp, ": ", conditionMessage(e)))
  }

  message("Plots written for: ", comp)
}
message("mageck_plots.R complete. Outputs in: ", opt$outdir)
