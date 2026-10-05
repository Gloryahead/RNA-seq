#!/usr/bin/env Rscript
# pathway_analysis.R — GO ORA + KEGG ORA organised into GO-ALL / GO-UP / GO-DOWN
# Output structure per {outdir}:
#   GO-ALL/  — all significant DEGs (up + down)
#   GO-UP/   — upregulated DEGs only
#   GO-DOWN/ — downregulated DEGs only
# Each directory contains:
#   output-GO-Enrichment.{csv,txt}
#   output-KEGG.{csv,txt}
#   output-{GO,BiologicalProcess,CellularComponent,MolecularFunction}-{dotplot,barplot,cnetplot,emapplot}.pdf
#   output-kegg-{dotplot,barplot}.pdf

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(DOSE)
})

opt_list <- list(
  make_option("--degs",       type="character"),
  make_option("--organism",   type="character", default="mouse"),
  make_option("--pval",       type="double",    default=0.05),
  make_option("--qval",       type="double",    default=0.2),
  make_option("--msigdb_gmt", type="character", default=""),
  make_option("--run_go",     type="logical",   default=TRUE),
  make_option("--run_kegg",   type="logical",   default=TRUE),
  make_option("--run_gsea",   type="logical",   default=TRUE),
  make_option("--outdir",     type="character", default="results/pathways"),
  make_option("--label",      type="character", default="comparison")
)
opt <- parse_args(OptionParser(option_list=opt_list))

# ── Organism database ─────────────────────────────────────────────────
org_db <- switch(opt$organism,
  human = { suppressPackageStartupMessages(library(org.Hs.eg.db)); org.Hs.eg.db },
  mouse = { suppressPackageStartupMessages(library(org.Mm.eg.db)); org.Mm.eg.db },
  rat   = { suppressPackageStartupMessages(library(org.Rn.eg.db)); org.Rn.eg.db },
  stop("Unknown organism: ", opt$organism, ". Use human, mouse, or rat.")
)
kegg_organism <- switch(opt$organism, human="hsa", mouse="mmu", rat="rno")

# ── Load DEG table ────────────────────────────────────────────────────
degs <- fread(opt$degs, data.table=FALSE)
stopifnot("gene" %in% colnames(degs), "FDR" %in% colnames(degs))
lfc_col <- intersect(c("log2FC","log2FoldChange","logFC"), colnames(degs))[1]

gene_symbols <- degs$gene
from_type <- if (any(startsWith(na.omit(gene_symbols)[seq_len(min(5, sum(!is.na(gene_symbols))))], "ENS")))
  "ENSEMBL" else "SYMBOL"
message("Detected gene ID type: ", from_type)
entrez_map <- suppressMessages(
  bitr(gene_symbols, fromType=from_type, toType="ENTREZID", OrgDb=org_db)
)
degs <- merge(degs, entrez_map, by.x="gene", by.y=from_type, all.x=FALSE)

sig_all  <- degs$ENTREZID[!is.na(degs$FDR) & degs$FDR < opt$pval]
sig_up   <- degs$ENTREZID[!is.na(degs$FDR) & degs$FDR < opt$pval & degs[[lfc_col]] > 0]
sig_down <- degs$ENTREZID[!is.na(degs$FDR) & degs$FDR < opt$pval & degs[[lfc_col]] < 0]
universe <- degs$ENTREZID
fc_vec   <- setNames(degs[[lfc_col]], degs$ENTREZID)

# Symbol-keyed fc_vec for cnetplot — enrichGO(readable=TRUE) uses gene symbols,
# newer enrichplot versions require foldChange names to match the readable names.
sym_map <- suppressMessages(tryCatch(
  bitr(degs$ENTREZID, fromType="ENTREZID", toType="SYMBOL", OrgDb=org_db),
  error=function(e) NULL
))
fc_vec_sym <- if (!is.null(sym_map)) {
  tmp <- merge(data.frame(ENTREZID=degs$ENTREZID, fc=degs[[lfc_col]]),
               sym_map, by="ENTREZID", all.x=FALSE)
  v <- setNames(tmp$fc, tmp$SYMBOL)
  v[!duplicated(names(v))]
} else fc_vec

message(length(sig_all), " significant genes (",
        length(sig_up), " up, ", length(sig_down), " down)")

# ── Helpers ───────────────────────────────────────────────────────────
save_plot <- function(p, path, w=8, h=6) {
  tryCatch(ggsave(path, p, width=w, height=h),
           error=function(e) message("Plot failed: ", conditionMessage(e)))
}

# Ensure sentinel text files always exist so Snakemake never fails on missing output
ensure_file <- function(path, msg) {
  if (!file.exists(path)) writeLines(msg, path)
}

# ── Core analysis function ────────────────────────────────────────────
run_analysis <- function(gene_set, subdir) {
  dir.create(subdir, showWarnings=FALSE, recursive=TRUE)

  # GO ─────────────────────────────────────────────────────────────────
  if (isTRUE(opt$run_go)) {
    ont_specs <- list(
      list(ont="ALL", prefix="GO",                title="GO All Ontologies"),
      list(ont="BP",  prefix="BiologicalProcess", title="GO Biological Process"),
      list(ont="CC",  prefix="CellularComponent", title="GO Cellular Component"),
      list(ont="MF",  prefix="MolecularFunction", title="GO Molecular Function")
    )

    for (spec in ont_specs) {
      pfx <- spec$prefix

      if (length(gene_set) < 5) {
        message("Too few genes for GO ", spec$ont, " — skipping")
        if (spec$ont == "ALL") {
          writeLines("No significant GO terms (fewer than 5 genes)",
                     file.path(subdir, "output-GO-Enrichment.txt"))
          write.csv(data.frame(), file.path(subdir, "output-GO-Enrichment.csv"), row.names=FALSE)
        }
        next
      }

      ego <- tryCatch(
        suppressMessages(
          enrichGO(gene=gene_set, universe=universe, OrgDb=org_db,
                   ont=spec$ont, pAdjustMethod="BH",
                   pvalueCutoff=opt$pval, qvalueCutoff=opt$qval, readable=TRUE)
        ),
        error=function(e) { message("GO ", spec$ont, " error: ", conditionMessage(e)); NULL }
      )

      if (spec$ont == "ALL") {
        if (!is.null(ego) && nrow(ego) > 0) {
          df <- as.data.frame(ego)
          write.csv(df, file.path(subdir, "output-GO-Enrichment.csv"), row.names=FALSE)
          write.table(df, file.path(subdir, "output-GO-Enrichment.txt"),
                      sep="\t", row.names=FALSE, quote=FALSE)
        } else {
          writeLines("No significant GO terms found",
                     file.path(subdir, "output-GO-Enrichment.txt"))
          write.csv(data.frame(), file.path(subdir, "output-GO-Enrichment.csv"), row.names=FALSE)
        }
      }

      if (!is.null(ego) && nrow(ego) > 0) {
        tryCatch({
          p <- dotplot(ego, showCategory=20, title=spec$title) + theme_classic(base_size=11)
          save_plot(p, file.path(subdir, paste0("output-", pfx, "-dotplot.pdf")))
        }, error=function(e) message("dotplot failed (", spec$ont, "): ", conditionMessage(e)))

        tryCatch({
          p_bar <- barplot(ego, showCategory=20, title=spec$title) + theme_classic(base_size=10)
          save_plot(p_bar, file.path(subdir, paste0("output-", pfx, "-barplot.pdf")), w=10, h=7)
        }, error=function(e) message("barplot failed (", spec$ont, "): ", conditionMessage(e)))

        tryCatch({
          p_cnet <- cnetplot(ego, showCategory=6, foldChange=fc_vec_sym, circular=FALSE)
          save_plot(p_cnet, file.path(subdir, paste0("output-", pfx, "-cnetplot.pdf")), w=12, h=10)
        }, error=function(e) message("cnetplot failed (", spec$ont, "): ", conditionMessage(e)))

        tryCatch({
          ego2   <- pairwise_termsim(ego)
          p_emap <- emapplot(ego2, showCategory=30)
          save_plot(p_emap, file.path(subdir, paste0("output-", pfx, "-emapplot.pdf")), w=12, h=10)
        }, error=function(e) message("emapplot failed (", spec$ont, "): ", conditionMessage(e)))
      }
    }
  }

  # KEGG ───────────────────────────────────────────────────────────────
  if (isTRUE(opt$run_kegg)) {
    if (length(gene_set) < 5) {
      writeLines("No significant KEGG terms (fewer than 5 genes)",
                 file.path(subdir, "output-KEGG.txt"))
      write.csv(data.frame(), file.path(subdir, "output-KEGG.csv"), row.names=FALSE)
    } else {
      ekegg <- tryCatch({
        ek <- enrichKEGG(gene=gene_set, organism=kegg_organism,
                         pvalueCutoff=opt$pval, qvalueCutoff=opt$qval)
        if (!is.null(ek) && nrow(ek) > 0) setReadable(ek, OrgDb=org_db, keyType="ENTREZID") else ek
      }, error=function(e) { message("KEGG error: ", conditionMessage(e)); NULL })

      if (!is.null(ekegg) && nrow(ekegg) > 0) {
        df_kegg <- as.data.frame(ekegg)
        write.csv(df_kegg, file.path(subdir, "output-KEGG.csv"), row.names=FALSE)
        write.table(df_kegg, file.path(subdir, "output-KEGG.txt"),
                    sep="\t", row.names=FALSE, quote=FALSE)

        tryCatch({
          p <- dotplot(ekegg, showCategory=20, title="KEGG") + theme_classic(base_size=11)
          save_plot(p, file.path(subdir, "output-kegg-dotplot.pdf"))
        }, error=function(e) message("KEGG dotplot failed: ", conditionMessage(e)))

        tryCatch({
          p_bar <- barplot(ekegg, showCategory=20, title="KEGG") + theme_classic(base_size=10)
          save_plot(p_bar, file.path(subdir, "output-kegg-barplot.pdf"), w=10, h=7)
        }, error=function(e) message("KEGG barplot failed: ", conditionMessage(e)))
      } else {
        writeLines("No significant KEGG terms found", file.path(subdir, "output-KEGG.txt"))
        write.csv(data.frame(), file.path(subdir, "output-KEGG.csv"), row.names=FALSE)
      }
    }
  }

  ensure_file(file.path(subdir, "output-GO-Enrichment.txt"), "GO analysis not run")
  ensure_file(file.path(subdir, "output-KEGG.txt"),          "KEGG analysis not run")
}

# ── Run for ALL, UP, DOWN ─────────────────────────────────────────────
run_analysis(sig_all,  file.path(opt$outdir, "GO-ALL"))
run_analysis(sig_up,   file.path(opt$outdir, "GO-UP"))
run_analysis(sig_down, file.path(opt$outdir, "GO-DOWN"))

# ── GSEA (ranked, GO BP) — written to root outdir as before ──────────
if (isTRUE(opt$run_gsea)) {
  ranked_list <- setNames(degs[[lfc_col]], degs$ENTREZID)
  ranked_list <- sort(ranked_list[!is.na(ranked_list)], decreasing=TRUE)
  ranked_list <- ranked_list[!duplicated(names(ranked_list))]

  if (length(ranked_list) >= 10) {
    message("Running GSEA (GO BP)...")
    gsea_res <- tryCatch(
      gseGO(geneList=ranked_list, OrgDb=org_db, ont="BP",
            minGSSize=10, maxGSSize=500,
            pvalueCutoff=opt$pval, pAdjustMethod="BH", verbose=FALSE),
      error=function(e) { message("GSEA error: ", conditionMessage(e)); NULL }
    )
    if (!is.null(gsea_res) && nrow(gsea_res) > 0) {
      fwrite(as.data.frame(gsea_res),
             file.path(opt$outdir, paste0(opt$label, "_GSEA_results.tsv")), sep="\t")
      tryCatch({
        p_dot <- dotplot(gsea_res, showCategory=15, split=".sign") +
                 facet_grid(.~.sign) + theme_classic(base_size=10)
        save_plot(p_dot, file.path(opt$outdir, paste0(opt$label, "_GSEA_dotplot.pdf")), w=10, h=7)
      }, error=function(e) message("GSEA dotplot failed: ", conditionMessage(e)))
    }
  }
}

message("Pathway analysis complete. Results in: ", opt$outdir)
