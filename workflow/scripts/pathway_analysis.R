#!/usr/bin/env Rscript
# pathway_analysis.R — GO ORA, KEGG ORA, and GSEA (Part 5)
# Input:  DEG results TSV (from deg_analysis.R)
# Output: GO/KEGG/GSEA result tables + dot plots
# Usage:  Called by Snakemake rule 03_pathways; standalone:
#   Rscript workflow/scripts/pathway_analysis.R \
#     --degs results/deg/KO_vs_WT_DEG_results.tsv \
#     --organism mouse --outdir results/pathways --label KO_vs_WT

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
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

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
# Standardise column names (output from deg_analysis.R)
stopifnot("gene" %in% colnames(degs), "FDR" %in% colnames(degs))
lfc_col <- intersect(c("log2FC","log2FoldChange","logFC"), colnames(degs))[1]

# Map gene symbols → Entrez IDs
gene_symbols <- degs$gene
entrez_map <- suppressMessages(
  bitr(gene_symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb=org_db)
)
degs <- merge(degs, entrez_map, by.x="gene", by.y="SYMBOL", all.x=FALSE)

# Significant genes for ORA
sig_genes   <- degs$ENTREZID[!is.na(degs$FDR) & degs$FDR < opt$pval]
universe    <- degs$ENTREZID

# Ranked gene list for GSEA (by signed -log10 FDR or log2FC)
ranked_list <- setNames(degs[[lfc_col]], degs$ENTREZID)
ranked_list <- sort(ranked_list[!is.na(ranked_list)], decreasing=TRUE)

message(length(sig_genes), " significant genes for ORA; ",
        length(ranked_list), " genes for GSEA")

save_table <- function(obj, path) {
  if (!is.null(obj) && nrow(as.data.frame(obj)) > 0)
    fwrite(as.data.frame(obj), path, sep="\t")
}

save_plot <- function(p, path, w=8, h=6) {
  tryCatch(ggsave(path, p, width=w, height=h), error=function(e) message("Plot failed: ", e))
}

# ── GO ORA ───────────────────────────────────────────────────────────
if (isTRUE(opt$run_go) && length(sig_genes) >= 5) {
  message("Running GO ORA...")
  ego <- enrichGO(gene         = sig_genes,
                  universe     = universe,
                  OrgDb        = org_db,
                  ont          = "BP",
                  pAdjustMethod= "BH",
                  pvalueCutoff = opt$pval,
                  qvalueCutoff = opt$qval,
                  readable     = TRUE)
  save_table(ego, file.path(opt$outdir, paste0(opt$label, "_GO_results.tsv")))
  if (!is.null(ego) && nrow(ego) > 0) {
    p <- dotplot(ego, showCategory=20, title=paste("GO BP:", opt$label)) +
         theme_classic(base_size=11)
    save_plot(p, file.path(opt$outdir, paste0(opt$label, "_GO_dotplot.pdf")))
  }
}

# ── KEGG ORA ──────────────────────────────────────────────────────────
if (isTRUE(opt$run_kegg) && length(sig_genes) >= 5) {
  message("Running KEGG ORA...")
  ekegg <- enrichKEGG(gene         = sig_genes,
                      organism     = kegg_organism,
                      pvalueCutoff = opt$pval,
                      qvalueCutoff = opt$qval)
  ekegg <- setReadable(ekegg, OrgDb=org_db, keyType="ENTREZID")
  save_table(ekegg, file.path(opt$outdir, paste0(opt$label, "_KEGG_results.tsv")))
  if (!is.null(ekegg) && nrow(ekegg) > 0) {
    p <- dotplot(ekegg, showCategory=20, title=paste("KEGG:", opt$label)) +
         theme_classic(base_size=11)
    save_plot(p, file.path(opt$outdir, paste0(opt$label, "_KEGG_dotplot.pdf")))
  }
}

# ── GSEA (GO Biological Process) ──────────────────────────────────────
if (isTRUE(opt$run_gsea) && length(ranked_list) >= 10) {
  message("Running GSEA (GO BP)...")
  gsea_res <- gseGO(geneList      = ranked_list,
                    OrgDb         = org_db,
                    ont           = "BP",
                    minGSSize     = 10,
                    maxGSSize     = 500,
                    pvalueCutoff  = opt$pval,
                    pAdjustMethod = "BH",
                    verbose       = FALSE)
  save_table(gsea_res, file.path(opt$outdir, paste0(opt$label, "_GSEA_results.tsv")))
  if (!is.null(gsea_res) && nrow(gsea_res) > 0) {
    p_dot  <- dotplot(gsea_res, showCategory=15, split=".sign") +
              facet_grid(.~.sign) + theme_classic(base_size=10)
    save_plot(p_dot, file.path(opt$outdir, paste0(opt$label, "_GSEA_dotplot.pdf")), w=10, h=7)
    p_ridge <- ridgeplot(gsea_res, showCategory=20) + theme_classic(base_size=10)
    save_plot(p_ridge, file.path(opt$outdir, paste0(opt$label, "_GSEA_ridge.pdf")), w=9, h=8)
  }
}

# ── MSigDB GSEA (optional, Part 5) ────────────────────────────────────
if (nchar(opt$msigdb_gmt) > 0 && file.exists(opt$msigdb_gmt)) {
  message("Running GSEA with MSigDB gene sets: ", opt$msigdb_gmt)
  msigdb_sets <- read.gmt(opt$msigdb_gmt)
  ranked_sym   <- setNames(degs[[lfc_col]], degs$gene)
  ranked_sym   <- sort(ranked_sym[!is.na(ranked_sym)], decreasing=TRUE)
  msigdb_gsea  <- GSEA(ranked_sym, TERM2GENE=msigdb_sets,
                        pvalueCutoff=opt$pval, verbose=FALSE)
  save_table(msigdb_gsea, file.path(opt$outdir, paste0(opt$label, "_MSigDB_GSEA_results.tsv")))
}

message("Pathway analysis complete. Results in: ", opt$outdir)
