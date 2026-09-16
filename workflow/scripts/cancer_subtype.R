#!/usr/bin/env Rscript
# cancer_subtype.R — PAM50 subtyping, GSVA, CIBERSORTx prep
# Tutorial ref: Parts 5/6/11 of ngs101.com RNA-seq series
# --mode: pam50 | gsva | cibersortx_prep

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(pheatmap)
  library(DESeq2)
  library(GSVA)
})

option_list <- list(
  make_option("--mode",     type="character", help="pam50 | gsva | cibersortx_prep"),
  make_option("--counts",   type="character", help="Raw counts TSV (featureCounts)"),
  make_option("--samples",  type="character", help="Sample metadata TSV"),
  make_option("--outdir",   type="character", default="results/cancer_subtype"),
  make_option("--organism", type="character", default="human")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

# ── Load data ─────────────────────────────────────────────────────────
meta <- read.table(opt$samples, sep="\t", header=TRUE, stringsAsFactors=FALSE)
rownames(meta) <- meta$sample

raw <- read.table(opt$counts, sep="\t", header=TRUE, skip=1, check.names=FALSE)
count_cols <- grep("\\.bam$", colnames(raw), value=TRUE)
mat <- as.matrix(raw[, count_cols])
colnames(mat) <- sub(".*/(.*)\\.bam$", "\\1", colnames(mat))
rownames(mat) <- raw[[1]]
mat <- mat[, meta$sample]

# log2 TPM (approx from VST)
dds <- DESeqDataSetFromMatrix(mat, colData=meta, design=~group)
dds <- estimateSizeFactors(dds)
expr <- assay(varianceStabilizingTransformation(dds, blind=TRUE))

# ── PAM50 mode ────────────────────────────────────────────────────────
if (opt$mode == "pam50") {
  if (!requireNamespace("PAM50", quietly=TRUE)) {
    stop("PAM50 package not installed. Run: devtools::install_github('ccchang0111/PAM50')")
  }
  library(PAM50)

  # PAM50 expects a matrix of genes x samples with HGNC gene symbols
  pam50_res  <- intrinsic.cluster.predict(
    sbt.model   = pam50.robust,
    data        = expr,
    annot       = data.frame(Gene.Symbol=rownames(expr)),
    do.mapping  = TRUE,
    verbose     = FALSE
  )
  subtypes_df <- data.frame(
    sample  = colnames(expr),
    subtype = pam50_res$subtype
  )
  write.table(subtypes_df, file.path(opt$outdir, "pam50_subtypes.tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)

  # Heatmap of PAM50 gene expression
  pam50_genes <- intersect(pam50.robust$centroids.map$Gene.Symbol, rownames(expr))
  pdf(file.path(opt$outdir, "pam50_heatmap.pdf"), width=10, height=8)
  pheatmap(
    expr[pam50_genes, ],
    annotation_col = data.frame(
      Subtype = subtypes_df$subtype,
      Group   = meta$group,
      row.names = meta$sample
    ),
    scale  = "row",
    show_rownames = FALSE,
    main   = "PAM50 Gene Expression"
  )
  dev.off()
  message("PAM50 subtyping complete → ", file.path(opt$outdir, "pam50_subtypes.tsv"))
}

# ── GSVA mode ─────────────────────────────────────────────────────────
if (opt$mode == "gsva") {
  if (!requireNamespace("msigdbr", quietly=TRUE)) {
    stop("msigdbr not installed. Run: install.packages('msigdbr')")
  }
  library(msigdbr)

  species_map <- list(human="Homo sapiens", mouse="Mus musculus", rat="Rattus norvegicus")
  species_name <- species_map[[opt$organism]]

  # Hallmark + KEGG gene sets
  hallmark <- msigdbr(species=species_name, category="H")
  kegg     <- msigdbr(species=species_name, category="C2", subcategory="CP:KEGG")
  gsets_df <- rbind(hallmark, kegg)
  gene_sets <- split(gsets_df$gene_symbol, gsets_df$gs_name)

  gsva_res <- gsva(
    expr = expr,
    gset.idx.list = gene_sets,
    method   = "gsva",
    kcdf     = "Gaussian",
    verbose  = FALSE
  )
  write.table(gsva_res, file.path(opt$outdir, "gsva_scores.tsv"),
              sep="\t", quote=FALSE)

  # Top variable gene sets heatmap
  gsva_var <- apply(gsva_res, 1, var)
  top_sets <- head(names(sort(gsva_var, decreasing=TRUE)), 50)
  pdf(file.path(opt$outdir, "gsva_heatmap.pdf"), width=14, height=12)
  pheatmap(
    gsva_res[top_sets, ],
    annotation_col = data.frame(Group=meta$group, row.names=meta$sample),
    show_rownames  = TRUE,
    fontsize_row   = 6,
    main           = "GSVA — Top 50 Variable Gene Sets"
  )
  dev.off()
  message("GSVA complete → ", file.path(opt$outdir, "gsva_scores.tsv"))
}

# ── CIBERSORTx prep mode ──────────────────────────────────────────────
if (opt$mode == "cibersortx_prep") {
  # CIBERSORTx requires a matrix: GeneSymbol | Sample1 | Sample2 | ...
  out <- data.frame(GeneSymbol=rownames(expr), expr, check.names=FALSE)
  out_file <- file.path(opt$outdir, "cibersortx_mixture.tsv")
  write.table(out, out_file, sep="\t", row.names=FALSE, quote=FALSE)
  message(
    "CIBERSORTx mixture matrix written to: ", out_file, "\n",
    "Upload at: https://cibersortx.stanford.edu\n",
    "Select 'Impute Cell Fractions' and upload the mixture file."
  )
}
