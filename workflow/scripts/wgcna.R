#!/usr/bin/env Rscript
# wgcna.R — WGCNA co-expression network analysis
# Tutorial ref: WGCNA tutorial (ngs101.com)
# Outputs: module assignments, hub genes, module-trait correlations, STRINGdb plot

suppressPackageStartupMessages({
  library(optparse)
  library(WGCNA)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(STRINGdb)
})

option_list <- list(
  make_option("--counts",   type="character", help="Raw counts matrix (featureCounts TSV)"),
  make_option("--samples",  type="character", help="Sample metadata TSV"),
  make_option("--outdir",   type="character", default="results/networks/wgcna"),
  make_option("--threads",  type="integer",   default=8),
  make_option("--organism", type="character", default="human",
              help="human | mouse | rat"),
  make_option("--min_module_size", type="integer", default=30),
  make_option("--soft_power",      type="integer", default=0,
              help="0 = auto-detect via scale-free fit")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

enableWGCNAThreads(nThreads=opt$threads)

# ── Load and normalize ────────────────────────────────────────────────
meta <- read.table(opt$samples, sep="\t", header=TRUE, stringsAsFactors=FALSE)
rownames(meta) <- meta$sample

raw <- read.table(opt$counts, sep="\t", header=TRUE, skip=1, check.names=FALSE)
count_cols <- grep("\\.bam$", colnames(raw), value=TRUE)
mat <- as.matrix(raw[, count_cols])
colnames(mat) <- sub(".*/(.*)\\.bam$", "\\1", colnames(mat))
rownames(mat) <- raw[[1]]
mat <- mat[, meta$sample]

# VST normalization
dds <- DESeqDataSetFromMatrix(mat, colData=meta, design=~group)
dds <- estimateSizeFactors(dds)
datExpr0 <- t(assay(varianceStabilizingTransformation(dds, blind=TRUE)))

# Filter genes: remove low-variance genes (keep top 75%)
vars <- apply(datExpr0, 2, var)
datExpr <- datExpr0[, vars >= quantile(vars, 0.25)]
message(ncol(datExpr), " genes after variance filtering")

# Check for good genes/samples
gsg <- goodSamplesGenes(datExpr, verbose=0)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  message("Removed ", sum(!gsg$goodSamples), " bad samples, ",
          sum(!gsg$goodGenes), " bad genes")
}

# ── Soft-threshold selection ──────────────────────────────────────────
if (opt$soft_power == 0) {
  powers <- c(seq(1,10), seq(12,20,2))
  sft    <- pickSoftThreshold(datExpr, powerVector=powers, verbose=0)
  # Choose smallest power where scale-free R^2 >= 0.85
  idx    <- which(sft$fitIndices$SFT.R.sq >= 0.85)
  soft_power <- if (length(idx)) powers[idx[1]] else 6
  message("Auto soft power = ", soft_power)
} else {
  soft_power <- opt$soft_power
}

# ── Network construction and module detection ─────────────────────────
net <- blockwiseModules(
  datExpr,
  power            = soft_power,
  TOMType          = "unsigned",
  minModuleSize    = opt$min_module_size,
  reassignThreshold= 0,
  mergeCutHeight   = 0.25,
  numericLabels    = FALSE,
  verbose          = 0
)
moduleColors <- net$colors
MEs          <- net$MEs
message(length(unique(moduleColors)) - 1, " modules detected")

# ── Module-trait correlation ──────────────────────────────────────────
# Encode group as numeric (one-hot for each level)
traits <- model.matrix(~0 + group, data=meta)
rownames(traits) <- meta$sample
nSamples <- nrow(datExpr)
moduleTraitCor <- cor(MEs, traits, use="p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

pdf(file.path(opt$outdir, "module_trait_correlations.pdf"), width=10, height=8)
labeledHeatmap(
  Matrix    = moduleTraitCor,
  xLabels   = colnames(traits),
  yLabels   = rownames(moduleTraitCor),
  ySymbols  = rownames(moduleTraitCor),
  colorLabels = FALSE,
  colors    = blueWhiteRed(50),
  textMatrix= paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep=""),
  setStdMargins = FALSE,
  cex.text  = 0.7,
  zlim      = c(-1,1),
  main      = "Module-Trait Correlations"
)
dev.off()

# ── Hub gene identification ────────────────────────────────────────────
geneModuleMembership <- cor(datExpr, MEs, use="p")
colnames(geneModuleMembership) <- paste0("MM.", colnames(MEs))

# Top 10 hub genes per module (highest module membership)
hub_list <- lapply(unique(moduleColors), function(m) {
  genes <- names(moduleColors)[moduleColors == m]
  mm_col <- paste0("MM.ME", m)
  if (!mm_col %in% colnames(geneModuleMembership)) return(NULL)
  mm_vals <- geneModuleMembership[genes, mm_col]
  top10 <- head(sort(mm_vals, decreasing=TRUE), 10)
  data.frame(module=m, gene=names(top10), module_membership=top10)
})
hub_df <- do.call(rbind, hub_list)
write.table(hub_df, file.path(opt$outdir, "hub_genes.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

# ── Module assignment table ───────────────────────────────────────────
assign_df <- data.frame(
  gene   = names(moduleColors),
  module = moduleColors
)
write.table(assign_df, file.path(opt$outdir, "module_assignments.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

# ── STRINGdb network for the most correlated module ────────────────────
top_module <- names(which.max(abs(moduleTraitCor[,1])))
top_genes  <- names(moduleColors)[moduleColors == sub("ME","", top_module)]

string_version <- "11.5"
species_map <- list(human="9606", mouse="10090", rat="10116")
species_id  <- species_map[[opt$organism]]

tryCatch({
  string_db <- STRINGdb$new(version=string_version, species=as.integer(species_id),
                            score_threshold=400)
  mapped    <- string_db$map(data.frame(gene=top_genes), "gene", removeUnmappedRows=TRUE)
  if (nrow(mapped) > 2) {
    pdf(file.path(opt$outdir, paste0("stringdb_", top_module, ".pdf")), width=10, height=8)
    string_db$plot_network(mapped$STRING_id)
    dev.off()
  }
}, error=function(e) message("STRINGdb skipped: ", conditionMessage(e)))

message("wgcna.R complete. Outputs in: ", opt$outdir)
