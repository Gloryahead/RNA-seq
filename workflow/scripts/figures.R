#!/usr/bin/env Rscript
# figures.R — Publication-ready figures (Part 4)
# Volcano (EnhancedVolcano), DEG heatmaps (all, top100, top50, up, down)
# Input:  DEG results TSV, DGE/DESeq2 RDS object
# Usage:  Called by Snakemake rule publication_figures

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(EnhancedVolcano)
  library(SummarizedExperiment)
})

opt_list <- list(
  make_option("--degs",             type="character"),
  make_option("--rds",              type="character"),
  make_option("--treatment",        type="character"),
  make_option("--control",          type="character"),
  make_option("--organism",         type="character", default="mouse"),
  make_option("--alpha",            type="double",    default=0.05),
  make_option("--lfc",              type="double",    default=1.0),
  make_option("--out_volcano",      type="character"),
  make_option("--out_heatmap",      type="character"),
  make_option("--out_heatmap_top100",  type="character", default=NULL),
  make_option("--out_heatmap_top50",   type="character", default=NULL),
  make_option("--out_heatmap_up",      type="character", default=NULL),
  make_option("--out_heatmap_down",    type="character", default=NULL),
  make_option("--out_volcano_html",     type="character", default=NULL),
  make_option("--out_volcano_labeled", type="character", default=NULL),
  make_option("--n_label",             type="integer",   default=20,
              help="Top N up + N down DEGs to label [default %default]")
)
opt <- parse_args(OptionParser(option_list=opt_list))

degs    <- fread(opt$degs, data.table=FALSE)
lfc_col <- intersect(c("log2FC","log2FoldChange","logFC"), colnames(degs))[1]
p_col   <- intersect(c("FDR","adj.P.Val","padj"), colnames(degs))[1]
# fread reads all-NA columns as character; coerce back to numeric
degs[[p_col]]   <- suppressWarnings(as.numeric(degs[[p_col]]))
degs[[lfc_col]] <- suppressWarnings(as.numeric(degs[[lfc_col]]))

# ── Add gene symbols if missing (Ensembl IDs detected) ────────────────
if (!"gene_symbol" %in% colnames(degs) &&
    any(startsWith(na.omit(degs$gene)[seq_len(min(5, sum(!is.na(degs$gene))))], "ENS"))) {
  org_db <- tryCatch({
    switch(opt$organism,
      human = { suppressPackageStartupMessages(library(org.Hs.eg.db)); org.Hs.eg.db },
      mouse = { suppressPackageStartupMessages(library(org.Mm.eg.db)); org.Mm.eg.db },
      rat   = { suppressPackageStartupMessages(library(org.Rn.eg.db)); org.Rn.eg.db },
      NULL)
  }, error=function(e) NULL)
  if (!is.null(org_db)) {
    suppressPackageStartupMessages(library(AnnotationDbi))
    syms <- mapIds(org_db, keys=degs$gene, column="SYMBOL",
                   keytype="ENSEMBL", multiVals="first")
    degs$gene_symbol <- ifelse(is.na(syms), degs$gene, syms)
    message("Annotated ", sum(!is.na(syms)), "/", nrow(degs), " gene symbols")
  }
}

# ── EnhancedVolcano (NGS101 style + gene labels) ──────────────────────
label_col <- if ("gene_symbol" %in% colnames(degs)) "gene_symbol" else "gene"

n_up   <- sum(degs[[lfc_col]] >  opt$lfc & degs[[p_col]] < opt$alpha, na.rm=TRUE)
n_down <- sum(degs[[lfc_col]] < -opt$lfc & degs[[p_col]] < opt$alpha, na.rm=TRUE)

all_na_pval <- all(is.na(degs[[p_col]]))

if (all_na_pval) {
  # n=1 per group: no p-values — produce a fold-change rank plot instead
  message("n=1 per group — generating fold-change rank plots (p-values not available).")

  make_fc_plot <- function(n_label=0, path=NULL) {
    if (is.null(path)) return(invisible(NULL))
    d <- degs[order(-degs[[lfc_col]]), ]
    d$rank <- seq_len(nrow(d))
    d$dir  <- ifelse(d[[lfc_col]] >  opt$lfc, "Up",
               ifelse(d[[lfc_col]] < -opt$lfc, "Down", "NS"))
    n_up <- sum(d$dir == "Up")
    n_dn <- sum(d$dir == "Down")
    p <- ggplot(d, aes(rank, .data[[lfc_col]], color=dir)) +
      geom_point(size=0.5, alpha=0.75) +
      scale_color_manual(
        values = c(Up="firebrick", Down="steelblue", NS="grey60"),
        labels = c(Up   = paste0("Up |LFC|>", opt$lfc, " (n=", n_up, ")"),
                   Down = paste0("Down |LFC|>", opt$lfc, " (n=", n_dn, ")"),
                   NS   = "NS")) +
      geom_hline(yintercept=c(opt$lfc, -opt$lfc), linetype=2,
                 color="grey40", linewidth=0.4) +
      labs(x   = "Gene rank (ordered by log₂FC)",
           y   = paste0("log₂FC (", opt$treatment, " / ", opt$control, ")"),
           title    = paste("Fold Change:", opt$treatment, "vs", opt$control),
           subtitle = paste0("n=1 per group — p-values not available | ",
                             "Up: ", n_up, "  Down: ", n_dn),
           color = NULL) +
      theme_classic(base_size=12) +
      theme(legend.position="top",
            plot.subtitle=element_text(size=9, color="grey50"))
    if (n_label > 0) {
      top_genes <- d[d$dir != "NS", ]
      top_genes <- top_genes[order(-abs(top_genes[[lfc_col]])),
                             ][seq_len(min(n_label * 2, nrow(top_genes))), ]
      if (nrow(top_genes) > 0)
        p <- p + geom_text_repel(data=top_genes, aes(label=.data[[label_col]]),
                                 size=2.5, max.overlaps=20, color="grey20")
    }
    pdf(path, width=22/2.54, height=18/2.54)
    print(p)
    dev.off()
    message("Fold-change plot written: ", path)
  }

  make_fc_plot(n_label=0,         path=opt$out_volcano)
  make_fc_plot(n_label=opt$n_label, path=opt$out_volcano_labeled)

} else {
  pdf(opt$out_volcano, width=22/2.54, height=18/2.54)
  print(EnhancedVolcano(degs,
    lab              = rep("", nrow(degs)),
    x                = lfc_col,
    y                = p_col,
    title            = "Volcano plot",
    subtitle         = "EnhancedVolcano",
    caption          = paste0("total = ", nrow(degs), " variables"),
    pCutoff          = opt$alpha,
    FCcutoff         = opt$lfc,
    legendPosition   = "top",
    colAlpha         = 1,
    pointSize        = 0.3,
    drawConnectors   = FALSE
  ))
  dev.off()

  # ── Labeled volcano (top N up + N down annotated) ─────────────────────
  if (!is.null(opt$out_volcano_labeled)) {
    top_up <- degs[!is.na(degs[[p_col]]) & degs[[p_col]] < opt$alpha &
                     degs[[lfc_col]] >  opt$lfc, ]
    top_dn <- degs[!is.na(degs[[p_col]]) & degs[[p_col]] < opt$alpha &
                     degs[[lfc_col]] < -opt$lfc, ]
    top_up <- top_up[order(top_up[[p_col]]), ][seq_len(min(opt$n_label, nrow(top_up))), ]
    top_dn <- top_dn[order(top_dn[[p_col]]), ][seq_len(min(opt$n_label, nrow(top_dn))), ]
    select_genes <- c(top_up[[label_col]], top_dn[[label_col]])

    pdf(opt$out_volcano_labeled, width=22/2.54, height=18/2.54)
    print(EnhancedVolcano(degs,
      lab              = degs[[label_col]],
      selectLab        = select_genes,
      x                = lfc_col,
      y                = p_col,
      title            = "Volcano plot",
      subtitle         = "EnhancedVolcano",
      caption          = paste0("total = ", nrow(degs), " variables"),
      pCutoff          = opt$alpha,
      FCcutoff         = opt$lfc,
      labSize          = 3,
      legendPosition   = "top",
      colAlpha         = 1,
      pointSize        = 0.3,
      drawConnectors   = TRUE,
      max.overlaps     = Inf
    ))
    dev.off()
    message("Labeled volcano written: ", opt$out_volcano_labeled)
  }
}

# ── Heatmap (NGS101 style) ────────────────────────────────────────────
obj <- readRDS(opt$rds)
if (inherits(obj, "list") && "vsd" %in% names(obj)) {
  norm_mat  <- assay(obj$vsd)
  col_group <- obj$dds$group
} else if (inherits(obj, "list") && "norm_counts" %in% names(obj)) {
  # n=1 fallback: list(dge, norm_counts, res)
  norm_mat  <- obj$norm_counts
  col_group <- obj$dge$samples$group
} else if (inherits(obj, "DESeqDataSet")) {
  norm_mat  <- assay(obj)
  col_group <- obj$group
} else if (inherits(obj, "DGEList")) {
  norm_mat  <- edgeR::cpm(obj, log=TRUE)
  col_group <- obj$samples$group
} else {
  norm_mat  <- obj$v$E
  col_group <- obj$v$targets$group
}

# Significant DEGs ordered by FDR then |LFC|; fall back to top |LFC| when none pass
sig_degs <- degs[!is.na(degs[[p_col]]) &
                  degs[[p_col]] < opt$alpha &
                  abs(degs[[lfc_col]]) >= opt$lfc, ]
sig_degs <- sig_degs[order(sig_degs[[p_col]], -abs(sig_degs[[lfc_col]])), ]

if (nrow(sig_degs) == 0) {
  message("No significant DEGs — using top ", min(200, nrow(degs)),
          " genes by |LFC| for heatmaps.")
  sig_degs <- degs[order(-abs(degs[[lfc_col]])), ][seq_len(min(200, nrow(degs))), ]
}
sig_genes_all <- sig_degs$gene

full_mat <- norm_mat[rownames(norm_mat) %in% sig_genes_all, , drop=FALSE]

# Apply gene symbol row names
if ("gene_symbol" %in% colnames(degs)) {
  sym_map <- setNames(degs$gene_symbol, degs$gene)
  rownames(full_mat) <- ifelse(is.na(sym_map[rownames(full_mat)]),
                                rownames(full_mat), sym_map[rownames(full_mat)])
}

# Order columns by group (control first)
if (!is.null(col_group)) {
  col_order <- order(col_group)
  full_mat  <- full_mat[, col_order, drop=FALSE]
  col_group <- col_group[col_order]
  col_annot <- data.frame(Treatment=as.character(col_group),
                           row.names=colnames(full_mat))
} else {
  col_annot <- NULL
}

make_heatmap <- function(mat, path, title="Differential Expression Heatmap") {
  if (nrow(mat) < 2) {
    message("Fewer than 2 genes for heatmap: ", path, " — writing placeholder.")
    pdf(path); plot.new(); dev.off()
    return(invisible(NULL))
  }
  # z-score row scaling with 2 samples forces all genes to exactly ±0.707 (constant).
  # For 2 samples: mean-center rows instead so color reflects actual expression magnitude.
  hm_scale <- if (all_na_pval || ncol(mat) <= 2) "none" else "row"
  if (ncol(mat) <= 2) mat <- mat - rowMeans(mat)
  # Scale height: ~0.12 inch per row, min 4 inches, max 20 inches
  h <- max(4, min(20, nrow(mat) * 0.12 + 2))
  pdf(path, width=6, height=h)
  pheatmap(mat,
           scale          = hm_scale,
           cluster_rows   = TRUE,
           cluster_cols   = FALSE,
           show_rownames  = nrow(mat) <= 80,
           show_colnames  = TRUE,
           color          = colorRampPalette(c("navy","white","red"))(100),
           annotation_col = col_annot,
           main           = title)
  dev.off()
  message("Heatmap written (", nrow(mat), " genes): ", path)
}

# All significant DEGs
make_heatmap(full_mat, opt$out_heatmap,
             paste0("All DE Genes (n=", nrow(full_mat), ")"))

# Top 100 by expression variance (matches core: no p-value filter applied)
if (!is.null(opt$out_heatmap_top100)) {
  gene_vars  <- apply(norm_mat, 1, var)
  top_idx    <- order(-gene_vars)[seq_len(min(100, nrow(norm_mat)))]
  top100_mat <- norm_mat[top_idx, , drop=FALSE]
  if ("gene_symbol" %in% colnames(degs)) {
    sym_map <- setNames(degs$gene_symbol, degs$gene)
    rownames(top100_mat) <- ifelse(is.na(sym_map[rownames(top100_mat)]),
                                    rownames(top100_mat), sym_map[rownames(top100_mat)])
  }
  if (!is.null(col_group)) top100_mat <- top100_mat[, order(col_group), drop=FALSE]
  make_heatmap(top100_mat, opt$out_heatmap_top100,
               paste0("Top ", nrow(top100_mat), " Genes by Expression Variance"))
}

# Top 50 by expression variance
if (!is.null(opt$out_heatmap_top50)) {
  gene_vars  <- apply(norm_mat, 1, var)
  top_idx    <- order(-gene_vars)[seq_len(min(50, nrow(norm_mat)))]
  top50_mat  <- norm_mat[top_idx, , drop=FALSE]
  if ("gene_symbol" %in% colnames(degs)) {
    sym_map <- setNames(degs$gene_symbol, degs$gene)
    rownames(top50_mat) <- ifelse(is.na(sym_map[rownames(top50_mat)]),
                                   rownames(top50_mat), sym_map[rownames(top50_mat)])
  }
  if (!is.null(col_group)) top50_mat <- top50_mat[, order(col_group), drop=FALSE]
  make_heatmap(top50_mat, opt$out_heatmap_top50,
               paste0("Top ", nrow(top50_mat), " Genes by Expression Variance"))
}

# ── Separate up- and down-regulated heatmaps (top 100 by |LFC| or p-value) ──
make_directional_heatmap <- function(direction, path) {
  if (is.null(path)) return(invisible(NULL))
  if (all_na_pval) {
    # n=1: filter by LFC direction, sort by |LFC|
    dir_degs <- if (direction == "up")
      degs[!is.na(degs[[lfc_col]]) & degs[[lfc_col]] >  opt$lfc, ]
    else
      degs[!is.na(degs[[lfc_col]]) & degs[[lfc_col]] < -opt$lfc, ]
    dir_degs <- dir_degs[order(-abs(dir_degs[[lfc_col]])), ]
  } else {
    # Normal: filter by FDR + LFC direction, sort by smallest p-value
    dir_degs <- if (direction == "up")
      degs[!is.na(degs[[p_col]]) & degs[[p_col]] < opt$alpha &
             degs[[lfc_col]] >  opt$lfc, ]
    else
      degs[!is.na(degs[[p_col]]) & degs[[p_col]] < opt$alpha &
             degs[[lfc_col]] < -opt$lfc, ]
    dir_degs <- dir_degs[order(dir_degs[[p_col]]), ]
  }
  dir_degs <- dir_degs[seq_len(min(100, nrow(dir_degs))), ]
  dir_mat  <- norm_mat[rownames(norm_mat) %in% dir_degs$gene, , drop=FALSE]
  if ("gene_symbol" %in% colnames(degs)) {
    sym_map <- setNames(degs$gene_symbol, degs$gene)
    rownames(dir_mat) <- ifelse(is.na(sym_map[rownames(dir_mat)]),
                                rownames(dir_mat), sym_map[rownames(dir_mat)])
  }
  if (!is.null(col_group)) {
    col_ord  <- order(col_group)
    dir_mat  <- dir_mat[, col_ord, drop=FALSE]
  }
  dir_label <- if (direction == "up") "Up-regulated" else "Down-regulated"
  make_heatmap(dir_mat, path,
               paste0("Top ", nrow(dir_mat), " ", dir_label,
                      " Genes (smallest p-value)"))
}

make_directional_heatmap("up",   opt$out_heatmap_up)
make_directional_heatmap("down", opt$out_heatmap_down)

# ── Interactive HTML volcano ──────────────────────────────────────────
if (!is.null(opt$out_volcano_html)) {
  tryCatch({
    # Build plot_df — for n=1 (no p-values) use LFC direction only; y = |LFC|
    if (all_na_pval) {
      plot_df <- degs[!is.na(degs[[lfc_col]]), ]
      plot_df$sym     <- plot_df[[label_col]]
      plot_df$lfc_val <- round(plot_df[[lfc_col]], 4)
      plot_df$fdr_val <- NA_real_
      plot_df$nlp     <- round(abs(plot_df[[lfc_col]]), 4)  # |LFC| as y-axis
      plot_df$status  <- ifelse(plot_df[[lfc_col]] >  opt$lfc, "up",
                          ifelse(plot_df[[lfc_col]] < -opt$lfc, "down", "ns"))
      y_label   <- "|log₂FC|"
      has_fdr   <- "false"
    } else {
      plot_df <- degs[!is.na(degs[[lfc_col]]) & !is.na(degs[[p_col]]), ]
      plot_df$sym     <- plot_df[[label_col]]
      plot_df$lfc_val <- round(plot_df[[lfc_col]], 4)
      plot_df$fdr_val <- plot_df[[p_col]]
      plot_df$nlp     <- round(-log10(pmax(plot_df[[p_col]], 1e-300)), 4)
      plot_df$status  <- ifelse(
        plot_df[[lfc_col]] >  opt$lfc & plot_df[[p_col]] < opt$alpha, "up",
        ifelse(plot_df[[lfc_col]] < -opt$lfc & plot_df[[p_col]] < opt$alpha, "down", "ns"))
      y_label   <- "-log₁₀(FDR)"
      has_fdr   <- "true"
    }

    esc <- function(x) gsub('"', '\\\\"', gsub("\\\\", "\\\\\\\\", as.character(x)))
    rows <- paste(apply(plot_df[, c("sym","lfc_val","nlp","fdr_val","status")], 1, function(r) {
      p_str <- if (is.na(r[4])) "null" else r[4]
      sprintf('{"g":"%s","x":%s,"y":%s,"p":%s,"s":"%s"}',
              esc(r[1]), r[2], r[3], p_str, r[5])
    }), collapse=",")

    fdr_line  <- round(-log10(opt$alpha), 4)
    title_str <- paste0(opt$treatment, " vs ", opt$control)

    html <- sprintf('<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>%s Volcano</title>
<style>
*{box-sizing:border-box;margin:0;padding:0}
body{font-family:system-ui,sans-serif;background:#f8f8f8;color:#222;padding:16px}
h2{font-size:1.1rem;margin-bottom:8px;text-align:center}
#controls{display:flex;gap:10px;align-items:center;flex-wrap:wrap;margin-bottom:8px}
#search{padding:4px 8px;border:1px solid #ccc;border-radius:4px;font-size:.9rem;width:180px}
.legend-dot{display:inline-block;width:10px;height:10px;border-radius:50%%;margin-right:4px}
.leg{display:flex;align-items:center;font-size:.85rem;gap:4px;cursor:pointer;
  padding:3px 7px;border-radius:4px;border:1.5px solid transparent;user-select:none;transition:opacity .15s}
.leg:hover{background:rgba(0,0,0,.06)}
.leg.hidden{opacity:.35;text-decoration:line-through}
#tooltip{position:fixed;background:rgba(0,0,0,.82);color:#fff;padding:6px 10px;border-radius:5px;
  font-size:.8rem;pointer-events:none;display:none;white-space:nowrap;z-index:999}
svg text{font-family:system-ui,sans-serif}
.axis path,.axis line{stroke:#999}
.axis text{font-size:11px;fill:#444}
#chart{background:#fff;border:1px solid #e0e0e0;border-radius:6px;display:block;margin:0 auto}
</style>
</head>
<body>
<h2>%s</h2>
<div id="controls">
  <input id="search" placeholder="Search gene…" oninput="highlight(this.value)">
  <span class="leg" data-cat="up" onclick="toggleCat(this,this.dataset.cat)" title="Click to show/hide">
    <span class="legend-dot" style="background:#e53935"></span><span id="leg-up"></span></span>
  <span class="leg" data-cat="down" onclick="toggleCat(this,this.dataset.cat)" title="Click to show/hide">
    <span class="legend-dot" style="background:#1e88e5"></span><span id="leg-dn"></span></span>
  <span class="leg" data-cat="ns" onclick="toggleCat(this,this.dataset.cat)" title="Click to show/hide">
    <span class="legend-dot" style="background:#888"></span>Non-significant</span>
  <span style="font-size:.8rem;color:#666;margin-left:auto">Click legend to filter · Scroll to zoom · Drag to pan</span>
</div>
<svg id="chart"></svg>
<div id="tooltip"></div>
<script src="https://cdnjs.cloudflare.com/ajax/libs/d3/7.9.0/d3.min.js"></script>
<script>
const RAW = [%s];
const FDR_CUT = %s, LFC_CUT = %s, HAS_FDR = %s;
const Y_LABEL = "%s";

const W = Math.min(window.innerWidth - 40, 900),
      H = Math.round(W * 0.72),
      M = {top:20, right:30, bottom:50, left:58};
const IW = W - M.left - M.right, IH = H - M.top - M.bottom;

const svg = d3.select("#chart").attr("width",W).attr("height",H);
const g = svg.append("g").attr("transform",`translate(${M.left},${M.top})`);
const clip = svg.append("defs").append("clipPath").attr("id","clip")
  .append("rect").attr("width",IW).attr("height",IH);
const plot = g.append("g").attr("clip-path","url(#clip)");

const xExt = d3.extent(RAW, d=>d.x), yExt = [0, d3.max(RAW, d=>d.y)*1.05];
const xPad = (xExt[1]-xExt[0])*0.05;
let x = d3.scaleLinear().domain([xExt[0]-xPad, xExt[1]+xPad]).range([0,IW]);
let y = d3.scaleLinear().domain(yExt).range([IH,0]);

const xAxis = d3.axisBottom(x).ticks(8);
const yAxis = d3.axisLeft(y).ticks(8);
const gx = g.append("g").attr("class","axis").attr("transform",`translate(0,${IH})`).call(xAxis);
const gy = g.append("g").attr("class","axis").call(yAxis);

g.append("text").attr("x",IW/2).attr("y",IH+40).attr("text-anchor","middle")
  .attr("font-size","13px").text("log₂ Fold Change");
g.append("text").attr("transform","rotate(-90)").attr("x",-IH/2).attr("y",-44)
  .attr("text-anchor","middle").attr("font-size","13px").text(Y_LABEL);

const refLine = plot.append("g").attr("class","refs");
const drawRefs = (cx,cy) => {
  refLine.selectAll("*").remove();
  [LFC_CUT,-LFC_CUT].forEach(v => refLine.append("line")
    .attr("x1",cx(v)).attr("x2",cx(v)).attr("y1",0).attr("y2",IH)
    .attr("stroke","#aaa").attr("stroke-dasharray","4,3").attr("stroke-width",1));
  if (HAS_FDR) refLine.append("line")
    .attr("x1",0).attr("x2",IW).attr("y1",cy(FDR_CUT)).attr("y2",cy(FDR_CUT))
    .attr("stroke","#aaa").attr("stroke-dasharray","4,3").attr("stroke-width",1);
};
drawRefs(x,y);

const COL = {up:"#e53935", down:"#1e88e5", ns:"#888"};
const fill = d => COL[d.s] || "#888";
const op   = d => d.s==="ns" ? 0.35 : 0.85;

const circles = plot.selectAll("circle").data(RAW).join("circle")
  .attr("cx", d=>x(d.x)).attr("cy", d=>y(d.y))
  .attr("r", d=>d.s==="ns"?2:3)
  .attr("fill", fill).attr("fill-opacity", op)
  .attr("stroke","none");

// Top 10 up + 10 down labels — sort by p-value if available, else by |LFC|
const sortFn = HAS_FDR ? (a,b)=>a.p-b.p : (a,b)=>Math.abs(b.x)-Math.abs(a.x);
const topUp = RAW.filter(d=>d.s==="up").sort(sortFn).slice(0,10);
const topDn = RAW.filter(d=>d.s==="down").sort(sortFn).slice(0,10);
const labeled = [...topUp, ...topDn];
const labels = plot.selectAll("text.gene").data(labeled).join("text")
  .attr("class","gene")
  .attr("x", d=>x(d.x)+5).attr("y", d=>y(d.y)+4)
  .attr("font-size","8px").attr("fill", fill).text(d=>d.g)
  .attr("pointer-events","none");

// Legend counts
const nUp = RAW.filter(d=>d.s==="up").length;
const nDn = RAW.filter(d=>d.s==="down").length;
document.getElementById("leg-up").textContent = `Upregulated (n=${nUp})`;
document.getElementById("leg-dn").textContent = `Downregulated (n=${nDn})`;

// Tooltip
const tip = d3.select("#tooltip");
circles.on("mouseover", (e,d) => {
  const pStr = HAS_FDR && d.p !== null
    ? `<br>FDR: ${d.p < 0.001 ? d.p.toExponential(2) : d.p.toFixed(4)}`
    : "";
  tip.style("display","block")
    .html(`<b>${d.g}</b><br>LFC: ${d.x}<br>${Y_LABEL}: ${d.y}${pStr}`);
}).on("mousemove", e => {
  tip.style("left",(e.clientX+14)+"px").style("top",(e.clientY-30)+"px");
}).on("mouseout", () => tip.style("display","none"));

// Zoom
const zoom = d3.zoom().scaleExtent([0.4,40]).on("zoom", ({transform}) => {
  const nx = transform.rescaleX(x), ny = transform.rescaleY(y);
  gx.call(xAxis.scale(nx)); gy.call(yAxis.scale(ny));
  circles.attr("cx", d=>nx(d.x)).attr("cy", d=>ny(d.y));
  labels.attr("x", d=>nx(d.x)+5).attr("y", d=>ny(d.y)+4);
  drawRefs(nx, ny);
});
svg.call(zoom);

// Category filter
const hidden = new Set();
function toggleCat(el, cat) {
  if (hidden.has(cat)) { hidden.delete(cat); el.classList.remove("hidden"); }
  else                 { hidden.add(cat);    el.classList.add("hidden"); }
  applyVisibility();
}
function applyVisibility() {
  circles.attr("display", d => hidden.has(d.s) ? "none" : null);
  labels.attr("display",  d => hidden.has(d.s) ? "none" : null);
}

// Search
function highlight(q) {
  const lq = q.toLowerCase().trim();
  circles.attr("stroke", d => lq && d.g.toLowerCase().includes(lq) ? "#000" : "none")
         .attr("stroke-width", 1.5)
         .attr("r", d => (lq && d.g.toLowerCase().includes(lq)) ? 5 : (d.s==="ns"?2:3));
}
</script>
</body>
</html>', title_str, title_str, rows, fdr_line, opt$lfc, has_fdr, y_label)

    writeLines(html, opt$out_volcano_html)
    message("Interactive plot written: ", opt$out_volcano_html)
  }, error=function(e) message("Interactive plot skipped: ", conditionMessage(e)))
}

message("Figures written.")
