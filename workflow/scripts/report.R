#!/usr/bin/env Rscript
# report.R — HTML pipeline report generator (Part 20)
# Produces a single self-contained HTML with embedded plots and tables.
# Requires: magick, data.table, base64enc (all available in rnaseq_r_env)
# Usage: Called by Snakemake rule html_report

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(base64enc)
})

opt_list <- list(
  make_option("--label",          type="character"),
  make_option("--treatment",      type="character"),
  make_option("--control",        type="character"),
  make_option("--degs",           type="character"),
  make_option("--multiqc",        type="character", default=""),
  make_option("--pca",            type="character", default=""),
  make_option("--ma",             type="character", default=""),
  make_option("--sample_dist",    type="character", default=""),
  make_option("--volcano",        type="character", default=""),
  make_option("--volcano_html",   type="character", default=""),
  make_option("--heatmap",        type="character", default=""),
  make_option("--heatmap_top100", type="character", default=""),
  make_option("--heatmap_top50",  type="character", default=""),
  make_option("--pathways_dir",   type="character", default=""),
  make_option("--output",         type="character")
)
opt <- parse_args(OptionParser(option_list=opt_list))

# ── Helpers ─────────────────────────────────────────────────────────────

pdf_to_b64_png <- function(path, dpi=150) {
  if (!nzchar(path) || !file.exists(path)) return(NULL)
  tryCatch({
    img <- magick::image_read_pdf(path, density=dpi)
    img <- magick::image_convert(img, "png")
    tmp <- tempfile(fileext=".png")
    on.exit(unlink(tmp), add=TRUE)
    magick::image_write(img, tmp)
    paste0("data:image/png;base64,", base64enc::base64encode(tmp))
  }, error=function(e) { message("Could not convert ", path, ": ", e$message); NULL })
}

img_tag <- function(b64, alt="", extra_style="max-width:100%;border-radius:6px;") {
  if (is.null(b64)) return('<p style="color:#888;font-style:italic">Plot not available</p>')
  sprintf('<img src="%s" alt="%s" style="%s"/>', b64, alt, extra_style)
}

read_tsv_html <- function(path, n=20) {
  if (!nzchar(path) || !file.exists(path)) return('<p style="color:#888">Table not available</p>')
  tryCatch({
    dt <- fread(path, data.table=FALSE, nrows=n)
    nums <- sapply(dt, is.numeric)
    dt[, nums] <- lapply(dt[, nums, drop=FALSE], function(x) round(x, 4))
    rows <- apply(dt, 1, function(r) {
      paste0("<tr>", paste0("<td>", r, "</td>", collapse=""), "</tr>")
    })
    hdr <- paste0("<th>", colnames(dt), "</th>", collapse="")
    paste0('<div style="overflow-x:auto"><table class="dtable">',
           '<thead><tr>', hdr, '</tr></thead><tbody>',
           paste(rows, collapse=""), '</tbody></table></div>')
  }, error=function(e) paste0('<p style="color:#888">Could not read table: ', e$message, '</p>'))
}

section <- function(id, title, icon, content) {
  sprintf('
  <section id="%s">
    <h2>%s %s</h2>
    %s
  </section>', id, icon, title, content)
}

subsection <- function(title, content) {
  sprintf('<div class="subsec"><h3>%s</h3>%s</div>', title, content)
}

fig_grid <- function(...) {
  items <- list(...)
  cells <- sapply(items, function(x) sprintf('<div class="fig-cell">%s</div>', x))
  paste0('<div class="fig-grid">', paste(cells, collapse=""), '</div>')
}

# ── Load magick only if available ────────────────────────────────────────
has_magick <- requireNamespace("magick", quietly=TRUE)
if (!has_magick) {
  message("WARNING: magick not installed — PDF plots will not be embedded. ",
          "Install with: install.packages('magick')")
  pdf_to_b64_png <- function(path, dpi=150) NULL
}

# ── Pre-convert plots ─────────────────────────────────────────────────────
message("Converting plots to PNG...")
b64 <- list(
  pca         = pdf_to_b64_png(opt$pca),
  ma          = pdf_to_b64_png(opt$ma),
  sample_dist = pdf_to_b64_png(opt$sample_dist),
  volcano     = pdf_to_b64_png(opt$volcano),
  heatmap     = pdf_to_b64_png(opt$heatmap),
  hm100       = pdf_to_b64_png(opt$heatmap_top100),
  hm50        = pdf_to_b64_png(opt$heatmap_top50)
)

# Pathway plots from pathways_dir
pw_plots <- list()
if (nzchar(opt$pathways_dir) && dir.exists(opt$pathways_dir)) {
  pw_names <- c(
    go_dot    = "GO-ALL/output-GO-dotplot.pdf",
    go_bar    = "GO-ALL/output-GO-barplot.pdf",
    go_cnet   = "GO-ALL/output-GO-cnetplot.pdf",
    go_emap   = "GO-ALL/output-GO-emapplot.pdf",
    go_cc_dot = "GO-ALL/output-CellularComponent-dotplot.pdf",
    go_cc_bar = "GO-ALL/output-CellularComponent-barplot.pdf",
    go_mf_dot = "GO-ALL/output-MolecularFunction-dotplot.pdf",
    go_mf_bar = "GO-ALL/output-MolecularFunction-barplot.pdf",
    go_up_dot = "GO-UP/output-GO-dotplot.pdf",
    go_dn_dot = "GO-DOWN/output-GO-dotplot.pdf",
    kegg_dot  = "GO-ALL/output-kegg-dotplot.pdf",
    gsea_dot  = paste0(opt$label, "_GSEA_dotplot.pdf"),
    gsea_rdg  = paste0(opt$label, "_GSEA_ridge.pdf")
  )
  for (nm in names(pw_names)) {
    fpath <- file.path(opt$pathways_dir, pw_names[[nm]])
    pw_plots[[nm]] <- pdf_to_b64_png(fpath)
  }
}

# ── DEG summary stats ─────────────────────────────────────────────────────
degs <- tryCatch(fread(opt$degs, data.table=FALSE), error=function(e) NULL)
lfc_col <- if (!is.null(degs)) intersect(c("log2FC","log2FoldChange","logFC"), colnames(degs))[1] else NULL
p_col   <- if (!is.null(degs)) intersect(c("FDR","adj.P.Val","padj"), colnames(degs))[1] else NULL

deg_stats <- list(total=0, up=0, down=0)
deg_table_html <- '<p style="color:#888">DEG table not available</p>'
if (!is.null(degs) && !is.null(p_col) && !is.null(lfc_col)) {
  sig <- degs[!is.na(degs[[p_col]]) & degs[[p_col]] < 0.05, ]
  deg_stats$total <- nrow(sig)
  deg_stats$up    <- sum(sig[[lfc_col]] > 1, na.rm=TRUE)
  deg_stats$down  <- sum(sig[[lfc_col]] < -1, na.rm=TRUE)

  sym_col <- if ("gene_symbol" %in% colnames(degs)) "gene_symbol" else "gene"
  show_cols <- intersect(c(sym_col, lfc_col, p_col, "baseMean","lfcSE","stat","pvalue"),
                          colnames(degs))
  top_degs <- sig[order(sig[[p_col]]), show_cols, drop=FALSE][seq_len(min(20, nrow(sig))), ]
  nums <- sapply(top_degs, is.numeric)
  top_degs[, nums] <- lapply(top_degs[, nums, drop=FALSE], function(x) signif(x, 4))
  rows_html <- apply(top_degs, 1, function(r) paste0("<tr><td>", paste(r, collapse="</td><td>"), "</td></tr>"))
  hdr_html  <- paste0("<th>", colnames(top_degs), "</th>", collapse="")
  deg_table_html <- paste0(
    '<div style="overflow-x:auto"><table class="dtable">',
    '<thead><tr>', hdr_html, '</tr></thead><tbody>',
    paste(rows_html, collapse=""), '</tbody></table></div>'
  )
}

# ── Volcano HTML embed ─────────────────────────────────────────────────────
volcano_iframe <- ""
if (nzchar(opt$volcano_html) && file.exists(opt$volcano_html)) {
  vcontent <- tryCatch(paste(readLines(opt$volcano_html), collapse="\n"), error=function(e) "")
  if (nzchar(vcontent)) {
    vb64 <- base64enc::base64encode(charToRaw(vcontent))
    volcano_iframe <- paste0(
      '<iframe src="data:text/html;base64,', vb64,
      '" style="width:100%;height:600px;border:none;border-radius:8px;" title="Interactive Volcano Plot"></iframe>')
  }
}

# ── Build HTML ────────────────────────────────────────────────────────────
now_str <- format(Sys.time(), "%Y-%m-%d %H:%M")

# Helper: inline value safely (no sprintf, avoids % interpretation issues)
V <- function(x) as.character(x)

volcano_section <- if (nzchar(volcano_iframe)) {
  paste0('<div class="subsec"><h3>Interactive Volcano (D3)</h3>', volcano_iframe, '</div>')
} else {
  ""
}

html <- paste0(
'<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1.0"/>
<title>RNA-Seq Report: ', V(opt$label), '</title>
<style>
:root {
  --sidebar-bg: #1a3a5c;
  --sidebar-hover: #2a5a8c;
  --sidebar-active: #4a9fd4;
  --accent: #4a9fd4;
  --accent2: #e74c3c;
  --bg: #f7f9fc;
  --card-bg: #ffffff;
  --text: #2c3e50;
  --muted: #6c757d;
  --border: #dee2e6;
  font-size: 15px;
}
* { box-sizing: border-box; margin: 0; padding: 0; }
body { display: flex; min-height: 100vh; font-family: "Segoe UI", system-ui, sans-serif;
       background: var(--bg); color: var(--text); }

/* Sidebar */
#sidebar {
  width: 240px; flex-shrink: 0; background: var(--sidebar-bg); color: #fff;
  position: fixed; top: 0; left: 0; height: 100vh; overflow-y: auto;
  display: flex; flex-direction: column; z-index: 100;
}
#sidebar .logo { padding: 24px 20px 16px; border-bottom: 1px solid rgba(255,255,255,0.1); }
#sidebar .logo h1 { font-size: 1.1rem; font-weight: 700; line-height: 1.3; }
#sidebar .logo p  { font-size: 0.75rem; opacity: 0.7; margin-top: 4px; }
#sidebar nav { flex: 1; padding: 12px 0; }
#sidebar nav a {
  display: block; padding: 9px 20px; color: rgba(255,255,255,0.85);
  text-decoration: none; font-size: 0.88rem; border-left: 3px solid transparent;
  transition: background 0.15s, border-color 0.15s;
}
#sidebar nav a:hover   { background: var(--sidebar-hover); border-left-color: var(--accent); }
#sidebar nav a.active  { background: var(--sidebar-hover); border-left-color: var(--accent); color:#fff; font-weight:600; }
#sidebar nav .nav-grp  { padding: 12px 20px 4px; font-size: 0.7rem; text-transform: uppercase;
                          letter-spacing: 0.1em; opacity: 0.5; }
#sidebar .meta { padding: 14px 20px; font-size: 0.72rem; opacity: 0.55;
                 border-top: 1px solid rgba(255,255,255,0.1); }

/* Main */
#main { margin-left: 240px; padding: 0 32px 48px; max-width: 1200px; width: 100%; }
section { padding-top: 40px; }
section + section { border-top: 1px solid var(--border); margin-top: 12px; }
h2 { font-size: 1.5rem; font-weight: 700; color: var(--sidebar-bg); margin-bottom: 20px;
     padding-bottom: 8px; border-bottom: 2px solid var(--accent); }
h3 { font-size: 1.05rem; font-weight: 600; color: var(--text); margin-bottom: 12px; margin-top: 4px; }

/* Cards / subsections */
.subsec { background: var(--card-bg); border: 1px solid var(--border); border-radius: 8px;
          padding: 20px; margin-bottom: 20px; box-shadow: 0 1px 4px rgba(0,0,0,0.05); }

/* Stat tiles */
.stat-row { display: flex; gap: 16px; flex-wrap: wrap; margin-bottom: 24px; }
.stat-tile { flex: 1; min-width: 130px; background: var(--card-bg); border: 1px solid var(--border);
             border-radius: 10px; padding: 18px 16px; text-align: center;
             box-shadow: 0 1px 4px rgba(0,0,0,0.06); }
.stat-tile .num { font-size: 2.2rem; font-weight: 800; color: var(--sidebar-bg); line-height: 1; }
.stat-tile .num.up   { color: #c0392b; }
.stat-tile .num.down { color: #2980b9; }
.stat-tile .lbl { font-size: 0.8rem; color: var(--muted); margin-top: 4px; }

/* Figure grid */
.fig-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
            gap: 20px; margin-top: 8px; }
.fig-cell img { max-width: 100%; border: 1px solid var(--border); border-radius: 6px; }
.fig-cell p.cap { font-size: 0.78rem; color: var(--muted); margin-top: 6px; text-align: center; }

/* Tables */
.dtable { width: 100%; border-collapse: collapse; font-size: 0.82rem; }
.dtable th { background: var(--sidebar-bg); color: #fff; padding: 8px 10px;
             text-align: left; font-weight: 600; white-space: nowrap; }
.dtable td { padding: 6px 10px; border-bottom: 1px solid var(--border); }
.dtable tr:nth-child(even) td { background: #f2f6fa; }
.dtable tr:hover td { background: #e8f0f8; }

/* Banner */
.banner { background: linear-gradient(135deg, var(--sidebar-bg) 0%, #2a6496 100%);
          color: #fff; padding: 28px 32px; margin: 0 -32px 0; }
.banner h1 { font-size: 1.8rem; font-weight: 800; }
.banner p  { opacity: 0.8; margin-top: 6px; font-size: 0.95rem; }
.badge { display: inline-block; padding: 3px 10px; border-radius: 20px;
         font-size: 0.75rem; font-weight: 600; margin-right: 6px; }
.badge-trt { background: rgba(231,76,60,0.25); color: #ffc; border: 1px solid rgba(231,76,60,0.4); }
.badge-ctl { background: rgba(74,159,212,0.25); color: #cef; border: 1px solid rgba(74,159,212,0.4); }

@media (max-width: 768px) {
  #sidebar { width: 200px; }
  #main { margin-left: 200px; padding: 0 16px 40px; }
  .fig-grid { grid-template-columns: 1fr; }
}
@media (max-width: 560px) {
  #sidebar { display: none; }
  #main { margin-left: 0; }
}
</style>
</head>
<body>

<nav id="sidebar">
  <div class="logo">
    <h1>RNA-Seq Pipeline</h1>
    <p>Analysis Report</p>
  </div>
  <nav>
    <div class="nav-grp">Overview</div>
    <a href="#overview">&#128200; Summary</a>
    <div class="nav-grp">QC</div>
    <a href="#qc-section">&#9989; Quality Control</a>
    <div class="nav-grp">Differential Expression</div>
    <a href="#deg-section">&#128202; DEG Analysis</a>
    <a href="#volcano-section">&#127776; Volcano Plot</a>
    <a href="#heatmap-section">&#127863; Heatmaps</a>
    <div class="nav-grp">Pathway Enrichment</div>
    <a href="#go-bp-section">&#127758; GO Biological Process</a>
    <a href="#go-cc-mf-section">&#128084; GO CC / MF</a>
    <a href="#go-dir-section">&#8593;&#8595; Up/Down Pathways</a>
    <a href="#kegg-section">&#128202; KEGG</a>
    <a href="#gsea-section">&#128288; GSEA</a>
    <div class="nav-grp">Data</div>
    <a href="#deg-table-section">&#128203; DEG Table</a>
    <a href="#methods-section">&#128196; Methods</a>
  </nav>
  <div class="meta">Generated ', V(now_str), '<br/>NGS101 Pipeline</div>
</nav>

<div id="main">

<div class="banner">
  <h1>', V(opt$label), '</h1>
  <p>
    <span class="badge badge-trt">Treatment: ', V(opt$treatment), '</span>
    <span class="badge badge-ctl">Control: ', V(opt$control), '</span>
  </p>
  <p style="margin-top:10px;opacity:0.7;font-size:0.85rem;">
    RNA-Seq differential expression and pathway enrichment analysis
  </p>
</div>

<!-- OVERVIEW -->
<section id="overview">
  <h2>&#128200; Summary</h2>
  <div class="stat-row">
    <div class="stat-tile">
      <div class="num">', V(deg_stats$total), '</div>
      <div class="lbl">Total DE genes<br/>(FDR &lt; 0.05)</div>
    </div>
    <div class="stat-tile">
      <div class="num up">', V(deg_stats$up), '</div>
      <div class="lbl">Upregulated<br/>(log2FC &gt; 1)</div>
    </div>
    <div class="stat-tile">
      <div class="num down">', V(deg_stats$down), '</div>
      <div class="lbl">Downregulated<br/>(log2FC &lt; -1)</div>
    </div>
  </div>
  <div class="subsec">
    <h3>Pipeline Overview</h3>
    <p style="line-height:1.7;color:var(--muted)">
      This report summarises the NGS101 RNA-seq pipeline run for the comparison
      <strong>', V(opt$treatment), ' vs ', V(opt$control), '</strong>.
      Reads were trimmed with Trim Galore, aligned to the reference genome with STAR,
      and counts were quantified with featureCounts.
      Differential expression was performed with DESeq2 (Wald test, BH-adjusted FDR).
      Pathway enrichment used clusterProfiler (GO ORA, KEGG ORA, and GSEA on GO BP terms).
    </p>
  </div>
</section>

<!-- QC -->
<section id="qc-section">
  <h2>&#9989; Quality Control</h2>
  <div class="subsec">
    <h3>PCA Plot</h3>
    ', img_tag(b64$pca, "PCA"), '
  </div>
  <div class="subsec">
    <h3>Sample Distance Matrix</h3>
    ', img_tag(b64$sample_dist, "Sample Distance"), '
  </div>
  <div class="subsec">
    <h3>MA Plot</h3>
    ', img_tag(b64$ma, "MA Plot"), '
  </div>
</section>

<!-- DEG -->
<section id="deg-section">
  <h2>&#128202; Differential Expression Analysis</h2>
  <div class="subsec">
    <h3>Overview</h3>
    <p style="color:var(--muted);line-height:1.7">
      DESeq2 was used with a Wald test. Log2-fold-change was shrunk with the
      <code>apeglm</code> estimator. Genes with mean count &lt; 10 across all
      samples were pre-filtered. Significance threshold: FDR &lt; 0.05,
      |log2FC| &gt; 1.
    </p>
  </div>
</section>

<!-- VOLCANO -->
<section id="volcano-section">
  <h2>&#127776; Volcano Plot</h2>
  <div class="subsec">
    <h3>Static Volcano (EnhancedVolcano)</h3>
    ', img_tag(b64$volcano, "Volcano"), '
    <p class="cap" style="font-size:0.8rem;color:var(--muted);margin-top:8px">
      Top 10 most significant upregulated and downregulated genes are labelled.
      Red: |log2FC| &gt; 1 &amp; FDR &lt; 0.05. Blue: FDR &lt; 0.05 only. Grey: not significant.
    </p>
  </div>
  ', volcano_section, '
</section>

<!-- HEATMAPS -->
<section id="heatmap-section">
  <h2>&#127863; Heatmaps</h2>
  <div class="subsec">
    <h3>All DE Genes</h3>
    ', img_tag(b64$heatmap, "All DE genes heatmap"), '
  </div>
  <div class="subsec">
    <h3>Top 100 DE Genes</h3>
    ', img_tag(b64$hm100, "Top 100 heatmap"), '
  </div>
  <div class="subsec">
    <h3>Top 50 DE Genes</h3>
    ', img_tag(b64$hm50, "Top 50 heatmap"), '
  </div>
</section>

<!-- GO BP -->
<section id="go-bp-section">
  <h2>&#127758; GO Biological Process</h2>
  <div class="fig-grid">
    <div class="fig-cell">', img_tag(pw_plots[["go_dot"]], "GO BP dotplot"), '<p class="cap">Dotplot</p></div>
    <div class="fig-cell">', img_tag(pw_plots[["go_bar"]], "GO BP barplot"), '<p class="cap">Barplot</p></div>
  </div>
  <div class="fig-grid">
    <div class="fig-cell">', img_tag(pw_plots[["go_cnet"]], "GO BP cnetplot"), '<p class="cap">Network plot (cnetplot)</p></div>
    <div class="fig-cell">', img_tag(pw_plots[["go_emap"]], "GO BP emapplot"), '<p class="cap">Enrichment map (emapplot)</p></div>
  </div>
</section>

<!-- GO CC / MF -->
<section id="go-cc-mf-section">
  <h2>&#128084; GO Cellular Component &amp; Molecular Function</h2>
  <div class="fig-grid">
    <div class="fig-cell">', img_tag(pw_plots[["go_cc_dot"]], "GO CC dotplot"), '<p class="cap">CC Dotplot</p></div>
    <div class="fig-cell">', img_tag(pw_plots[["go_cc_bar"]], "GO CC barplot"), '<p class="cap">CC Barplot</p></div>
  </div>
  <div class="fig-grid">
    <div class="fig-cell">', img_tag(pw_plots[["go_mf_dot"]], "GO MF dotplot"), '<p class="cap">MF Dotplot</p></div>
    <div class="fig-cell">', img_tag(pw_plots[["go_mf_bar"]], "GO MF barplot"), '<p class="cap">MF Barplot</p></div>
  </div>
</section>

<!-- UP / DOWN -->
<section id="go-dir-section">
  <h2>&#8593;&#8595; Directional GO Analysis</h2>
  <div class="fig-grid">
    <div class="fig-cell">', img_tag(pw_plots[["go_up_dot"]], "Upregulated GO dotplot"), '<p class="cap">Upregulated genes — GO BP dotplot</p></div>
    <div class="fig-cell">', img_tag(pw_plots[["go_dn_dot"]], "Downregulated GO dotplot"), '<p class="cap">Downregulated genes — GO BP dotplot</p></div>
  </div>
</section>

<!-- KEGG -->
<section id="kegg-section">
  <h2>&#128202; KEGG Pathway Enrichment</h2>
  <div class="subsec">', img_tag(pw_plots[["kegg_dot"]], "KEGG dotplot"), '</div>
</section>

<!-- GSEA -->
<section id="gsea-section">
  <h2>&#128288; Gene Set Enrichment Analysis (GSEA)</h2>
  <div class="fig-grid">
    <div class="fig-cell">', img_tag(pw_plots[["gsea_dot"]], "GSEA dotplot"), '<p class="cap">GSEA dotplot (split by direction)</p></div>
    <div class="fig-cell">', img_tag(pw_plots[["gsea_rdg"]], "GSEA ridgeplot"), '<p class="cap">GSEA ridgeplot</p></div>
  </div>
</section>

<!-- DEG TABLE -->
<section id="deg-table-section">
  <h2>&#128203; Top DE Genes</h2>
  <div class="subsec">
    <h3>Top 20 significant genes (sorted by FDR)</h3>
    ', deg_table_html, '
  </div>
</section>

<!-- METHODS -->
<section id="methods-section">
  <h2>&#128196; Methods</h2>
  <div class="subsec">
    <h3>Analysis Pipeline</h3>
    <p style="line-height:1.8;color:var(--muted)">
      <strong>Quality control:</strong> FastQC v0.12 + MultiQC v1.22.<br/>
      <strong>Trimming:</strong> Trim Galore v0.6 (Cutadapt backend).<br/>
      <strong>Alignment:</strong> STAR v2.7 (two-pass mode) against the reference genome.<br/>
      <strong>Quantification:</strong> featureCounts (Subread v2.0) on gene-level annotation.<br/>
      <strong>Differential expression:</strong> DESeq2 v1.42 (Wald test; apeglm LFC shrinkage).<br/>
      <strong>Gene annotation:</strong> Bioconductor org.Mm.eg.db / org.Hs.eg.db.<br/>
      <strong>Pathway enrichment:</strong> clusterProfiler v4 — GO ORA (BP/CC/MF), KEGG ORA, and GSEA on ranked log2FC.<br/>
      <strong>Visualisation:</strong> ggplot2, EnhancedVolcano, pheatmap, enrichplot, D3.js v7.
    </p>
  </div>
</section>

</div><!-- /main -->

<script>
const sections = document.querySelectorAll("section[id]");
const links    = document.querySelectorAll("#sidebar nav a");
const observer = new IntersectionObserver(entries => {
  entries.forEach(e => {
    if (e.isIntersecting) {
      links.forEach(l => l.classList.remove("active"));
      const active = document.querySelector(`#sidebar nav a[href="#${e.target.id}"]`);
      if (active) active.classList.add("active");
    }
  });
}, { threshold: 0.2 });
sections.forEach(s => observer.observe(s));
</script>
</body>
</html>'
)

writeLines(html, opt$output)
message("Report written to: ", opt$output)
