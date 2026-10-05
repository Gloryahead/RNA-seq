#!/usr/bin/env python3
"""
mapt_summary_plots.py — Gene-wide summary bar plots from mapt_complete_analysis.py TSV outputs.

Reads the TSV files produced by mapt_complete_analysis.py and generates one
summary value per sample per panel — no per-intron breakdown.

PANELS
  B  Gene-wide % unspliced          — from _aggregate.tsv
  C  Mean IR ratio across introns   — from _ir_ratio.tsv (weighted by coverage)
  D  PSI per alternative exon       — from _psi.tsv (5 exons × samples)
  E  3R / 4R ratio + RPM            — from _3r4r.tsv
  F  Mean RPKM of constitutive exons— from _exon_rpkm.tsv (exons 1,4,5,7,9,11,12,13)

USAGE
  python mapt_summary_plots.py --prefix MAPT_set1
  python mapt_summary_plots.py --prefix MAPT_set1 --out MAPT_set1_summary
  python mapt_summary_plots.py --prefix MAPT_set1 MAPT_set2 --labels Set1 Set2
"""
import argparse, math, sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# ── Constitutive exons (not alternatively spliced) ────────────────────────────
CONSTITUTIVE = {"exon 1","exon 4","exon 5","exon 7","exon 9",
                "exon 11","exon 12","exon 13"}
ALT_ORDER    = ["exon 2","exon 3","exon 4a","exon 6","exon 10"]

# ── Colour tokens ─────────────────────────────────────────────────────────────
SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, AXIS                 = "#e1e0d9", "#c3c2b7"
_PAL = ["#2a78d6","#eb6834","#1baf7a","#9b59b6","#e67e22","#1abc9c"]

def sample_color(s, labels):
    d = {"Control":"#2a78d6","Disease":"#eb6834","Remedy":"#1baf7a"}
    return d.get(s, _PAL[labels.index(s) % len(_PAL)])


# ── Helpers ───────────────────────────────────────────────────────────────────
def wilson(k, n, z=1.96):
    if n <= 0 or not math.isfinite(k) or not math.isfinite(n):
        return float("nan"), float("nan")
    k = max(0, min(k, n))
    p = k / n; d = 1 + z*z/n
    c = (p + z*z/(2*n)) / d
    h = z * math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / d
    return max(0., c-h), min(1., c+h)

def sem(vals):
    vals = [v for v in vals if math.isfinite(v)]
    if len(vals) < 2: return 0.
    return np.std(vals, ddof=1) / math.sqrt(len(vals))

def style(ax, ylabel="", title="", subtitle=""):
    ax.set_facecolor(SURFACE)
    for sp in ("top","right"):   ax.spines[sp].set_visible(False)
    for sp in ("left","bottom"): ax.spines[sp].set_color(AXIS)
    ax.tick_params(colors=INK2, labelsize=8.5, length=3, color=AXIS)
    ax.yaxis.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)
    if ylabel: ax.set_ylabel(ylabel, fontsize=9, color=INK2)
    if title:
        ax.set_title(title, loc="left", fontsize=10.5, color=INK,
                     fontweight="bold", pad=2)
    if subtitle:
        ax.annotate(subtitle, xy=(0, 1), xycoords="axes fraction",
                    xytext=(0, 26), textcoords="offset points",
                    fontsize=8, color=MUTED, va="bottom")


# ── Panel builders ────────────────────────────────────────────────────────────
def panel_B(ax, agg_df, S, colors):
    """Gene-wide % unspliced — one bar per sample."""
    style(ax, ylabel="% unspliced (gene-wide)",
          title="B  Gene-wide % unspliced",
          subtitle="(EI + IE) / (2 × junction + EI + IE)")
    sub = agg_df.set_index("sample").reindex(S)
    xs  = range(len(S))
    bars = ax.bar(xs, sub["pct_unspliced"].fillna(0), width=0.55,
                  color=colors, edgecolor=SURFACE, lw=1.2, zorder=2)
    ax.errorbar(xs, sub["pct_unspliced"],
                yerr=[(sub["pct_unspliced"] - sub["pct_ci_low"]).clip(lower=0),
                       (sub["pct_ci_high"] - sub["pct_unspliced"]).clip(lower=0)],
                fmt="none", ecolor=INK2, elinewidth=1, capsize=3, zorder=3)
    ax.set_xticks(xs); ax.set_xticklabels(S)
    top = sub["pct_ci_high"].dropna().max() if not sub["pct_ci_high"].dropna().empty else 1
    ax.set_ylim(0, max(top * 1.25, 0.5))
    for bar, s in zip(bars, S):
        v = sub.loc[s, "pct_unspliced"]
        if math.isfinite(v):
            ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+top*0.03,
                    f"{v:.2f}%", ha="center", va="bottom", fontsize=8, color=INK2)


def panel_C(ax, ir_df, S, colors):
    """Mean IR ratio across all introns, weighted by total coverage."""
    style(ax, ylabel="Mean IR ratio (weighted)",
          title="C  Mean intron retention ratio",
          subtitle="mid / (mid + junction), coverage-weighted across all tau introns")
    rows = []
    for s in S:
        sub = ir_df[ir_df["sample"] == s].copy()
        sub["total"] = sub["junction_reads"] + sub["mid_reads"]
        sub = sub[sub["total"] > 0]
        if sub.empty:
            rows.append(dict(sample=s, mean_ir=float("nan"), err=0.))
            continue
        w      = sub["total"]
        wmean  = (sub["IR_ratio"] * w).sum() / w.sum()
        wvar   = (w * (sub["IR_ratio"] - wmean)**2).sum() / w.sum()
        wse    = math.sqrt(wvar / len(sub))
        rows.append(dict(sample=s, mean_ir=wmean, err=wse))
    df = pd.DataFrame(rows).set_index("sample").reindex(S)
    xs   = range(len(S))
    bars = ax.bar(xs, df["mean_ir"].fillna(0), width=0.55,
                  color=colors, edgecolor=SURFACE, lw=1.2, zorder=2)
    ax.errorbar(xs, df["mean_ir"], yerr=df["err"].fillna(0),
                fmt="none", ecolor=INK2, elinewidth=1, capsize=3, zorder=3)
    ax.set_xticks(xs); ax.set_xticklabels(S)
    top = (df["mean_ir"] + df["err"]).dropna().max() if not df.empty else 0.01
    ax.set_ylim(0, max(top * 1.25, 0.001))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.4f}"))
    for bar, (_, row) in zip(bars, df.iterrows()):
        if math.isfinite(row["mean_ir"]):
            ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+top*0.03,
                    f"{row['mean_ir']:.4f}", ha="center", va="bottom",
                    fontsize=8, color=INK2)


def panel_D(ax, psi_df, S, colors):
    """PSI per alternative exon — 5 exon groups, bars per sample."""
    style(ax, ylabel="PSI — % transcripts including exon",
          title="D  Alternative exon PSI",
          subtitle="inc / (inc + skip) × 100   |   exons 2, 3, 4a, 6, 10")
    feats = ALT_ORDER
    bw    = 0.7 / len(S)
    off   = {s: (i-(len(S)-1)/2)*bw for i, s in enumerate(S)}
    top   = 0
    for s, col in zip(S, colors):
        sub = psi_df[psi_df["sample"]==s].set_index("exon").reindex(feats)
        xs  = [i + off[s] for i in range(len(feats))]
        ax.bar(xs, sub["PSI"].fillna(0), width=bw,
               color=col, edgecolor=SURFACE, lw=1.2, zorder=2, label=s)
        ax.errorbar(xs, sub["PSI"],
                    yerr=[(sub["PSI"]-sub["ci_low"]).clip(lower=0),
                           (sub["ci_high"]-sub["PSI"]).clip(lower=0)],
                    fmt="none", ecolor=INK2, elinewidth=0.9, capsize=2, zorder=3)
        t = sub["ci_high"].dropna().max() if not sub["ci_high"].dropna().empty else 0
        top = max(top, t)
    ax.set_xticks(range(len(feats)))
    ax.set_xticklabels([f.replace("exon ","Exon ") for f in feats])
    ax.set_ylim(0, max(top*1.2, 1))
    ax.legend(frameon=False, fontsize=8, loc="upper right")


def panel_E(ax1, ax2, rr_df, S, colors):
    """3R / 4R PSI and RPM."""
    for ax_, col_, title_, ylabel_ in [
        (ax1, "PSI",
         "E  3R vs 4R tau  (%)",
         "% of MAPT transcripts"),
        (ax2, "RPM",
         "E  3R / 4R  (RPM)",
         "Reads per million"),
    ]:
        style(ax_, ylabel=ylabel_, title=title_,
              subtitle="exon 10 PSI — 4R = included, 3R = skipped" if col_ == "PSI" else "")
        bw  = 0.7 / len(S)
        off = {s: (i-(len(S)-1)/2)*bw for i, s in enumerate(S)}
        top = 0
        for s, col in zip(S, colors):
            sub = rr_df[rr_df["sample"]==s].set_index("isoform").reindex(["4R","3R"])
            xs  = [i + off[s] for i in range(2)]
            ax_.bar(xs, sub[col_].fillna(0), width=bw,
                    color=col, edgecolor=SURFACE, lw=1.2, zorder=2)
            if col_ == "PSI":
                ax_.errorbar(xs, sub["PSI"],
                             yerr=[(sub["PSI"]-sub["ci_low"]).clip(lower=0),
                                    (sub["ci_high"]-sub["PSI"]).clip(lower=0)],
                             fmt="none", ecolor=INK2, elinewidth=0.9,
                             capsize=2, zorder=3)
            t = sub[col_].dropna().max() if not sub[col_].dropna().empty else 0
            top = max(top, t)
        ax_.set_xticks([0,1]); ax_.set_xticklabels(["4R","3R"])
        ax_.set_ylim(0, top*1.25)


def panel_F(ax, rpkm_df, S, colors):
    """Mean RPKM of constitutive exons — one bar per sample."""
    style(ax, ylabel="Mean RPKM (constitutive exons)",
          title="F  Overall MAPT expression",
          subtitle="Mean RPKM ± SEM across constitutive exons 1, 4, 5, 7, 9, 11, 12, 13")
    rows = []
    for s in S:
        sub  = rpkm_df[(rpkm_df["sample"]==s) & (rpkm_df["exon"].isin(CONSTITUTIVE))]
        vals = sub["RPKM"].dropna().tolist()
        mean = float(np.mean(vals)) if vals else float("nan")
        se   = sem(vals)
        rows.append(dict(sample=s, mean_rpkm=mean, err=se))
    df   = pd.DataFrame(rows).set_index("sample").reindex(S)
    xs   = range(len(S))
    bars = ax.bar(xs, df["mean_rpkm"].fillna(0), width=0.55,
                  color=colors, edgecolor=SURFACE, lw=1.2, zorder=2)
    ax.errorbar(xs, df["mean_rpkm"], yerr=df["err"].fillna(0),
                fmt="none", ecolor=INK2, elinewidth=1, capsize=3, zorder=3)
    ax.set_xticks(xs); ax.set_xticklabels(S)
    top = (df["mean_rpkm"] + df["err"]).dropna().max() if not df.empty else 1
    ax.set_ylim(0, max(top * 1.25, 0.1))
    for bar, (_, row) in zip(bars, df.iterrows()):
        if math.isfinite(row["mean_rpkm"]):
            ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+top*0.03,
                    f"{row['mean_rpkm']:.3f}", ha="center", va="bottom",
                    fontsize=8, color=INK2)


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--prefix", nargs="+", required=True,
                    help="Prefix(es) used with --out in mapt_complete_analysis.py "
                         "e.g. MAPT_set1  or  MAPT_set1 MAPT_set2")
    ap.add_argument("--labels", nargs="*", default=None,
                    help="Display labels for each prefix (default: use prefix)")
    ap.add_argument("--out", default=None,
                    help="Output prefix for PNG/PDF (default: first prefix + _summary)")
    a = ap.parse_args()

    prefixes = a.prefix
    labels   = a.labels if a.labels else prefixes
    out      = a.out if a.out else prefixes[0] + "_summary"

    if len(labels) != len(prefixes):
        sys.exit("--labels count must match --prefix count")

    # ── Load TSVs ─────────────────────────────────────────────────────────────
    def load(suffix):
        frames = []
        for pfx, lbl in zip(prefixes, labels):
            p = Path(f"{pfx}{suffix}")
            if not p.exists():
                sys.exit(f"File not found: {p}\n"
                         f"Run mapt_complete_analysis.py --out {pfx} first.")
            df = pd.read_csv(p, sep="\t")
            if len(prefixes) > 1:
                df["sample"] = lbl + " · " + df["sample"].astype(str)
            frames.append(df)
        return pd.concat(frames, ignore_index=True)

    agg_df  = load("_aggregate.tsv")
    ir_df   = load("_ir_ratio.tsv")
    psi_df  = load("_psi.tsv")
    rr_df   = load("_3r4r.tsv")
    rpkm_df = load("_exon_rpkm.tsv")

    S      = list(dict.fromkeys(agg_df["sample"].tolist()))
    colors = [sample_color(s, S) for s in S]

    # ── Figure ────────────────────────────────────────────────────────────────
    plt.rcParams.update({"font.family": "DejaVu Sans"})
    fig = plt.figure(figsize=(14, 18), facecolor=SURFACE)
    gs  = fig.add_gridspec(4, 3,
                           height_ratios=[1, 1.2, 1.0, 1.0],
                           width_ratios=[1.1, 1.1, 0.9],
                           hspace=0.80, wspace=0.40)

    ax_B  = fig.add_subplot(gs[0, 0])
    ax_C  = fig.add_subplot(gs[0, 1:])
    ax_D  = fig.add_subplot(gs[1, :])
    ax_E1 = fig.add_subplot(gs[2, 0])
    ax_E2 = fig.add_subplot(gs[2, 1])
    ax_F  = fig.add_subplot(gs[3, :2])
    fig.add_subplot(gs[3, 2]).set_visible(False)

    panel_B(ax_B,  agg_df,  S, colors)
    panel_C(ax_C,  ir_df,   S, colors)
    panel_D(ax_D,  psi_df,  S, colors)
    panel_E(ax_E1, ax_E2, rr_df, S, colors)
    panel_F(ax_F,  rpkm_df, S, colors)

    title_str = " · ".join(labels) if len(labels) > 1 else labels[0]
    fig.suptitle(f"MAPT — Gene-level summary  ({title_str})",
                 fontsize=13, fontweight="bold", color=INK, y=0.99)
    fig.text(0.01, 0.002,
             "B: gene-wide aggregate from _aggregate.tsv.  "
             "C: coverage-weighted mean IR across all tau introns.  "
             "D: per-exon PSI; Wilson 95% CI.  "
             "E: exon 10 PSI.  "
             "F: mean RPKM of constitutive exons ± SEM across exons.",
             fontsize=7, color=MUTED, va="bottom")

    fig.savefig(f"{out}.png", dpi=300, facecolor=SURFACE, bbox_inches="tight")
    fig.savefig(f"{out}.pdf", facecolor=SURFACE, bbox_inches="tight")
    print(f"Saved: {out}.png  {out}.pdf")

    # ── Console summary ───────────────────────────────────────────────────────
    print("\n== Gene-wide % unspliced ==")
    print(agg_df[["sample","pct_unspliced","pct_ci_low","pct_ci_high"]].to_string(index=False))

    print("\n== Alternative exon PSI ==")
    wide = psi_df.pivot(index="exon", columns="sample", values="PSI").reindex(ALT_ORDER)
    print(wide.map(lambda v: f"{v:.2f}" if v==v else "—").to_string())

    print("\n== 3R / 4R ==")
    wide2 = rr_df.pivot(index="isoform", columns="sample", values="PSI")
    print(wide2.map(lambda v: f"{v:.2f}" if v==v else "—").to_string())


if __name__ == "__main__":
    main()
