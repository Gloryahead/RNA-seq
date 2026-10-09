#!/usr/bin/env python3
"""
Convert Ensembl/GENCODE GTF to BED12 format for RSeQC infer_experiment.py.
Usage: python gtf_to_bed12.py annotation.gtf > annotation.bed

One BED12 row per transcript, with exons as blocks. Handles GTF files where
transcript features may appear before or after their exon features.
"""
import sys
import re
from collections import defaultdict

def parse_attrs(s):
    return {m.group(1): m.group(2) for m in re.finditer(r'(\w+)\s+"([^"]+)"', s)}

# tid -> {chrom, strand, tx_start, tx_end, exons}
tx = defaultdict(lambda: {
    "chrom": None, "strand": None, "tx_start": None, "tx_end": None, "exons": []
})

with open(sys.argv[1]) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        chrom, _, feat, start, end, _, strand, _, attrs = f
        if feat not in ("exon", "transcript"):
            continue
        a = parse_attrs(attrs)
        tid = a.get("transcript_id", "")
        if not tid:
            continue
        s0 = int(start) - 1  # GTF 1-based → BED 0-based
        e0 = int(end)
        t = tx[tid]
        t["chrom"]  = chrom
        t["strand"] = strand
        if feat == "transcript":
            t["tx_start"] = s0
            t["tx_end"]   = e0
        else:
            t["exons"].append((s0, e0))

for tid, t in tx.items():
    exons = sorted(t["exons"])
    if not exons or not t["chrom"]:
        continue
    t_start = t["tx_start"] if t["tx_start"] is not None else exons[0][0]
    t_end   = t["tx_end"]   if t["tx_end"]   is not None else exons[-1][1]
    n       = len(exons)
    sizes   = ",".join(str(e - s) for s, e in exons) + ","
    starts  = ",".join(str(s - t_start) for s, e in exons) + ","
    print("\t".join([
        t["chrom"], str(t_start), str(t_end), tid, "0", t["strand"],
        str(t_start), str(t_end), "0", str(n), sizes, starts
    ]))
