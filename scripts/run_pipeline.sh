#!/usr/bin/env bash
# run_pipeline.sh
# Convenience wrapper to launch the full Snakemake pipeline on SLURM.
# Run this after setup is complete (environments built, indices made).
#
# Usage:
#   bash scripts/run_pipeline.sh                    # full pipeline, conda mode
#   bash scripts/run_pipeline.sh --singularity      # use Apptainer containers
#   bash scripts/run_pipeline.sh --dry-run          # check DAG without running
#   bash scripts/run_pipeline.sh --unlock           # unlock stale directory
#   bash scripts/run_pipeline.sh --report           # generate HTML report

set -euo pipefail

MODE="conda"          # conda | singularity
DRY_RUN=false
UNLOCK=false
REPORT=false
JOBS=100
BIND_PATHS="/ref,/data,/scratch"

for arg in "$@"; do
  case "$arg" in
    --singularity)  MODE="singularity" ;;
    --dry-run|-n)   DRY_RUN=true ;;
    --unlock)       UNLOCK=true ;;
    --report)       REPORT=true ;;
  esac
done

# ── Sanity checks ─────────────────────────────────────────────────────
[[ -f "workflow/Snakefile" ]] || {
  echo "ERROR: workflow/Snakefile not found. Are you in the repo root?"
  exit 1
}
[[ -f "config/config.yaml" ]] || {
  echo "ERROR: config/config.yaml not found."
  exit 1
}
[[ -f "config/samples.tsv" ]] || {
  echo "ERROR: config/samples.tsv not found."
  exit 1
}

if $UNLOCK; then
  snakemake --unlock --snakefile workflow/Snakefile
  echo "Directory unlocked."
  exit 0
fi

# ── Build Snakemake command ───────────────────────────────────────────
CMD=(
  snakemake
    --snakefile workflow/Snakefile
    --profile   workflow/profiles/slurm
    --jobs      "${JOBS}"
    --rerun-incomplete
    --keep-going
    --latency-wait 60
    --printshellcmds
)

case "${MODE}" in
  conda)
    CMD+=(--use-conda --conda-frontend mamba)
    ;;
  singularity)
    [[ -f "containers/rnaseq_core.sif" ]] || {
      echo "ERROR: SIF files not built. Run: bash containers/build_all_sifs.sh"
      exit 1
    }
    CMD+=(
      --use-singularity
      --singularity-args "--bind ${BIND_PATHS}"
    )
    ;;
esac

$DRY_RUN  && CMD+=(--dry-run)
$REPORT   && CMD+=(--report results/snakemake_report.html)

echo "=== Launching Snakemake pipeline ==="
echo "  Mode    : ${MODE}"
echo "  Dry run : ${DRY_RUN}"
echo "  Command : ${CMD[*]}"
echo ""

"${CMD[@]}"

if ! $DRY_RUN && ! $REPORT; then
  echo ""
  echo "=== Pipeline finished ==="
  echo "  Run MultiQC to aggregate final QC:"
  echo "  bash scripts/06_multiqc.sh"
fi
