#!/usr/bin/env bash
# new_run.sh
# Create a new tracked RNA-seq analysis run.
#
# Usage:
#   bash scripts/new_run.sh --name "mouse_pancreas" --notes "KRAS vs KRAS+SPIB"
#   bash scripts/new_run.sh --name "human_tcga" --organism human --notes "TCGA BRCA cohort"
#
# What it does:
#   1. Assigns a unique run ID: YYYY-MM-DD_<name>
#   2. Creates results/<run_id>/ on xdisk
#   3. Snapshots current config/config.yaml into runs/<run_id>/
#   4. Snapshots current config/samples.tsv into runs/<run_id>/
#   5. Registers the run in runs/manifest.tsv
#   6. Prints the sbatch command to launch the pipeline for this run

set -euo pipefail

SCRIPT_BASE="${SCRIPT_BASE:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
DATA_BASE="${DATA_BASE:-/xdisk/haining/maarowosegbe/RNA-seq}"
MANIFEST="${SCRIPT_BASE}/runs/manifest.tsv"

# ── Parse arguments ───────────────────────────────────────────────────
NAME=""
NOTES=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --name)     NAME="$2";     shift 2 ;;
    --notes)    NOTES="$2";   shift 2 ;;
    *) echo "Unknown argument: $1"; exit 1 ;;
  esac
done

if [[ -z "${NAME}" ]]; then
  echo "ERROR: --name is required"
  echo "Usage: bash scripts/new_run.sh --name <short_name> [--notes \"description\"]"
  exit 1
fi

# ── Build run ID ──────────────────────────────────────────────────────
DATE=$(date +%Y-%m-%d)
RUN_ID="${DATE}_${NAME}"

# Guard against duplicate run IDs
if grep -q "^${RUN_ID}" "${MANIFEST}" 2>/dev/null; then
  echo "ERROR: Run ID '${RUN_ID}' already exists in manifest."
  echo "       Use a different --name or check runs/manifest.tsv"
  exit 1
fi

# ── Read metadata from config ─────────────────────────────────────────
CONFIG="${SCRIPT_BASE}/config/config.yaml"
ORGANISM=$(grep '^organism:' "${CONFIG}" | awk '{print $2}' | tr -d '"')
SAMPLES_TSV="${SCRIPT_BASE}/config/samples.tsv"
SAMPLE_LIST=$(tail -n +2 "${SAMPLES_TSV}" | grep -v '^#' | awk -F'\t' '{print $1}' | paste -sd ',' -)
COMPARISONS=$(grep -A5 '^comparisons:' "${CONFIG}" | grep '^\s*-' | tr -d ' -[]"' | paste -sd ',' -)

# ── Create run snapshot directory ────────────────────────────────────
RUN_DIR="${SCRIPT_BASE}/runs/${RUN_ID}"
mkdir -p "${RUN_DIR}"
cp "${CONFIG}"      "${RUN_DIR}/config.yaml"
cp "${SAMPLES_TSV}" "${RUN_DIR}/samples.tsv"

# ── Create results directory on xdisk ────────────────────────────────
RESULTS_DIR="${DATA_BASE}/results/${RUN_ID}"
mkdir -p "${RESULTS_DIR}"

# ── Register in manifest ──────────────────────────────────────────────
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "${RUN_ID}" "${DATE}" "${ORGANISM}" "${SAMPLE_LIST}" \
  "${COMPARISONS}" "created" "${NOTES}" \
  >> "${MANIFEST}"

# ── Write a run-specific launcher script ─────────────────────────────
LAUNCHER="${RUN_DIR}/run.sh"
cat > "${LAUNCHER}" <<LAUNCHER
#!/usr/bin/env bash
# Auto-generated launcher for run: ${RUN_ID}
# Created: ${DATE}
# Organism: ${ORGANISM}
# Notes: ${NOTES}
#
# Submit with: sbatch ${LAUNCHER}

#SBATCH --job-name=${RUN_ID}
#SBATCH --output=${DATA_BASE}/logs/slurm/${RUN_ID}-%j.out
#SBATCH --error=${DATA_BASE}/logs/slurm/${RUN_ID}-%j.err
#SBATCH --account=haining
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=240:00:00

SCRIPT_BASE=/home/u11/maarowosegbe/RNA-seq
DATA_BASE=/xdisk/haining/maarowosegbe/RNA-seq
RUN_ID=${RUN_ID}
MANIFEST=\${SCRIPT_BASE}/runs/manifest.tsv

mkdir -p "\${DATA_BASE}/logs/slurm"

set +eu; source ~/.bashrc; set -eu
eval "\$(mamba-haining shell hook --shell=bash)"
mamba-haining activate snakemake_env

# Mark run as running in manifest
sed -i "s|^\${RUN_ID}\t.*\tcreated\t|\${RUN_ID}\t${DATE}\t${ORGANISM}\t${SAMPLE_LIST}\t${COMPARISONS}\trunning\t|" "\${MANIFEST}"

echo "========== Run \${RUN_ID} started: \$(date) =========="

snakemake \\
  --use-conda \\
  --profile \${SCRIPT_BASE}/workflow/profiles/slurm \\
  --config outdir="\${DATA_BASE}/results/\${RUN_ID}" \\
           samples="\${SCRIPT_BASE}/runs/\${RUN_ID}/samples.tsv" \\
  --configfile "\${SCRIPT_BASE}/runs/\${RUN_ID}/config.yaml" \\
  --cores 2 \\
  --rerun-incomplete \\
  2>"\${DATA_BASE}/logs/slurm/\${RUN_ID}_snakemake.log"

EXIT_CODE=\$?

# Update manifest status
if [[ \${EXIT_CODE} -eq 0 ]]; then
  STATUS="done"
else
  STATUS="failed"
fi
sed -i "s|^\${RUN_ID}\t.*\trunning\t|\${RUN_ID}\t${DATE}\t${ORGANISM}\t${SAMPLE_LIST}\t${COMPARISONS}\t\${STATUS}\t|" "\${MANIFEST}"

echo "========== Run \${RUN_ID} \${STATUS}: \$(date) =========="
exit \${EXIT_CODE}
LAUNCHER
chmod +x "${LAUNCHER}"

# ── Print summary ─────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "  New run created: ${RUN_ID}"
echo "============================================================"
echo "  Organism    : ${ORGANISM}"
echo "  Samples     : ${SAMPLE_LIST}"
echo "  Comparisons : ${COMPARISONS}"
echo "  Notes       : ${NOTES}"
echo ""
echo "  Config snapshot : runs/${RUN_ID}/config.yaml"
echo "  Sample snapshot : runs/${RUN_ID}/samples.tsv"
echo "  Results dir     : ${RESULTS_DIR}/"
echo "  Launcher        : runs/${RUN_ID}/run.sh"
echo ""
echo "  To launch:"
echo "    sbatch runs/${RUN_ID}/run.sh"
echo "============================================================"
