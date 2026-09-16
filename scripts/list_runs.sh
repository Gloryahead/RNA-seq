#!/usr/bin/env bash
# list_runs.sh
# Display all tracked RNA-seq runs and their status.
#
# Usage:
#   bash scripts/list_runs.sh              # show all runs
#   bash scripts/list_runs.sh --status done    # filter by status
#   bash scripts/list_runs.sh --organism mouse # filter by organism

set -euo pipefail

SCRIPT_BASE="${SCRIPT_BASE:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
MANIFEST="${SCRIPT_BASE}/runs/manifest.tsv"

FILTER_STATUS=""
FILTER_ORGANISM=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --status)   FILTER_STATUS="$2";   shift 2 ;;
    --organism) FILTER_ORGANISM="$2"; shift 2 ;;
    *) echo "Unknown argument: $1"; exit 1 ;;
  esac
done

if [[ ! -f "${MANIFEST}" ]]; then
  echo "No manifest found at ${MANIFEST}"
  exit 1
fi

# Count runs
TOTAL=$(tail -n +2 "${MANIFEST}" | grep -v '^$' | wc -l)

echo ""
echo "═══════════════════════════════════════════════════════════════"
echo "  RNA-seq Run Tracker  |  Total runs: ${TOTAL}"
echo "═══════════════════════════════════════════════════════════════"
printf "  %-35s %-8s %-10s %-10s %s\n" "RUN ID" "STATUS" "ORGANISM" "DATE" "NOTES"
echo "  ---------------------------------------------------------------"

while IFS=$'\t' read -r run_id date organism samples comparisons status notes; do
  [[ "${run_id}" == "run_id" ]] && continue  # skip header
  [[ -z "${run_id}" ]] && continue            # skip empty lines

  # Apply filters
  [[ -n "${FILTER_STATUS}"   && "${status}"   != "${FILTER_STATUS}"   ]] && continue
  [[ -n "${FILTER_ORGANISM}" && "${organism}" != "${FILTER_ORGANISM}" ]] && continue

  # Color-code status
  case "${status}" in
    done)    STATUS_FMT="\033[0;32m${status}\033[0m"    ;;  # green
    running) STATUS_FMT="\033[0;33m${status}\033[0m"   ;;  # yellow
    failed)  STATUS_FMT="\033[0;31m${status}\033[0m"   ;;  # red
    *)       STATUS_FMT="${status}"                      ;;
  esac

  printf "  %-35s %-8b %-10s %-10s %s\n" \
    "${run_id}" "${STATUS_FMT}" "${organism}" "${date}" "${notes}"

done < "${MANIFEST}"

echo ""
echo "  To see details for a run:"
echo "    cat runs/<run_id>/config.yaml"
echo "    cat runs/<run_id>/samples.tsv"
echo ""
echo "  Results live at:"
echo "    /xdisk/haining/maarowosegbe/RNA-seq/results/<run_id>/"
echo "═══════════════════════════════════════════════════════════════"
echo ""
