#!/usr/bin/env bash
# setup_all_envs.sh
# Creates all micromamba environments for the RNA-seq pipeline.
# Each YAML file defines its own environment name.
# Continues past individual env failures — run again to retry any that failed.
#
# Usage:
#   bash environments/setup_all_envs.sh
#
# Run on an interactive node (not login node):
#   interactive -a haining -n 8 --mem=32gb -t 04:00:00

set -uo pipefail

export MAMBA_ROOT_PREFIX=/groups/haining/maarowosegbe/micromamba
export MAMBA_EXE=/opt/ohpc/pub/apps/micromamba/2.0.2-2/bin/micromamba
MC="${MAMBA_EXE}"

ENVDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

failed=()

for yml in "${ENVDIR}"/*.yml; do
    # tr -d '\r' strips Windows carriage returns if files were edited on Windows
    name=$(grep '^name:' "${yml}" | awk '{print $2}' | tr -d '\r')
    echo ""
    echo "=== ${name} ($(basename "${yml}")) ==="
    if "${MC}" env list 2>/dev/null | awk '{print $1}' | grep -qx "${name}"; then
        echo "  already exists, skipping"
    elif "${MC}" env create -f "${yml}" --yes; then
        echo "  done"
    else
        echo "  FAILED — skipping"
        echo "  Retry:  ${MC} env create -f ${yml} --yes"
        failed+=("${name}")
    fi
done

echo ""
if [[ ${#failed[@]} -gt 0 ]]; then
    echo "=== Build complete — ${#failed[@]} environment(s) FAILED ==="
    printf '  - %s\n' "${failed[@]}"
    exit 1
else
    echo "=== All environments built successfully ==="
fi
