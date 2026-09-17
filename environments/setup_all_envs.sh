#!/usr/bin/env bash
# setup_all_envs.sh
# Creates all conda environments for the RNA-seq pipeline.
# Each YAML file defines its own environment name.
#
# Usage:
#   bash environments/setup_all_envs.sh
#
# Run on an interactive node (not login node):
#   interactive -a haining -n 8 --mem=32gb -t 04:00:00

set -euo pipefail

export MAMBA_ROOT_PREFIX=/groups/haining/maarowosegbe/micromamba
export MAMBA_EXE=/opt/ohpc/pub/apps/micromamba/2.0.2-2/bin/micromamba
MC="${MAMBA_EXE}"

ENVDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

for yml in "${ENVDIR}"/*.yml; do
    name=$(grep '^name:' "${yml}" | awk '{print $2}')
    echo ""
    echo "=== ${name} ($(basename "${yml}")) ==="
    if [[ -d "${MAMBA_ROOT_PREFIX}/envs/${name}" ]]; then
        echo "  already exists, skipping"
    else
        "${MC}" env create -f "${yml}" --yes
        echo "  done"
    fi
done

echo ""
echo "=== All environments built ==="
