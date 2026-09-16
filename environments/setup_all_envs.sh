#!/usr/bin/env bash
# setup_all_envs.sh
# Creates all conda environments for the NGS101 RNA-seq tutorial pipeline.
# Run from the environments/ directory.
# Usage: bash setup_all_envs.sh [--r-extras]
#   --r-extras   Also run R install scripts after each relevant env is built.
#
# Prerequisites:
#   - Miniforge (https://github.com/conda-forge/miniforge) or Anaconda/Miniconda
#   - mamba installed: conda install -n base -c conda-forge mamba
# Estimated time: 30-90 minutes depending on network speed.
# Estimated disk: ~20-40 GB across all environments.

set -euo pipefail

# UA HPC: point micromamba at the haining group prefix
export MAMBA_ROOT_PREFIX=/groups/haining/maarowosegbe/micromamba
export MAMBA_EXE=/opt/ohpc/pub/apps/micromamba/2.0.2-2/bin/micromamba

R_EXTRAS=false
[[ "${1:-}" == "--r-extras" ]] && R_EXTRAS=true

ENVDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_DIR="${ENVDIR}/r_scripts"

build_env() {
  local yml="$1"
  local name
  name=$(grep '^name:' "${ENVDIR}/${yml}" | awk '{print $2}')
  echo ""
  echo "══════════════════════════════════════════════════════"
  echo "  Building: ${name}  (${yml})"
  echo "══════════════════════════════════════════════════════"
  if micromamba env list | grep -q "^${name} "; then
    echo "  → Already exists; skipping. To rebuild: micromamba env remove -n ${name}"
  else
    micromamba env create -f "${ENVDIR}/${yml}"
    echo "  ✓ ${name} created"
  fi
}

run_r_script() {
  local env_name="$1"
  local script="$2"
  echo "  → Running R install script: ${script}"
  micromamba run -n "${env_name}" Rscript "${SCRIPT_DIR}/${script}"
}

# ── Core environments ─────────────────────────────────────────────────────────
build_env 01_rnaseq_core.yml
build_env 02_rnaseq_r.yml
$R_EXTRAS && run_r_script rnaseq_r_env install_core_r_extras.R

# ── Specialized analysis environments ────────────────────────────────────────
build_env 03_cancer_subtype_r.yml
$R_EXTRAS && run_r_script cancer_subtype_env install_cancer_subtype_r.R

build_env 04_spladder.yml
build_env 05_reditools2.yml
build_env 06_suppa2.yml
build_env 07_circrna.yml
build_env 08_ciri3.yml
$R_EXTRAS && run_r_script ciri3_env install_ciri3_r.R

build_env 09_mirna.yml
build_env 10_starfusion.yml
build_env 11_fusioncatcher.yml
build_env 12_esviritu.yml
build_env 13_virtus2.yml

build_env 14_mageck.yml
$R_EXTRAS && run_r_script mageck_env install_mageck_r.R

build_env 15_wgcna_r.yml
build_env 16_genie3_masterreg_r.yml
$R_EXTRAS && run_r_script genie3_masterreg_env install_genie3_masterreg_r.R

echo ""
echo "══════════════════════════════════════════════════════"
echo "  All environments built successfully."
echo "══════════════════════════════════════════════════════"
echo "  NOTE — manual post-install steps required:"
echo "    05_reditools2: git clone + pip install (Python 2.7, see YAML)"
echo "    06_suppa2:     git clone comprna/SUPPA + pip install -e ."
echo "    08_ciri3:      download CIRI3 jar from GitHub releases"
echo "    10_starfusion: download CTAT genome lib (~30 GB) from Broad"
echo "    11_fusioncatcher: run fusioncatcher-build for each organism"
echo "    12_esviritu:   esviritu download_db -o ~/esviritu_db"
echo "    13_virtus2:    git clone yyoshiaki/VIRTUS2"
echo "    Part 7 (CIBERSORTx): web-only — register at cibersortx.stanford.edu"
