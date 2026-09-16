#!/usr/bin/env bash
# build_all_sifs.sh
# Builds all Apptainer SIF containers from .def files.
# Run this ONCE before launching the Snakemake pipeline with --use-singularity.
# Requires: apptainer >= 1.2 installed (module load apptainer on most HPCs)
#
# Usage:
#   cd <repo_root>
#   bash containers/build_all_sifs.sh
#
# To build a single SIF:
#   apptainer build containers/rnaseq_core.sif containers/rnaseq_core.def
#
# After all SIFs are built, run the pipeline:
#   snakemake --use-singularity \
#             --singularity-args "--bind /ref,/data,/scratch" \
#             --profile workflow/profiles/slurm
#
# SIF sizes (approximate):
#   rnaseq_core.sif   ~2.5 GB
#   rnaseq_r.sif      ~3.5 GB
#   mirna.sif         ~1.5 GB
#   mageck.sif        ~2.0 GB
#   (total ~20-30 GB for all environments)

set -euo pipefail

CONTAINER_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "${CONTAINER_DIR}")"

module load apptainer 2>/dev/null || true  # HPC module; skip if not needed

build_sif() {
  local def="$1"
  local sif="$2"
  local name
  name=$(basename "${sif}" .sif)

  echo ""
  echo "═══════════════════════════════════════════════════"
  echo "  Building: ${name}.sif"
  echo "═══════════════════════════════════════════════════"

  if [[ -f "${sif}" ]]; then
    echo "  → Already exists; skipping. Remove to rebuild."
    return
  fi

  # Build from repo root so %files paths resolve correctly
  (cd "${REPO_ROOT}" && apptainer build "${sif}" "${def}")
  echo "  ✓ ${sif} built ($(du -sh "${sif}" | cut -f1))"
}

# ── Core environments (build first — others may depend on them) ────────
build_sif containers/rnaseq_core.def containers/rnaseq_core.sif
build_sif containers/rnaseq_r.def    containers/rnaseq_r.sif

# ── Specialized environments ──────────────────────────────────────────
# Each .def file follows the same pattern as rnaseq_core.def:
# Bootstrap: docker / From: mambaorg/micromamba:1.5.8
# %files: YAML → /opt/<name>.yml
# %post:  micromamba install -y -n base -f /opt/<name>.yml
# Create the missing .def files from their YAML counterparts as needed.

build_sif containers/cancer_subtype.def   containers/cancer_subtype.sif
build_sif containers/spladder.def         containers/spladder.sif
build_sif containers/reditools2.def       containers/reditools2.sif
build_sif containers/suppa2.def           containers/suppa2.sif
build_sif containers/circrna.def          containers/circrna.sif
build_sif containers/ciri3.def            containers/ciri3.sif
build_sif containers/mirna.def            containers/mirna.sif
build_sif containers/starfusion.def       containers/starfusion.sif
build_sif containers/fusioncatcher.def    containers/fusioncatcher.sif
build_sif containers/esviritu.def         containers/esviritu.sif
build_sif containers/virtus2.def          containers/virtus2.sif
build_sif containers/mageck.def           containers/mageck.sif
build_sif containers/wgcna.def            containers/wgcna.sif
build_sif containers/genie3_masterreg.def containers/genie3_masterreg.sif

echo ""
echo "═══════════════════════════════════════════════════"
echo "  All SIFs built. Now run:"
echo "  snakemake --use-singularity \\"
echo "    --singularity-args \"--bind /ref,/data\" \\"
echo "    --profile workflow/profiles/slurm"
echo "═══════════════════════════════════════════════════"
