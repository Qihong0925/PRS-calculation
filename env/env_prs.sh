#!/bin/bash
# Environment for the PRS pipeline. Site-specific settings come from config/local.sh
# (copy config/local.example.sh). Must be sourced, not executed.
_PRS_PKG=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
if [[ -f "$_PRS_PKG/config/local.sh" ]]; then
    source "$_PRS_PKG/config/local.sh"
else
    echo "WARNING: $_PRS_PKG/config/local.sh not found; using conda/tools already on PATH" >&2
fi

# Optional environment modules
if [[ -n "${PRS_MODULES:-}" ]]; then
    command -v module >/dev/null 2>&1 || { echo "ERROR: PRS_MODULES set but 'module' is unavailable" >&2; return 1 2>/dev/null || exit 1; }
    module --force purge 2>/dev/null || module purge
    for m in $PRS_MODULES; do module load "$m"; done
fi

# Clear python paths AFTER module loads (spack modules may re-export them)
unset PYTHONPATH
export PYTHONNOUSERSITE=1

if [[ -n "${CONDA_SH:-}" ]]; then
    source "$CONDA_SH"
elif command -v conda >/dev/null 2>&1; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
else
    echo "ERROR: conda not found; set CONDA_SH in config/local.sh" >&2; return 1 2>/dev/null || exit 1
fi
conda activate "${PRS_CONDA_ENV:-PRS}"

# Force the PRS env first on PATH (modules may pull in another python)
export PATH="$CONDA_PREFIX/bin:$PATH"
hash -r

# Fail loudly if the tools did not load (e.g. job submitted from inside a container)
for t in bcftools samtools bgzip tabix python; do
    command -v "$t" >/dev/null || { echo "ERROR: $t not found after loading env (submitted from a container? MODULEPATH=${MODULEPATH:-})" >&2; return 1 2>/dev/null || exit 1; }
done
