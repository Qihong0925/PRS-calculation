#!/bin/bash
# Step 02: validate resolved PRS sites against gnomAD (flags only, nothing is removed).
# Usage: 02_validate_sites.sh <sites_dir> <gnomad_dir> <gnomad_pattern> <outdir>
set -euo pipefail

SITES_DIR=$1
GNOMAD_DIR=$2
GNOMAD_PATTERN=$3
OUTDIR=$4
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

mkdir -p "$OUTDIR"
[[ -f "$SITES_DIR/score_plink2.tsv" ]] || { echo "Run step 01 first" >&2; exit 1; }

python "$SCRIPT_DIR/../py/validate_sites.py" \
    --score-plink2 "$SITES_DIR/score_plink2.tsv" \
    --gnomad-dir "$GNOMAD_DIR" \
    --gnomad-pattern "$GNOMAD_PATTERN" \
    --af-diff-threshold 0.1 \
    --threads "${SLURM_CPUS_PER_TASK:-4}" \
    --outdir "$OUTDIR"

echo "Step 02 done: $OUTDIR"
