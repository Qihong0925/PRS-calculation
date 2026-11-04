#!/bin/bash
# Step 01: build GRCh38 PRS sites and plink2 score file from a CanRisk PRS file (GRCh37).
# Usage: 01_make_sites.sh <canrisk_prs> <chain> <ref_fasta> <corrections_tsv> <pgs_check|none> <outdir>
set -euo pipefail

CANRISK=$1
CHAIN=$2
REF_FASTA=$3
CORRECTIONS=$4
PGS_CHECK=$5
OUTDIR=$6
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

for f in "$CANRISK" "$CHAIN" "${REF_FASTA}.fai" "$CORRECTIONS"; do
    [[ -f "$f" ]] || { echo "Missing input: $f" >&2; exit 1; }
done
mkdir -p "$OUTDIR"

python "$SCRIPT_DIR/../py/make_sites.py" \
    --canrisk "$CANRISK" \
    --chain "$CHAIN" \
    --ref-fasta "$REF_FASTA" \
    --corrections "$CORRECTIONS" \
    --pgs-check "$PGS_CHECK" \
    --outdir "$OUTDIR"

bgzip -f "$OUTDIR/sites.tsv"
tabix -f -s1 -b2 -e2 "$OUTDIR/sites.tsv.gz"

echo "Step 01 done: $OUTDIR"
