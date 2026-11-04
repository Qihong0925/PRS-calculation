#!/bin/bash
# SLURM array wrapper for step 03: task N genotypes data row N of the samples table.
# Samples table: tab-separated with header; columns 'sample' and 'cram' (resolved by name).
# Usage (via sbatch --array): 03_genotype_array.sh <pkg> <samples.tsv> <ref> <sites> <outdir> <min_mapq> <min_baseq>
PKG=$1; SAMPLES=$2; REF=$3; SITES=$4; OUTDIR=$5; MIN_MAPQ=$6; MIN_BASEQ=$7
source "$PKG/env/env_prs.sh" || exit 1
set -euo pipefail

read -r SAMPLE CRAM < <(awk -F'\t' -v n="$SLURM_ARRAY_TASK_ID" '
    NR==1 {for (i=1; i<=NF; i++) c[$i]=i; if (!("sample" in c) || !("cram" in c)) {print "ERR"; exit}; next}
    NR==n+1 {print $c["sample"], $c["cram"]}' "$SAMPLES")
[[ -n "${SAMPLE:-}" && "$SAMPLE" != "ERR" ]] || { echo "Bad samples table or task $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
[[ -f "$CRAM" ]] || { echo "Missing CRAM: $CRAM" >&2; exit 1; }

echo "Task $SLURM_ARRAY_TASK_ID: $SAMPLE $CRAM"
bash "$PKG/bin/03_genotype_sample.sh" "$SAMPLE" "$CRAM" "$REF" "$SITES" "$OUTDIR" "$MIN_MAPQ" "$MIN_BASEQ"
