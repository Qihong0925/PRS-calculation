#!/bin/bash
# Step 04: merge per-sample constrained calls into a cohort VCF and score each panel.
# Optional: set FORCE_EXPECTED=<tsv> (column 'site' CHROM:POS or 'ID') to treat those sites as
# untrusted in this batch; they are scored at expected dosage 2 x EAF for every sample.
# Usage: [FORCE_EXPECTED=file] 04_merge_score.sh <pkg> <samples.tsv> <gt_dir> <outdir> <panel_dir> [<panel_dir> ...]
set -euo pipefail
PKG=$1; SAMPLES=$2; GT_DIR=$3; OUTDIR=$4; shift 4

mkdir -p "$OUTDIR"
LIST="$OUTDIR/vcf_list.txt"
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; next} {print $c["sample"]}' "$SAMPLES" \
    | while read -r s; do echo "$GT_DIR/$s.vcf.gz"; done > "$LIST"

missing=0
while read -r f; do [[ -f "$f" && -f "$f.tbi" ]] || { echo "MISSING: $f" >&2; missing=$((missing+1)); }; done < "$LIST"
[[ $missing -eq 0 ]] || { echo "ERROR: $missing sample VCFs missing; finish step 03 first" >&2; exit 1; }

COHORT="$OUTDIR/cohort.prs_sites.vcf.gz"
bcftools merge -l "$LIST" -m none -Oz -o "$COHORT"
bcftools index -f -t "$COHORT"
echo "Cohort VCF: $(bcftools query -l "$COHORT" | wc -l) samples, $(bcftools view -H "$COHORT" | wc -l) sites with >=1 record"

FE_ARG=()
if [[ -n "${FORCE_EXPECTED:-}" ]]; then
    [[ -f "$FORCE_EXPECTED" ]] || { echo "ERROR: FORCE_EXPECTED file not found: $FORCE_EXPECTED" >&2; exit 1; }
    FE_ARG=(--force-expected "$FORCE_EXPECTED")
    echo "Force-expected sites: $FORCE_EXPECTED"
fi

for PANEL in "$@"; do
    name=$(basename "$PANEL")
    python "$PKG/py/score_prs.py" \
        --vcf "$COHORT" \
        --panel-dir "$PANEL" \
        --standardisation "$PKG/resources/canrisk_standardisation.tsv" \
        "${FE_ARG[@]}" \
        --out-prefix "$OUTDIR/$name"
done
