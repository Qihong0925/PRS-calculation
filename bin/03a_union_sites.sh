#!/bin/bash
# Step 03a: union of panel sites so each sample is genotyped once for all panels.
# Writes <out>.tsv.gz (CHROM POS REF,ALT; for bcftools call -T) and
#        <out>.regions.tsv.gz (CHROM START END covering each allele; for indexed mpileup -R)
# Usage: 03a_union_sites.sh <out_sites.tsv.gz> <panel_dir> [<panel_dir> ...]
set -euo pipefail
OUT=$1; shift
REGIONS="${OUT%.tsv.gz}.regions.tsv.gz"
TMP=$(mktemp)
for P in "$@"; do zcat "$P/sites.tsv.gz"; done | sort -t$'\t' -k1,1V -k2,2n -k3,3 -u > "$TMP"
dup=$(cut -f1,2 "$TMP" | uniq -d | wc -l)
if [[ $dup -ne 0 ]]; then
    echo "ERROR: $dup positions have different alleles across panels:" >&2
    cut -f1,2 "$TMP" | uniq -d >&2
    rm -f "$TMP"; exit 1
fi
bgzip -c "$TMP" > "$OUT"
tabix -f -s1 -b2 -e2 "$OUT"
awk -F'\t' 'BEGIN{OFS="\t"} {split($3,a,","); L=length(a[1]); if (length(a[2])>L) L=length(a[2]); print $1,$2,$2+L-1}' "$TMP" \
    | bgzip -c > "$REGIONS"
tabix -f -s1 -b2 -e3 "$REGIONS"
rm -f "$TMP"
echo "Union sites: $(zcat "$OUT" | wc -l) -> $OUT (+ regions file)"
