#!/bin/bash
# Step 00 (per batch): confirm every panel REF allele matches the batch reference FASTA.
# Usage: 00_check_reference.sh <ref_fasta> <panel_dir> [<panel_dir> ...]
set -euo pipefail
REF=$1; shift
[[ -f "$REF.fai" ]] || { echo "Missing $REF.fai" >&2; exit 1; }
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
fail=0
for PANEL in "$@"; do
    name=$(basename "$PANEL")
    zcat "$PANEL/sites.tsv.gz" | awk -F'\t' '{split($3,a,","); print $1":"$2"-"$2+length(a[1])-1"\t"a[1]}' > "$TMP/$name.tsv"
    cut -f1 "$TMP/$name.tsv" > "$TMP/$name.regions"
    samtools faidx "$REF" -r "$TMP/$name.regions" -n 100000 | grep -v '^>' | tr a-z A-Z > "$TMP/$name.seq"
    n_bad=$(paste "$TMP/$name.tsv" "$TMP/$name.seq" \
        | awk -F'\t' -v p="$name" '$2!=$3{print "MISMATCH\t"p"\t"$1"\tpanel="$2"\tref="$3 > "/dev/stderr"; c++} END{print c+0}')
    echo "$name: $(wc -l < "$TMP/$name.tsv") sites, $n_bad REF mismatches"
    [[ $n_bad -eq 0 ]] || fail=1
done
[[ $fail -eq 0 ]] || { echo "ERROR: panel REF does not match $REF" >&2; exit 1; }
