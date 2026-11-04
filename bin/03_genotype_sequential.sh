#!/bin/bash
# Step 03 (sequential): genotype samples one at a time, for CRAMs on tiered/archived storage
# where concurrent reads time out. For each sample: read the index first (forces recall),
# then genotype; retry after a wait on failure.
# Samples table: tab-separated with header; columns 'sample' and 'cram' (resolved by name).
# Usage: 03_genotype_sequential.sh <pkg> <samples.tsv> <ref> <sites> <outdir> <min_mapq> <min_baseq> [max_attempts] [wait_sec]
PKG=$1; SAMPLES=$2; REF=$3; SITES=$4; OUTDIR=$5; MIN_MAPQ=$6; MIN_BASEQ=$7
MAX_ATTEMPTS=${8:-4}; WAIT=${9:-600}
source "$PKG/env/env_prs.sh" || exit 1
set -uo pipefail

n_ok=0; n_fail=0
while read -r SAMPLE CRAM; do
    if [[ -f "$OUTDIR/$SAMPLE.gt.tsv" ]]; then echo "SKIP $SAMPLE (already done)"; continue; fi
    IDX=""; for x in "$CRAM.crai" "$CRAM.bai" "${CRAM%.*}.crai" "${CRAM%.*}.bai"; do [[ -f "$x" ]] && { IDX=$x; break; }; done
    [[ -n "$IDX" ]] || { echo "FAILED $SAMPLE (no index)"; n_fail=$((n_fail+1)); continue; }
    ok=0
    for attempt in $(seq 1 "$MAX_ATTEMPTS"); do
        echo "$(date '+%F %T') $SAMPLE attempt $attempt/$MAX_ATTEMPTS"
        if cat "$IDX" > /dev/null && \
           bash "$PKG/bin/03_genotype_sample.sh" "$SAMPLE" "$CRAM" "$REF" "$SITES" "$OUTDIR" "$MIN_MAPQ" "$MIN_BASEQ"; then
            echo "$(date '+%F %T') DONE $SAMPLE (attempt $attempt)"; ok=1; break
        fi
        rm -f "$OUTDIR/$SAMPLE".vcf.gz* "$OUTDIR/$SAMPLE".free.vcf.gz* "$OUTDIR/$SAMPLE".gt.tsv
        [[ $attempt -lt $MAX_ATTEMPTS ]] && { echo "$(date '+%F %T') failed; waiting ${WAIT}s"; sleep "$WAIT"; }
    done
    if [[ $ok -eq 1 ]]; then n_ok=$((n_ok+1)); else n_fail=$((n_fail+1)); echo "FAILED $SAMPLE"; fi
done < <(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; next} {print $c["sample"], $c["cram"]}' "$SAMPLES")

echo "Sequential run finished: $n_ok done, $n_fail failed"
[[ $n_fail -eq 0 ]]
