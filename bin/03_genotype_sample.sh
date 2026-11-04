#!/bin/bash
# Step 03: force-genotype one sample at PRS sites.
# Outputs in <outdir>:
#   <sample>.vcf.gz       calls constrained to panel alleles, ALT always kept (used for scoring)
#   <sample>.free.vcf.gz  unconstrained calls in the same regions (QC only, e.g. multi-allelic STRs)
#   <sample>.gt.tsv       per-site GT/DP/AD/GQ; target sites with no overlapping reads at all are
#                         added as GT=./. DP=0 (and listed in <sample>.nocov.tsv)
# A target site without a record is accepted only if samtools confirms zero overlapping reads;
# otherwise the read is treated as truncated and the sample fails.
# Usage: 03_genotype_sample.sh <sample> <cram|bam> <ref_fasta> <sites.tsv.gz> <outdir> <min_mapq> <min_baseq>
set -euo pipefail
SAMPLE=$1; CRAM=$2; REF=$3; SITES=$4; OUTDIR=$5; MIN_MAPQ=$6; MIN_BASEQ=$7
REGIONS="${SITES%.tsv.gz}.regions.tsv.gz"
[[ -f "$REGIONS" ]] || { echo "Missing regions file: $REGIONS (run 03a)" >&2; exit 1; }

cleanup_outputs() { rm -f "$OUTDIR/$SAMPLE".vcf.gz* "$OUTDIR/$SAMPLE".free.vcf.gz* "$OUTDIR/$SAMPLE".gt.tsv "$OUTDIR/$SAMPLE".nocov.tsv; }

mkdir -p "$OUTDIR"
TMP=$(mktemp -d -p "$OUTDIR" ".tmp_${SAMPLE}_XXXX"); trap 'rm -rf "$TMP"' EXIT
echo "$SAMPLE" > "$TMP/name.txt"

# Indexed random access (-R); fail on any htslib read error
bcftools mpileup -f "$REF" -R "$REGIONS" -q "$MIN_MAPQ" -Q "$MIN_BASEQ" -d 1000 \
    -a FORMAT/AD,FORMAT/DP -Ob -o "$TMP/pileup.bcf" "$CRAM" 2> "$TMP/mpileup.err"
cat "$TMP/mpileup.err" >&2
if grep -q -E "^\[E::|rror" "$TMP/mpileup.err"; then echo "ERROR: $SAMPLE read error during mpileup" >&2; exit 1; fi

# Constrained to panel REF,ALT; -A keeps the panel ALT even for 0/0 calls;
# norm trims the repeat context bcftools adds to indel alleles so they match the panel
bcftools call -m -C alleles -A -T "$SITES" -a GQ -Ou "$TMP/pileup.bcf" \
    | bcftools norm -f "$REF" -Oz -o "$TMP/constrained.vcf.gz"
bcftools reheader -s "$TMP/name.txt" -o "$OUTDIR/$SAMPLE.vcf.gz" "$TMP/constrained.vcf.gz"
bcftools index -f -t "$OUTDIR/$SAMPLE.vcf.gz"

# Sanity: every record must carry exactly the panel REF,ALT
n_bad=$(join -t$'\t' \
    <(bcftools query -f '%CHROM:%POS\t%REF,%ALT\n' "$OUTDIR/$SAMPLE.vcf.gz" | sort) \
    <(zcat "$SITES" | awk -F'\t' '{print $1":"$2"\t"$3}' | sort) \
    | awk -F'\t' '$2!=$3' | wc -l)
n_dot=$(bcftools query -f '%ALT\n' "$OUTDIR/$SAMPLE.vcf.gz" | grep -cx '\.' || true)
if [[ $n_bad -ne 0 || $n_dot -ne 0 ]]; then
    echo "ERROR: $SAMPLE has $n_bad records with non-panel alleles and $n_dot with ALT='.'" >&2
    cleanup_outputs; exit 1
fi

# Unconstrained (QC only)
bcftools call -m -a GQ -Ou "$TMP/pileup.bcf" \
    | bcftools norm -f "$REF" -Oz -o "$TMP/free.vcf.gz"
bcftools reheader -s "$TMP/name.txt" -o "$OUTDIR/$SAMPLE.free.vcf.gz" "$TMP/free.vcf.gz"
bcftools index -f -t "$OUTDIR/$SAMPLE.free.vcf.gz"

# Per-site table
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%DP\t%AD\t%GQ]\n' "$OUTDIR/$SAMPLE.vcf.gz" \
    | sed '1s/^# *//; 1s/\[[0-9]*\]//g; 1s/'"$SAMPLE"'://g' > "$OUTDIR/$SAMPLE.gt.tsv"

# Target sites with no record: accept only if there are truly no overlapping reads
n_target=$(zcat "$SITES" | wc -l)
zcat "$SITES" | awk -F'\t' 'NR==FNR{if (FNR>1) seen[$1":"$2]=1; next} !($1":"$2 in seen)' \
    "$OUTDIR/$SAMPLE.gt.tsv" - > "$TMP/missing_sites.tsv"
n_missing=$(wc -l < "$TMP/missing_sites.tsv")
printf "CHROM\tPOS\tREF\tALT\tusable_reads\n" > "$OUTDIR/$SAMPLE.nocov.tsv"
if [[ $n_missing -gt 0 ]]; then
    while IFS=$'\t' read -r chrom pos alleles; do
        ref=${alleles%%,*}; alt=${alleles#*,}
        end=$(( pos + ${#ref} - 1 ))
        # Same read filters as bcftools mpileup defaults + -q: skip unmapped/secondary/QC-fail/dup, low MAPQ, anomalous pairs
        n_reads=$(samtools view -c -T "$REF" -q "$MIN_MAPQ" -F 0x704 -e '!flag.paired || flag.proper_pair' "$CRAM" "$chrom:$pos-$end")
        if [[ $n_reads -gt 0 ]]; then
            echo "ERROR: $SAMPLE has no record at $chrom:$pos but $n_reads overlapping reads (truncated read?)" >&2
            cleanup_outputs; exit 1
        fi
        printf "%s\t%s\t%s\t%s\t0\n" "$chrom" "$pos" "$ref" "$alt" >> "$OUTDIR/$SAMPLE.nocov.tsv"
        printf "%s\t%s\t%s\t%s\t./.\t0\t.\t.\n" "$chrom" "$pos" "$ref" "$alt" >> "$OUTDIR/$SAMPLE.gt.tsv"
    done < "$TMP/missing_sites.tsv"
fi

n_rec=$(($(wc -l < "$OUTDIR/$SAMPLE.gt.tsv") - 1))
if [[ $n_rec -ne $n_target ]]; then
    echo "ERROR: $SAMPLE has $n_rec/$n_target sites in gt.tsv" >&2
    cleanup_outputs; exit 1
fi
n_called=$(awk -F'\t' 'NR>1 && $5!="./." && $5!="."' "$OUTDIR/$SAMPLE.gt.tsv" | wc -l)
mean_dp=$(awk -F'\t' 'NR>1 && $6!="."{s+=$6; n++} END{if(n) printf "%.1f", s/n; else print "NA"}' "$OUTDIR/$SAMPLE.gt.tsv")
echo -e "SUMMARY\t$SAMPLE\ttargets=$n_target\trecords=$((n_rec - n_missing))\tno_coverage=$n_missing\tcalled=$n_called\tmean_DP=$mean_dp"