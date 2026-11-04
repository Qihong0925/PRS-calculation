#!/usr/bin/env python3
"""Validate PRS sites against gnomAD genomes (flags only; nothing is removed).

match_status: exact / orientation_flip_candidate / site_present_allele_absent /
              absent_in_gnomad / no_gnomad_file
af_flag: |gnomAD NFE effect-allele freq - CanRisk EAF| > threshold
gnomAD records are cached; the cache is reused only if the queried positions are identical.
"""
import argparse
import os
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd

GNOMAD_FIELDS = ["CHROM", "POS", "REF", "ALT", "FILTER", "AF", "AF_nfe"]


def to_float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


def query_gnomad(vcf, chrom, positions, tmpdir):
    reg = os.path.join(tmpdir, f"{chrom}.regions.tsv")
    with open(reg, "w") as fh:
        for p in sorted(set(positions)):
            fh.write(f"{chrom}\t{p}\t{p}\n")
    cmd = ["bcftools", "query", "-R", reg,
           "-f", "%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/AF\t%INFO/AF_nfe\n", vcf]
    out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    rows = [l.split("\t") for l in out.splitlines() if l]
    return pd.DataFrame(rows, columns=GNOMAD_FIELDS)


def extract_gnomad(sc, gnomad_dir, pattern, threads, tmpdir):
    jobs = {}
    with ThreadPoolExecutor(max_workers=threads) as ex:
        for chrom, sub in sc.groupby("CHROM"):
            vcf = os.path.join(gnomad_dir, pattern.format(chrom=chrom))
            if not os.path.exists(vcf):
                continue
            jobs[ex.submit(query_gnomad, vcf, chrom, sub["POS"].tolist(), tmpdir)] = (chrom, len(sub))
        print(f"Querying gnomAD for {len(jobs)} chromosomes with {threads} threads", flush=True)
        t0, parts = time.time(), []
        for i, fut in enumerate(as_completed(jobs), 1):
            chrom, n = jobs[fut]
            parts.append(fut.result())
            print(f"  [{i}/{len(jobs)}] {chrom}: {n} sites done, {time.time() - t0:.0f}s elapsed", flush=True)
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=GNOMAD_FIELDS)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--score-plink2", required=True)
    ap.add_argument("--gnomad-dir", required=True)
    ap.add_argument("--gnomad-pattern", required=True)
    ap.add_argument("--af-diff-threshold", type=float, default=0.1)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    sc = pd.read_csv(args.score_plink2, sep="\t", dtype=str)
    if "EAF" not in sc.columns:
        sys.exit("score_plink2.tsv has no EAF column; re-run step 01")
    sc["CHROM"] = sc["ID"].str.split(":").str[0]
    sc["POS"] = sc["ID"].str.split(":").str[1].astype(int)
    sc["canrisk_eaf"] = sc["EAF"].map(to_float)

    missing_files = sorted(c for c in sc["CHROM"].unique()
                           if not os.path.exists(os.path.join(args.gnomad_dir, args.gnomad_pattern.format(chrom=c))))
    if missing_files:
        print(f"WARNING: no gnomAD file for {missing_files}", file=sys.stderr)

    cache = os.path.join(args.outdir, "gnomad_prs_sites.tsv")
    cache_pos = os.path.join(args.outdir, "gnomad_prs_sites.positions.tsv")
    want = set(zip(sc["CHROM"], sc["POS"]))
    use_cache = False
    if os.path.exists(cache) and os.path.exists(cache_pos):
        have = pd.read_csv(cache_pos, sep="\t", dtype={"CHROM": str, "POS": int})
        use_cache = set(zip(have["CHROM"], have["POS"])) == want
    if use_cache:
        print(f"Using cached gnomAD extract: {cache}", flush=True)
        gn = pd.read_csv(cache, sep="\t", dtype={"CHROM": str, "REF": str, "ALT": str, "FILTER": str})
    else:
        with tempfile.TemporaryDirectory(dir=args.outdir) as tmpdir:
            gn = extract_gnomad(sc, args.gnomad_dir, args.gnomad_pattern, args.threads, tmpdir)
        gn.to_csv(cache, sep="\t", index=False)
        pd.DataFrame(sorted(want), columns=["CHROM", "POS"]).to_csv(cache_pos, sep="\t", index=False)
        print(f"Cached gnomAD extract: {cache} ({len(gn)} records)", flush=True)
    if not gn.empty:
        gn["POS"] = gn["POS"].astype(int)
        gn["AF_nfe"] = gn["AF_nfe"].map(to_float)
        gn["AF"] = gn["AF"].map(to_float)
        gn["FILTER"] = gn["FILTER"].fillna("")

    thr = args.af_diff_threshold
    records = []
    for _, v in sc.iterrows():
        g = gn[(gn["CHROM"] == v["CHROM"]) & (gn["POS"] == v["POS"])]
        exact = g[(g["REF"] == v["REF"]) & (g["ALT"] == v["ALT"])]
        flip = g[(g["REF"] == v["ALT"]) & (g["ALT"] == v["REF"])]
        others = g[~g.index.isin(exact.index) & ~g.index.isin(flip.index)]

        if v["CHROM"] in missing_files:
            status = "no_gnomad_file"
        elif not exact.empty:
            status = "exact"
        elif not flip.empty:
            status = "orientation_flip_candidate"
        elif not g.empty:
            status = "site_present_allele_absent"
        else:
            status = "absent_in_gnomad"

        rec = dict(row_id=v["row_id"], ID=v["ID"], REF=v["REF"], ALT=v["ALT"],
                   EFFECT_ALLELE=v["EFFECT_ALLELE"], WEIGHT=v["WEIGHT"],
                   canrisk_eaf=v["canrisk_eaf"], match_status=status,
                   n_gnomad_records_at_pos=len(g),
                   other_alleles_at_pos=";".join(f"{r.REF}>{r.ALT}:{r.AF_nfe}" for r in others.itertuples()))

        if not exact.empty:
            e = exact.iloc[0]
            af = e["AF_nfe"]
            eaf = af if v["EFFECT_ALLELE"] == v["ALT"] else 1 - af
            rec.update(gnomad_filter=e["FILTER"], gnomad_AF_nfe=af, gnomad_eaf_nfe=eaf)
        elif not flip.empty:
            f = flip.iloc[0]
            af = f["AF_nfe"]   # frequency of our REF allele (gnomAD ALT)
            eaf = af if v["EFFECT_ALLELE"] == v["REF"] else 1 - af
            rec.update(gnomad_filter=f["FILTER"], gnomad_AF_nfe=af, gnomad_eaf_nfe=eaf)
        else:
            rec.update(gnomad_filter="", gnomad_AF_nfe=np.nan, gnomad_eaf_nfe=np.nan)

        diff = rec["gnomad_eaf_nfe"] - rec["canrisk_eaf"]
        rec["eaf_diff"] = diff
        flag = ""
        if not np.isnan(diff) and abs(diff) > thr:
            flag = "af_discordant"
            if abs((1 - rec["gnomad_eaf_nfe"]) - rec["canrisk_eaf"]) <= thr:
                flag = "af_suggests_effect_allele_swap"
        rec["af_flag"] = flag
        records.append(rec)

    res = pd.DataFrame(records)
    res.to_csv(os.path.join(args.outdir, "site_validation.tsv"), sep="\t", index=False)

    print(f"Variants checked : {len(res)}")
    print("match_status:")
    print(res["match_status"].value_counts().to_string())
    print(f"af_flag (|gnomAD NFE EAF - CanRisk EAF| > {thr}):")
    print(res["af_flag"].replace("", "none").value_counts().to_string())
    not_pass = (res["gnomad_filter"] != "") & (res["gnomad_filter"] != "PASS")
    print(f"gnomAD FILTER != PASS (among matched): {not_pass.sum()}")
    flagged = res[(res["match_status"] != "exact") | (res["af_flag"] != "") | not_pass]
    flagged.to_csv(os.path.join(args.outdir, "site_validation.flagged.tsv"), sep="\t", index=False)
    print(f"Flagged for review: {len(flagged)}  -> site_validation.flagged.tsv")


if __name__ == "__main__":
    main()
