#!/usr/bin/env python3
"""Step 05: sensitivity analysis and final per-sample PRS table for one batch.

Scores compared (z on the CanRisk scale):
  BAM        primary panel, directly from BAM calls (no batch-specific corrections)   prs_bam
  corrected  primary panel, batch corrections applied (--force-expected)  [primary]   prs_imputed
  flagged    corrected, plus validation-flagged sites set to 2xEAF                     prs_flagged_as_expected
  secondary  secondary panel (e.g. BRIDGES_306), corrected                             prs_imputed
Reports Pearson/Spearman correlation, max |dz|, and concordance of top-decile / top-quintile
membership. Writes a merged final table with both the direct-from-BAM and the corrected scores.
"""
import argparse

import numpy as np
import pandas as pd
from scipy.stats import norm, pearsonr, spearmanr


def z(raw, mu, sigma):
    return (raw - mu) / sigma


def compare(name, a, b, out):
    dz = b - a
    pa, pb = norm.cdf(a) * 100, norm.cdf(b) * 100
    row = dict(comparison=name, pearson=pearsonr(a, b)[0], spearman=spearmanr(a, b)[0],
               mean_dz=dz.mean(), max_abs_dz=np.abs(dz).max(),
               max_abs_dpercentile=np.abs(pb - pa).max())
    for thr, lab in [(90, "top10"), (80, "top20")]:
        ia, ib = pa >= thr, pb >= thr
        row[f"{lab}_n_A"], row[f"{lab}_n_B"] = int(ia.sum()), int(ib.sum())
        row[f"{lab}_discordant"] = int((ia != ib).sum())
    out.append(row)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--score-dir", required=True)
    ap.add_argument("--primary", default="BCAC_313")
    ap.add_argument("--secondary", default="BRIDGES_306")
    ap.add_argument("--standardisation", required=True)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    std = pd.read_csv(args.standardisation, sep="\t").set_index("prs")
    p = pd.read_csv(f"{args.score_dir}/{args.primary}.scores.tsv", sep="\t").set_index("sample")
    s = pd.read_csv(f"{args.score_dir}/{args.secondary}.scores.tsv", sep="\t").set_index("sample")
    s = s.reindex(p.index)
    if s.isna().all(axis=1).any():
        raise SystemExit("Sample sets differ between panels")
    for df, nm in [(p, args.primary), (s, args.secondary)]:
        if "prs_bam" not in df.columns:
            raise SystemExit(f"{nm}.scores.tsv has no prs_bam column; re-run step 04 with the current score_prs.py")

    mp, sp = std.loc[args.primary], std.loc[args.secondary]
    zA = z(p["prs_imputed"], mp["mu"], mp["sigma"])
    zA_bam = z(p["prs_bam"], mp["mu"], mp["sigma"])
    zB = z(p["prs_flagged_as_expected"], mp["mu"], mp["sigma"])
    zC = z(s["prs_imputed"], sp["mu"], sp["sigma"])
    zC_bam = z(s["prs_bam"], sp["mu"], sp["sigma"])
    zD = z(s["prs_flagged_as_expected"], sp["mu"], sp["sigma"])

    res = []
    compare(f"{args.primary}: direct-from-BAM vs corrected", zA_bam, zA, res)
    compare(f"{args.secondary}: direct-from-BAM vs corrected", zC_bam, zC, res)
    compare(f"{args.primary}: corrected vs flagged_as_expected", zA, zB, res)
    compare(f"{args.primary} vs {args.secondary} (corrected)", zA, zC, res)
    compare(f"{args.secondary}: corrected vs flagged_as_expected", zC, zD, res)
    res = pd.DataFrame(res)
    res.to_csv(f"{args.out_prefix}.sensitivity.tsv", sep="\t", index=False, float_format="%.4f")

    P, S = args.primary, args.secondary
    final = pd.DataFrame({
        "sample": p.index,
        f"{P}_n_genotyped": p["n_genotyped"].values,
        f"{P}_n_missing": p["n_missing"].values,
        f"{P}_n_forced_expected": p["n_forced_expected"].values,
        f"{P}_raw_bam": p["prs_bam"].values,
        f"{P}_z_bam": zA_bam.values,
        f"{P}_percentile_bam": norm.cdf(zA_bam.values) * 100,
        f"{P}_raw": p["prs_imputed"].values,
        f"{P}_z": zA.values,
        f"{P}_percentile": norm.cdf(zA.values) * 100,
        f"{P}_z_flagged_as_expected": zB.values,
        f"{S}_n_missing": s["n_missing"].values,
        f"{S}_n_forced_expected": s["n_forced_expected"].values,
        f"{S}_raw_bam": s["prs_bam"].values,
        f"{S}_z_bam": zC_bam.values,
        f"{S}_percentile_bam": norm.cdf(zC_bam.values) * 100,
        f"{S}_raw": s["prs_imputed"].values,
        f"{S}_z": zC.values,
        f"{S}_percentile": norm.cdf(zC.values) * 100,
    })
    final.to_csv(f"{args.out_prefix}.final_scores.tsv", sep="\t", index=False, float_format="%.5g")

    pd.set_option("display.width", 220)
    print(res.to_string(index=False))
    print(f"\nCohort mean z, corrected:        {P} {zA.mean():.3f} (sd {zA.std(ddof=1):.3f}); "
          f"{S} {zC.mean():.3f} (sd {zC.std(ddof=1):.3f})")
    print(f"Cohort mean z, direct from BAM:  {P} {zA_bam.mean():.3f} (sd {zA_bam.std(ddof=1):.3f}); "
          f"{S} {zC_bam.mean():.3f} (sd {zC_bam.std(ddof=1):.3f})")
    print(f"Samples in CanRisk top 10% ({P}, corrected): {(norm.cdf(zA) >= 0.9).sum()} / {len(zA)} "
          f"(expected {0.1*len(zA):.1f} in the general population)")
    print("\nfinal_scores.tsv: *_bam = direct from BAM calls; *_raw/_z/_percentile = corrected (primary)")


if __name__ == "__main__":
    main()
