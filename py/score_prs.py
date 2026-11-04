#!/usr/bin/env python3
"""Score a frozen PRS panel on a cohort VCF of force-genotyped panel sites.

Per sample:
  prs_observed             sum of weight x effect-allele dosage over genotyped sites only
  prs_imputed              missing sites filled with expected dosage 2 x EAF (CanRisk) [primary]
  prs_flagged_as_expected  as prs_imputed, but sites flagged in panel validation set to 2 x EAF
  prs_bam                  score directly from the BAM calls: --force-expected NOT applied
                           (only genuinely missing calls filled with 2 x EAF)
  z / percentile           standardised with CanRisk population mu/sigma (if available)
Optional --force-expected <tsv>: sites (column 'site' as CHROM:POS, or 'ID') whose calls are not
trusted in this batch (e.g. GRCh38 ALT-contig regions with a non ALT-aware alignment). Their dosage
is set to 2 x EAF in every sample and every score; they are counted in n_forced_expected.
Per site: cohort call rate and cohort effect-allele frequency vs CanRisk EAF.
Nothing is dropped silently: missing sites are counted and listed per sample.
"""
import argparse
import os
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm


def alt_dosage(gt):
    a = gt.replace("|", "/").split("/")
    if "." in a:
        return np.nan
    if any(x not in ("0", "1") for x in a):
        sys.exit(f"Unexpected genotype '{gt}' (multi-allelic?)")
    return float(sum(int(x) for x in a))


def read_genotypes(vcf):
    samples = subprocess.run(["bcftools", "query", "-l", vcf],
                             capture_output=True, text=True, check=True).stdout.split()
    out = subprocess.run(["bcftools", "query", "-f", "%CHROM:%POS:%REF:%ALT[\t%GT]\n", vcf],
                         capture_output=True, text=True, check=True).stdout
    ids, rows = [], []
    for line in out.splitlines():
        f = line.split("\t")
        ids.append(f[0])
        rows.append([alt_dosage(g) for g in f[1:]])
    if len(set(ids)) != len(ids):
        sys.exit("Duplicate variant IDs in cohort VCF")
    return samples, pd.DataFrame(rows, index=ids, columns=samples)


def read_force_expected(path, panel):
    if not path:
        return np.zeros(len(panel), dtype=bool), ""
    fe = pd.read_csv(path, sep="\t", dtype=str, comment="#")
    if "ID" in fe.columns:
        keys, by = set(fe["ID"]), panel["ID"]
    elif "site" in fe.columns:
        keys, by = set(fe["site"]), panel["ID"].str.split(":").str[:2].str.join(":")
    else:
        sys.exit("--force-expected file needs a column 'site' (CHROM:POS) or 'ID'")
    mask = by.isin(keys).to_numpy()
    return mask, ";".join(panel.loc[mask, "ID"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--panel-dir", required=True)
    ap.add_argument("--standardisation", default=None)
    ap.add_argument("--force-expected", default=None)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    name = os.path.basename(os.path.normpath(args.panel_dir))
    panel = pd.read_csv(os.path.join(args.panel_dir, "score_plink2.tsv"), sep="\t", dtype=str)
    for c in ["ID", "EFFECT_ALLELE", "WEIGHT", "EAF", "REF", "ALT"]:
        if c not in panel.columns:
            sys.exit(f"Panel score file missing column {c}")
    panel["WEIGHT"] = panel["WEIGHT"].astype(float)
    panel["EAF"] = panel["EAF"].astype(float)
    bad_ea = (panel["EFFECT_ALLELE"] != panel["REF"]) & (panel["EFFECT_ALLELE"] != panel["ALT"])
    if bad_ea.any():
        sys.exit(f"{bad_ea.sum()} panel rows have EFFECT_ALLELE not equal to REF or ALT")

    flag_file = os.path.join(args.panel_dir, "validation", "site_validation.flagged.tsv")
    flagged = set(pd.read_csv(flag_file, sep="\t", dtype=str)["ID"]) if os.path.exists(flag_file) else set()
    forced, forced_ids = read_force_expected(args.force_expected, panel)

    samples, alt = read_genotypes(args.vcf)
    alt = alt.reindex(panel["ID"])   # panel sites absent from VCF become all-missing rows

    w = panel["WEIGHT"].to_numpy()
    exp = 2 * panel["EAF"].to_numpy()
    ea_is_alt = (panel["EFFECT_ALLELE"] == panel["ALT"]).to_numpy()
    A = alt.to_numpy()                                   # n_var x n_samp, NaN = missing
    D = np.where(ea_is_alt[:, None], A, 2 - A)           # effect-allele dosage (observed)
    obs = ~np.isnan(D)
    D_use = np.where(forced[:, None], np.nan, D)         # untrusted calls removed
    obs_use = ~np.isnan(D_use)
    D_imp = np.where(obs_use, D_use, exp[:, None])
    is_flag = panel["ID"].isin(flagged).to_numpy()
    D_flag = np.where(is_flag[:, None], exp[:, None], D_imp)
    D_bam = np.where(obs, D, exp[:, None])               # direct from BAM calls, no forcing

    res = pd.DataFrame({
        "sample": samples,
        "prs": name,
        "n_variants": len(panel),
        "n_genotyped": obs_use.sum(axis=0),
        "n_missing": (~obs_use & ~forced[:, None]).sum(axis=0),
        "n_forced_expected": int(forced.sum()),
        "prs_observed": np.where(obs_use, D_use * w[:, None], 0).sum(axis=0),
        "prs_imputed": (D_imp * w[:, None]).sum(axis=0),
        "prs_flagged_as_expected": (D_flag * w[:, None]).sum(axis=0),
        "prs_bam": (D_bam * w[:, None]).sum(axis=0),
    })

    std = None
    if args.standardisation and os.path.exists(args.standardisation):
        s = pd.read_csv(args.standardisation, sep="\t")
        s = s[s["prs"] == name]
        if len(s) == 1:
            std = s.iloc[0]
    if std is not None:
        res["z_canrisk"] = (res["prs_imputed"] - std["mu"]) / std["sigma"]
        res["percentile_canrisk"] = norm.cdf(res["z_canrisk"]) * 100
        res["z_canrisk_bam"] = (res["prs_bam"] - std["mu"]) / std["sigma"]
        res["percentile_canrisk_bam"] = norm.cdf(res["z_canrisk_bam"]) * 100
    res["z_cohort"] = (res["prs_imputed"] - res["prs_imputed"].mean()) / res["prs_imputed"].std(ddof=1)

    res.to_csv(f"{args.out_prefix}.scores.tsv", sep="\t", index=False, float_format="%.6g")

    miss = [(samples[j], panel["ID"].iat[i]) for i, j in zip(*np.where(~obs & ~forced[:, None]))]
    pd.DataFrame(miss, columns=["sample", "ID"]).to_csv(f"{args.out_prefix}.missing_sites.tsv", sep="\t", index=False)

    n_obs = obs.sum(axis=1)
    cohort_eaf = np.divide(np.nansum(D, axis=1), 2 * n_obs, out=np.full(len(panel), np.nan), where=n_obs > 0)
    site = pd.DataFrame({
        "ID": panel["ID"],
        "flagged_in_validation": is_flag,
        "forced_expected": forced,
        "call_rate": obs.mean(axis=1),
        "cohort_eaf": cohort_eaf,
        "canrisk_eaf": panel["EAF"],
    })
    site["eaf_diff"] = site["cohort_eaf"] - site["canrisk_eaf"]
    site.to_csv(f"{args.out_prefix}.site_summary.tsv", sep="\t", index=False, float_format="%.4g")

    print(f"[{name}] samples={len(samples)} variants={len(panel)} flagged={is_flag.sum()} "
          f"forced_expected={forced.sum()}" + (f" ({forced_ids})" if forced.any() else ""))
    print(f"[{name}] missing per sample (excl. forced): median={np.median(res['n_missing']):.0f} max={res['n_missing'].max()}")
    low = site[(site["call_rate"] < 0.9) & ~site["forced_expected"]]
    print(f"[{name}] non-forced sites with call rate < 0.9: {len(low)}")
    if len(low):
        print(low.to_string(index=False))
    if std is not None:
        print(f"[{name}] z_canrisk corrected: mean={res['z_canrisk'].mean():.3f} sd={res['z_canrisk'].std(ddof=1):.3f}; "
              f"direct from BAM: mean={res['z_canrisk_bam'].mean():.3f} sd={res['z_canrisk_bam'].std(ddof=1):.3f}")


if __name__ == "__main__":
    main()
