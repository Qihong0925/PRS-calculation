#!/usr/bin/env python3
"""Tumour vs germline PRS: allele-specific loss (LOH / allelic imbalance) at PRS sites.

For every tumour and its matched normal, at each panel site:
  - the germline genotype and effect-allele dosage come from the normal (missing calls and
    --force-expected sites use the expected dosage 2 x EAF, as in the germline score);
  - at germline heterozygous sites with enough tumour depth, the tumour effect-allele fraction v
    (from allele depths) is classified using the tumour purity p:
        v >= LOH_VAF                          -> retain_effect  (any purity)
        v <= 1 - LOH_VAF                      -> lose_effect    (any purity)
        1 - NOCHANGE_VAF < v < NOCHANGE_VAF   -> no_change      (any purity)
        NOCHANGE_VAF <= v < LOH_VAF           -> retain_effect if v >= 0.5 + f*p/2, else no_change
        1 - LOH_VAF < v <= 1 - NOCHANGE_VAF   -> lose_effect   if v <= 0.5 - f*p/2, else no_change
    where f = --purity-fraction (fraction of the maximum shift expected from complete LOH, p/2).
  - tumour-cell effect-allele dosage = 2 (retain_effect), 0 (lose_effect), otherwise the germline
    dosage. Homozygous, missing, forced and low-depth sites keep the germline dosage.
Tumour-cell PRS = sum(weight x tumour-cell dosage); the germline PRS uses the germline dosage.
The risk allele at a site is the effect allele if weight > 0, otherwise the other allele.

Outputs (<out_prefix>.*):
  tumour_sites.tsv    one row per tumour x panel x site
  tumour_prs.tsv      one row per tumour x panel: germline, tumour-cell and directly called scores
  tumour_summary.tsv  cohort tests per panel (all tumours, and one tumour per patient)
  tumour_vs_germline.png
"""
import argparse
import os
import sys

import matplotlib
import numpy as np
import pandas as pd
from scipy.stats import binomtest, norm, wilcoxon

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

GT_COLS = ["CHROM", "POS", "REF", "ALT", "GT", "AD"]


def alt_dosage(gt):
    """Alt-allele count from a biallelic GT string; NaN if missing."""
    a = str(gt).replace("|", "/").split("/")
    if "." in a or len(a) != 2:
        return np.nan
    if any(x not in ("0", "1") for x in a):
        sys.exit(f"Unexpected genotype '{gt}' (multi-allelic?)")
    return float(int(a[0]) + int(a[1]))


def parse_ad(ad):
    """Allele depths 'ref,alt' -> (ref, alt); NaN when absent."""
    parts = str(ad).split(",")
    if len(parts) != 2 or "." in parts:
        return np.nan, np.nan
    return float(parts[0]), float(parts[1])


def read_gt(gt_dir, sample):
    path = os.path.join(gt_dir, f"{sample}.gt.tsv")
    if not os.path.exists(path):
        sys.exit(f"Missing genotype table: {path}")
    d = pd.read_csv(path, sep="\t", dtype=str)
    missing = [c for c in GT_COLS if c not in d.columns]
    if missing:
        sys.exit(f"{path} lacks columns {missing}")
    d["ID"] = d["CHROM"] + ":" + d["POS"] + ":" + d["REF"] + ":" + d["ALT"]
    if d["ID"].duplicated().any():
        sys.exit(f"Duplicate sites in {path}")
    ad = d["AD"].map(parse_ad)
    d["ad_ref"] = [x[0] for x in ad]
    d["ad_alt"] = [x[1] for x in ad]
    d["alt_dosage"] = d["GT"].map(alt_dosage)
    return d.set_index("ID")[["GT", "alt_dosage", "ad_ref", "ad_alt"]]


def read_force_expected(path, panel):
    if not path:
        return np.zeros(len(panel), dtype=bool)
    fe = pd.read_csv(path, sep="\t", dtype=str, comment="#")
    if "ID" in fe.columns:
        return panel["ID"].isin(set(fe["ID"])).to_numpy()
    if "site" in fe.columns:
        return panel["ID"].str.split(":").str[:2].str.join(":").isin(set(fe["site"])).to_numpy()
    sys.exit("--force-expected file needs a column 'site' (CHROM:POS) or 'ID'")


def classify(v, p, loh_vaf, nochange_vaf, frac):
    """Classify a germline-heterozygous site from tumour effect-allele fraction v and purity p."""
    hi_purity = 0.5 + frac * p / 2
    lo_purity = 0.5 - frac * p / 2
    if v >= loh_vaf:
        return "retain_effect"
    if v <= 1 - loh_vaf:
        return "lose_effect"
    if 1 - nochange_vaf < v < nochange_vaf:
        return "no_change"
    if v >= nochange_vaf:
        return "retain_effect" if v >= hi_purity else "no_change"
    return "lose_effect" if v <= lo_purity else "no_change"


def effect_fraction(ref, alt, ea_is_alt):
    dp = ref + alt
    with np.errstate(invalid="ignore", divide="ignore"):
        f = np.where(ea_is_alt, alt, ref) / dp
    return dp, f


def sites_for_tumour(row, panel, name, gt_dir, forced, args):
    """Site-level table for one tumour x panel."""
    n = read_gt(gt_dir, row["normal"]).reindex(panel["ID"])
    t = read_gt(gt_dir, row["tumour"]).reindex(panel["ID"])
    w = panel["WEIGHT"].to_numpy()
    eaf = panel["EAF"].to_numpy()
    ea_is_alt = (panel["EFFECT_ALLELE"] == panel["ALT"]).to_numpy()

    g_alt = n["alt_dosage"].to_numpy()
    g_eff = np.where(ea_is_alt, g_alt, 2 - g_alt)
    g_eff_used = np.where(forced | np.isnan(g_eff), 2 * eaf, g_eff)

    n_dp, n_v = effect_fraction(n["ad_ref"].to_numpy(), n["ad_alt"].to_numpy(), ea_is_alt)
    t_dp, t_v = effect_fraction(t["ad_ref"].to_numpy(), t["ad_alt"].to_numpy(), ea_is_alt)
    p = float(row["purity"])

    status, call = [], []
    for i in range(len(panel)):
        if forced[i]:
            s, c = "forced_expected", ""
        elif np.isnan(g_eff[i]):
            s, c = "germline_missing", ""
        elif g_eff[i] != 1:
            s, c = "germline_homozygous", ""
        elif np.isnan(t_dp[i]) or t_dp[i] == 0:
            s, c = "tumour_no_coverage", ""
        elif t_dp[i] < args.min_dp:
            s, c = "tumour_low_depth", ""
        else:
            s = "assessed"
            c = classify(t_v[i], p, args.loh_vaf, args.nochange_vaf, args.purity_fraction)
        status.append(s)
        call.append(c)
    call = np.array(call)
    t_eff = np.where(call == "retain_effect", 2.0, np.where(call == "lose_effect", 0.0, g_eff_used))

    risk_is_effect = w > 0
    retained = np.where(call == "retain_effect", np.where(risk_is_effect, "risk", "protective"),
                        np.where(call == "lose_effect", np.where(risk_is_effect, "protective", "risk"), ""))
    risk_vaf_t = np.where(risk_is_effect, t_v, 1 - t_v)
    risk_vaf_n = np.where(risk_is_effect, n_v, 1 - n_v)

    return pd.DataFrame({
        "tumour": row["tumour"], "normal": row["normal"], "purity": p, "panel": name,
        "ID": panel["ID"].to_numpy(), "weight": w, "effect_allele": panel["EFFECT_ALLELE"].to_numpy(),
        "risk_allele": np.where(risk_is_effect, panel["EFFECT_ALLELE"], panel["OTHER_ALLELE"]),
        "germline_GT": n["GT"].to_numpy(), "germline_effect_dosage": g_eff_used,
        "normal_DP": n_dp, "normal_effect_vaf": n_v,
        "tumour_GT": t["GT"].to_numpy(), "tumour_DP": t_dp, "tumour_effect_vaf": t_v,
        "normal_risk_vaf": risk_vaf_n, "tumour_risk_vaf": risk_vaf_t,
        "status": status, "call": call, "retained_allele": retained,
        "tumour_cell_effect_dosage": t_eff,
    })


def cohort_tests(prs, sites, label):
    """Cohort-level tests on one set of tumours."""
    out = []
    for name, d in prs.groupby("panel", sort=False):
        dz = d["delta_tumour_cell_z"].to_numpy()
        nz = dz[np.abs(dz) > 1e-12]
        s = sites[(sites["panel"] == name) & sites["tumour"].isin(d["tumour"])]
        a = s[s["status"] == "assessed"]
        n_risk = int((a["retained_allele"] == "risk").sum())
        n_prot = int((a["retained_allele"] == "protective").sum())
        shift = (a["tumour_risk_vaf"] - a["normal_risk_vaf"]).dropna().to_numpy()
        shift_nz = shift[np.abs(shift) > 1e-12]
        dcall = d["delta_called_z"].dropna().to_numpy()
        dcall_nz = dcall[np.abs(dcall) > 1e-12]
        out.append(dict(
            set=label, panel=name, n_tumours=len(d),
            mean_delta_tumour_cell_z=dz.mean(), median_delta_tumour_cell_z=np.median(dz),
            n_tumour_higher=int((dz > 1e-12).sum()), n_tumour_lower=int((dz < -1e-12).sum()),
            n_tumour_equal=int((np.abs(dz) <= 1e-12).sum()),
            wilcoxon_p_tumour_cell=wilcoxon(nz).pvalue if len(nz) > 0 else np.nan,
            sign_test_p_tumour_cell=binomtest(int((nz > 0).sum()), len(nz)).pvalue if len(nz) > 0 else np.nan,
            n_sites_assessed=len(a), n_loh_risk_retained=n_risk, n_loh_protective_retained=n_prot,
            binomial_p_risk_vs_protective=binomtest(n_risk, n_risk + n_prot).pvalue if n_risk + n_prot > 0 else np.nan,
            mean_risk_vaf_shift=shift.mean() if len(shift) else np.nan,
            wilcoxon_p_risk_vaf_shift=wilcoxon(shift_nz).pvalue if len(shift_nz) > 0 else np.nan,
            mean_delta_called_z=dcall.mean() if len(dcall) else np.nan,
            wilcoxon_p_called=wilcoxon(dcall_nz).pvalue if len(dcall_nz) > 0 else np.nan,
        ))
    return out


def plot(prs, primary, out):
    d = prs[prs["panel"] == primary]
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.5))
    sc = ax[0].scatter(d["germline_z"], d["tumour_cell_z"], c=d["purity"], cmap="viridis", s=18)
    lim = [np.nanmin(d[["germline_z", "tumour_cell_z"]].to_numpy()) - 0.2,
           np.nanmax(d[["germline_z", "tumour_cell_z"]].to_numpy()) + 0.2]
    ax[0].plot(lim, lim, color="#7a7a7a", lw=1, ls="--")
    ax[0].set_xlabel(f"Germline {primary} z")
    ax[0].set_ylabel(f"Tumour-cell {primary} z (purity-corrected)")
    fig.colorbar(sc, ax=ax[0], label="Purity")
    ax[1].hist(d["delta_tumour_cell_z"], bins=30, color="#2a6f97")
    ax[1].axvline(0, color="#7a7a7a", lw=1, ls="--")
    ax[1].set_xlabel("Tumour-cell z - germline z")
    ax[1].set_ylabel("Tumours")
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sheet", required=True, help="TSV with columns tumour, normal, purity")
    ap.add_argument("--gt-dir", required=True)
    ap.add_argument("--panel-dirs", nargs="+", required=True, help="first = primary")
    ap.add_argument("--germline-scores", required=True, help="score dir of the normals")
    ap.add_argument("--tumour-scores", required=True, help="score dir of the tumours (direct calls)")
    ap.add_argument("--standardisation", default=None)
    ap.add_argument("--force-expected", default=None)
    ap.add_argument("--min-dp", type=int, default=20)
    ap.add_argument("--loh-vaf", type=float, default=0.70)
    ap.add_argument("--nochange-vaf", type=float, default=0.55)
    ap.add_argument("--purity-fraction", type=float, default=0.8)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()
    if not 0.5 < args.nochange_vaf < args.loh_vaf <= 1:
        sys.exit("Need 0.5 < --nochange-vaf < --loh-vaf <= 1")

    sheet = pd.read_csv(args.sheet, sep="\t", dtype={"tumour": str, "normal": str})
    for c in ["tumour", "normal", "purity"]:
        if c not in sheet.columns:
            sys.exit(f"--sheet needs column '{c}'")
    std = None
    if args.standardisation and os.path.exists(args.standardisation):
        std = pd.read_csv(args.standardisation, sep="\t").set_index("prs")

    all_sites, all_prs = [], []
    for pdir in args.panel_dirs:
        name = os.path.basename(os.path.normpath(pdir))
        panel = pd.read_csv(os.path.join(pdir, "score_plink2.tsv"), sep="\t", dtype=str)
        panel["WEIGHT"] = panel["WEIGHT"].astype(float)
        panel["EAF"] = panel["EAF"].astype(float)
        forced = read_force_expected(args.force_expected, panel)
        germ = pd.read_csv(os.path.join(args.germline_scores, f"{name}.scores.tsv"), sep="\t",
                           dtype={"sample": str}).set_index("sample")
        tum = pd.read_csv(os.path.join(args.tumour_scores, f"{name}.scores.tsv"), sep="\t",
                          dtype={"sample": str}).set_index("sample")
        for _, row in sheet.iterrows():
            s = sites_for_tumour(row, panel, name, args.gt_dir, forced, args)
            all_sites.append(s)
            g_raw = float((s["weight"] * s["germline_effect_dosage"]).sum())
            t_raw = float((s["weight"] * s["tumour_cell_effect_dosage"]).sum())
            # Germline score must reproduce the pipeline's germline PRS for this normal
            ref = float(germ.loc[row["normal"], "prs_imputed"])
            if abs(g_raw - ref) > 1e-5 * max(1.0, abs(ref)) + 1e-6:
                sys.exit(f"Germline PRS mismatch for {row['normal']} {name}: {g_raw} vs {ref}")
            a = s[s["status"] == "assessed"]
            all_prs.append(dict(
                tumour=row["tumour"], normal=row["normal"], purity=float(row["purity"]), panel=name,
                germline_raw=g_raw, tumour_cell_raw=t_raw,
                tumour_called_raw=float(tum.loc[row["tumour"], "prs_imputed"]),
                n_sites=len(s), n_germline_het=int((s["germline_GT"].isin(["0/1", "0|1", "1|0"])).sum()),
                n_assessed=len(a),
                n_low_depth=int((s["status"] == "tumour_low_depth").sum()),
                n_no_coverage=int((s["status"] == "tumour_no_coverage").sum()),
                n_retain_effect=int((a["call"] == "retain_effect").sum()),
                n_lose_effect=int((a["call"] == "lose_effect").sum()),
                n_risk_retained=int((a["retained_allele"] == "risk").sum()),
                n_protective_retained=int((a["retained_allele"] == "protective").sum()),
            ))

    sites = pd.concat(all_sites, ignore_index=True)
    prs = pd.DataFrame(all_prs)
    mu = prs["panel"].map(std["mu"]) if std is not None else np.nan
    sigma = prs["panel"].map(std["sigma"]) if std is not None else np.nan
    for col in ["germline", "tumour_cell", "tumour_called"]:
        prs[f"{col}_z"] = (prs[f"{col}_raw"] - mu) / sigma
        prs[f"{col}_percentile"] = norm.cdf(prs[f"{col}_z"]) * 100
    prs["delta_tumour_cell_z"] = prs["tumour_cell_z"] - prs["germline_z"]
    prs["delta_called_z"] = prs["tumour_called_z"] - prs["germline_z"]

    # One tumour per patient (highest purity) for a sensitivity check of non-independence
    first = (prs.drop_duplicates("tumour").sort_values(["normal", "purity", "tumour"], ascending=[True, False, True])
             .drop_duplicates("normal")["tumour"])
    summ = cohort_tests(prs, sites, "all_tumours") + \
        cohort_tests(prs[prs["tumour"].isin(first)], sites, "one_tumour_per_patient")

    sites.to_csv(f"{args.out_prefix}.tumour_sites.tsv", sep="\t", index=False, float_format="%.4g")
    prs.to_csv(f"{args.out_prefix}.tumour_prs.tsv", sep="\t", index=False, float_format="%.6g")
    pd.DataFrame(summ).to_csv(f"{args.out_prefix}.tumour_summary.tsv", sep="\t", index=False, float_format="%.4g")
    primary = os.path.basename(os.path.normpath(args.panel_dirs[0]))
    plot(prs, primary, f"{args.out_prefix}.tumour_vs_germline.png")
    print(f"Tumours: {sheet.shape[0]}; panels: {prs['panel'].nunique()}; "
          f"LOH calls (risk/protective retained): {int(prs['n_risk_retained'].sum())}/{int(prs['n_protective_retained'].sum())}")


if __name__ == "__main__":
    main()
