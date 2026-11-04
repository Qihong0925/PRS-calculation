#!/usr/bin/env python3
"""Step 05b: cohort-level PRS summary tables and figures.

Outputs (<out_prefix>.*):
  cohort_summary.tsv     distribution statistics per panel
  percentile_bands.tsv   observed vs expected (general population) counts per CanRisk percentile band
  prs_distribution.png   z histogram vs population N(0,1)
  percentile_bands.png   observed/expected ratio per percentile band
  panel_concordance.png  primary vs secondary panel z
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm, pearsonr

BANDS = [0, 20, 40, 60, 80, 90, 95, 99, 100.0001]
LABELS = ["<20", "20-40", "40-60", "60-80", "80-90", "90-95", "95-99", "\u226599"]
COLORS = {"primary": "#2a6f97", "secondary": "#e07a3f", "ref": "#7a7a7a"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--final-scores", required=True)
    ap.add_argument("--primary", default="BCAC_313")
    ap.add_argument("--secondary", default="BRIDGES_306")
    ap.add_argument("--batch", required=True)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    d = pd.read_csv(args.final_scores, sep="\t")
    n = len(d)
    panels = [args.primary, args.secondary]
    expected = np.diff(BANDS) / 100 * n

    summ, band_rows = [], []
    for prs in panels:
        z, p = d[f"{prs}_z"], d[f"{prs}_percentile"]
        summ.append(dict(batch=args.batch, prs=prs, n=n, z_mean=z.mean(), z_median=z.median(),
                         z_sd=z.std(ddof=1), z_q25=z.quantile(.25), z_q75=z.quantile(.75),
                         z_min=z.min(), z_max=z.max(), percentile_median=p.median(),
                         n_top10=int((p >= 90).sum()), n_top20=int((p >= 80).sum())))
        obs = pd.cut(p, BANDS, labels=LABELS, right=False).value_counts().reindex(LABELS).values
        for lab, o, e in zip(LABELS, obs, expected):
            band_rows.append(dict(batch=args.batch, prs=prs, band=lab, observed_n=int(o),
                                  observed_pct=o / n * 100, expected_n=e, obs_exp=o / e))
    pd.DataFrame(summ).to_csv(f"{args.out_prefix}.cohort_summary.tsv", sep="\t", index=False, float_format="%.4f")
    bands = pd.DataFrame(band_rows)
    bands.to_csv(f"{args.out_prefix}.percentile_bands.tsv", sep="\t", index=False, float_format="%.4f")

    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})

    # 1. z distribution vs population
    fig, ax = plt.subplots(figsize=(6.5, 4))
    zp = d[f"{args.primary}_z"]
    ax.hist(zp, bins=np.arange(-3.5, 4.01, 0.25), density=True, color=COLORS["primary"], alpha=0.75,
            edgecolor="white", label=f"{args.batch} cases (n={n})")
    x = np.linspace(-4, 4.5, 300)
    ax.plot(x, norm.pdf(x), color=COLORS["ref"], lw=1.8, ls="--", label="General population N(0,1)")
    ax.axvline(0, color=COLORS["ref"], lw=0.8)
    ax.axvline(zp.mean(), color=COLORS["primary"], lw=1.5)
    ax.annotate(f"mean z = {zp.mean():.2f}", xy=(zp.mean(), ax.get_ylim()[1] * 0.92),
                xytext=(6, 0), textcoords="offset points", color=COLORS["primary"])
    ax.set_xlabel(f"{args.primary} PRS (z, CanRisk standardisation)")
    ax.set_ylabel("Density")
    ax.set_title(f"{args.batch}: PRS distribution vs general population", loc="left")
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout(); fig.savefig(f"{args.out_prefix}.prs_distribution.png", dpi=200); plt.close(fig)

    # 2. observed / expected per percentile band
    fig, ax = plt.subplots(figsize=(7, 4))
    xi = np.arange(len(LABELS)); w = 0.38
    for k, (prs, col) in enumerate(zip(panels, [COLORS["primary"], COLORS["secondary"]])):
        b = bands[bands.prs == prs]
        bars = ax.bar(xi + (k - 0.5) * w, b.obs_exp, w, color=col, label=prs)
        for r, v, o in zip(bars, b.obs_exp, b.observed_n):
            ax.text(r.get_x() + r.get_width() / 2, v + 0.05, str(o), ha="center", va="bottom", fontsize=8)
    ax.axhline(1, color=COLORS["ref"], lw=1, ls="--")
    ax.set_xticks(xi, LABELS)
    ax.set_xlabel("CanRisk population percentile band (%)")
    ax.set_ylabel("Observed / expected")
    ax.set_title(f"{args.batch}: enrichment by PRS percentile (numbers = observed n)", loc="left")
    ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(f"{args.out_prefix}.percentile_bands.png", dpi=200); plt.close(fig)

    # 3. panel concordance
    fig, ax = plt.subplots(figsize=(4.5, 4.5))
    za, zb = d[f"{args.primary}_z"], d[f"{args.secondary}_z"]
    ax.scatter(za, zb, s=14, color=COLORS["primary"], alpha=0.8)
    lim = [min(za.min(), zb.min()) - 0.3, max(za.max(), zb.max()) + 0.3]
    ax.plot(lim, lim, color=COLORS["ref"], lw=1, ls="--")
    q90 = norm.ppf(0.9)
    ax.axvline(q90, color=COLORS["ref"], lw=0.6, ls=":"); ax.axhline(q90, color=COLORS["ref"], lw=0.6, ls=":")
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel(f"{args.primary} z"); ax.set_ylabel(f"{args.secondary} z")
    ax.set_title(f"r = {pearsonr(za, zb)[0]:.3f}  (dotted: top 10% threshold)", loc="left", fontsize=10)
    fig.tight_layout(); fig.savefig(f"{args.out_prefix}.panel_concordance.png", dpi=200); plt.close(fig)

    print(pd.DataFrame(summ).round(2).to_string(index=False))
    print(f"\nWritten: {args.out_prefix}.cohort_summary.tsv, .percentile_bands.tsv, "
          f".prs_distribution.png, .percentile_bands.png, .panel_concordance.png")


if __name__ == "__main__":
    main()
