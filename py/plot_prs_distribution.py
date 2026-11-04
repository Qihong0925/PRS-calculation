#!/usr/bin/env python3
"""Plot cohort PRS z distribution (KDE curves) against the general-population N(0,1).

Shows: population N(0,1) curve; cohort KDE per panel; the population top-10% threshold (z = 1.28);
the cohort area above the threshold shaded (primary panel) and labelled with the observed
proportion of samples per panel. The legend (with N(mean, SD^2) per curve) sits below the plot.

Usage: plot_prs_distribution.py <final_scores.tsv> <batch> <out_prefix> [panel ...]
Default panels: BCAC_313 BRIDGES_306 (columns <panel>_z in final_scores.tsv); the first is primary.
Outputs <out_prefix>.prs_distribution.png (600 dpi) and .pdf.
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.stats import gaussian_kde, norm

final_scores, batch, out_prefix = sys.argv[1:4]
panels = sys.argv[4:] or ["BCAC_313", "BRIDGES_306"]
colors = ["#2a6f97", "#e07a3f", "#4f9d69", "#8e5ea2"]
ref_col = "#7a7a7a"
THRESH_LS = (0, (6, 2, 1, 2))

d = pd.read_csv(final_scores, sep="\t")
n = len(d)
x = np.linspace(-4, 4.5, 600)
q90 = norm.ppf(0.90)
tail = x >= q90
label = lambda p: p.replace("_", " ")

plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
fig, ax = plt.subplots(figsize=(7.5, 5.6))

# Curves
y_pop = norm.pdf(x)
ax.plot(x, y_pop, color=ref_col, lw=1.8, ls="--")
stats = []
ymax = y_pop.max()
for k, p in enumerate(panels):
    z = d[f"{p}_z"].dropna()
    y = gaussian_kde(z)(x)
    ymax = max(ymax, y.max())
    ax.plot(x, y, color=colors[k], lw=2.2)
    stats.append(dict(panel=p, mean=z.mean(), sd=z.std(ddof=1), top10=(z >= q90).mean() * 100,
                      n_top=int((z >= q90).sum()), color=colors[k], y=y))

# Shade the primary panel's area above the threshold; threshold line
ax.fill_between(x[tail], stats[0]["y"][tail], color=stats[0]["color"], alpha=0.25, lw=0)
ax.axvline(q90, color=ref_col, lw=1, ls=THRESH_LS)

# Percentage labels inside the plot, next to the shaded tail
tx, ty, dy = 2.35, ymax * 0.95, ymax * 0.075
ax.text(tx, ty, "Above threshold:", color="#333333", fontsize=8, va="top", fontweight="bold")
for k, s in enumerate(stats):
    ax.text(tx, ty - (k + 1) * dy, f"{label(s['panel'])}: {s['top10']:.0f}% ({s['n_top']}/{n})",
            color=s["color"], fontsize=8, va="top")
ax.text(tx, ty - (len(stats) + 1) * dy, "General population: 10%", color="#4d4d4d", fontsize=8, va="top")

# Rug: individual samples (primary panel)
z0 = d[f"{panels[0]}_z"]
ax.plot(z0, np.full(n, -0.014), "|", color=colors[0], ms=6, alpha=0.6)

ax.set_ylim(-0.028, ymax * 1.08)
ax.set_xlim(x.min(), x.max())
ax.set_xlabel("PRS (z, CanRisk population standardisation)")
ax.set_ylabel("Density")
ax.set_title(f"{batch} cases (n = {n}): PRS distribution vs general population", loc="left")

# Legend below the plot
handles = [Line2D([], [], color=ref_col, lw=1.8, ls="--", label="General population, N(0, 1)")]
handles += [Line2D([], [], color=s["color"], lw=2.2,
                   label=f"{batch}, {label(s['panel'])}, N({s['mean']:.2f}, {s['sd']:.2f}²)") for s in stats]
handles += [Line2D([], [], color=ref_col, lw=1, ls=THRESH_LS,
                   label=f"Population top-10% threshold (z = {q90:.2f}, both panels)"),
            Patch(facecolor=stats[0]["color"], alpha=0.25,
                  label=f"{batch}, {label(panels[0])}: area above threshold")]
fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.06, 0.035), ncol=2, frameon=False,
           fontsize=8, columnspacing=1.5, handlelength=2.5)
fig.text(0.07, 0.005, "N(mean, SD²). Cohort curves are kernel density estimates; "
         f"ticks below the axis are individual samples ({label(panels[0])}).",
         fontsize=7, color="#4d4d4d", va="bottom")

fig.tight_layout(rect=(0, 0.19, 1, 1))
fig.savefig(f"{out_prefix}.prs_distribution.png", dpi=600)
fig.savefig(f"{out_prefix}.prs_distribution.pdf")
print(f"Written: {out_prefix}.prs_distribution.png (600 dpi) / .pdf")
for s in stats:
    print(f"  {s['panel']}: mean z {s['mean']:.3f}, SD {s['sd']:.3f}, above threshold {s['top10']:.1f}% ({s['n_top']}/{n})")
