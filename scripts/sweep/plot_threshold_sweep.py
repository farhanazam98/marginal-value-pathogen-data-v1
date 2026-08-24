#!/usr/bin/env python3
"""Plot the bit-score threshold sweep: rho vs. snapshot year, one line per threshold.

Companion to `plot.py` (which draws the single-threshold headline curve). This
figure overlays every swept bit-score threshold so the effect of hit stringency
can be read off against database growth. Reads the same collected table and
writes one PNG per protein under `plots/`:

    python scripts/sweep/plot_threshold_sweep.py

Only the config baseline threshold has CI bands (baseline vs. the best-scoring
threshold) so the "is the improvement real?" comparison stays legible instead of
five overlapping bands. Anomalous cells are left in place -- a collapsed point is
a validation finding, not noise to hide.
"""

import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
IN_CSV = REPO_ROOT / "data" / "sweep_results.csv"
PROTEIN = Path(os.environ.get("PROTEIN_CONFIG", "config/spike.yaml")).stem
OUT_PNG = REPO_ROOT / "plots" / f"{PROTEIN}_threshold_sweep.png"

# Sequential (ColorBrewer Blues, darkened): threshold is ordinal, so a single-hue
# ramp reads "looser -> stricter" as light -> dark. Keyed by bits/residue.
THRESHOLD_COLORS = {
    0.1: "#9ecae1",
    0.2: "#6baed6",
    0.3: "#4292c6",
    0.4: "#2171b5",
    0.5: "#08306b",
}


def plot_assay(ax, sub, title):
    thresholds = sorted(sub["bitscore_per_residue"].unique())
    baseline = min(thresholds, key=lambda t: abs(t - 0.3))  # config default for spike
    best = max(thresholds)  # 0.5 is the standout column in the data
    for thr in thresholds:
        s = sub[sub["bitscore_per_residue"] == thr].sort_values("year")
        color = THRESHOLD_COLORS.get(thr, "#888888")
        if thr in (baseline, best):
            ax.fill_between(s["year"], s["bootstrap_ci_95_lo"], s["bootstrap_ci_95_hi"],
                            color=color, alpha=0.15, linewidth=0, zorder=1)
        ax.plot(s["year"], s["spearman_rho"], color=color, linewidth=2, zorder=2,
                label=f"{thr:g} bits/res")
        ax.scatter(s["year"], s["spearman_rho"], s=42, color=color,
                   edgecolor="#fcfcfb", linewidth=1.2, zorder=3)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("UniRef100 snapshot year")
    ax.set_xticks(sorted(sub["year"].unique()))
    ax.tick_params(axis="x", rotation=45)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#e1e0d9", linewidth=1, zorder=0)
    ax.set_axisbelow(True)


def main():
    df = pd.read_csv(IN_CSV)
    df = df[(df["status"] == "DONE") & (df["protein"] == PROTEIN)]
    if df.empty:
        raise SystemExit(f"no DONE rows for protein={PROTEIN!r} in {IN_CSV}")

    assays = sorted(df["dms_id"].unique())
    fig, axes = plt.subplots(1, len(assays), figsize=(6.2 * len(assays), 5),
                             sharey=True)
    axes = [axes] if len(assays) == 1 else axes
    for ax, assay in zip(axes, assays):
        plot_assay(ax, df[df["dms_id"] == assay], assay)
    axes[0].set_ylabel("Spearman's ρ (PSSM vs. DMS)")
    axes[-1].legend(frameon=False, loc="upper right", title="bit-score threshold")

    fig.suptitle(f"{PROTEIN}: PSSM accuracy vs. snapshot year, by bit-score threshold",
                 fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=150)
    print(f"Wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
