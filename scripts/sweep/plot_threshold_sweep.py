#!/usr/bin/env python3
"""Plot the bit-score threshold sweep: rho vs. snapshot year, one line per threshold.

Companion to `plot.py` (which draws the single-threshold headline curve). This
figure overlays every swept bit-score threshold so the effect of hit stringency
can be read off against database growth. Reads the same collected table and
writes one PNG per protein under `plots/`:

    python scripts/sweep/plot_threshold_sweep.py

Lines only, no CI bands: five overlapping bands turn the low-threshold cluster
into an unreadable haze, so the bootstrap CIs live in `data/sweep_results.csv`
instead. Anomalous cells are left in place -- a collapsed point is a validation
finding, not noise to hide.
"""

import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
IN_CSV = REPO_ROOT / "data" / "sweep_results.csv"
PROTEIN = Path(os.environ.get("PROTEIN_CONFIG", "config/spike.yaml")).stem
OUT_PNG = REPO_ROOT / "plots" / f"{PROTEIN}_threshold_sweep.png"
OUT_PNG_SELECTED = REPO_ROOT / "plots" / f"{PROTEIN}_threshold_sweep_everest_selected.png"

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
    for thr in thresholds:
        s = sub[sub["bitscore_per_residue"] == thr].sort_values("year")
        color = THRESHOLD_COLORS.get(thr, "#888888")
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


def plot_assay_selected(ax, sub, title):
    """Same rho-vs-year curves, but each point marked by whether EVEREST's
    alignment-selection floor (Neff/L >= 1.0, `03_weights.py`) would have kept
    it: filled = selected, hollow = rejected as too shallow. Lines are faded to
    context so the selected/rejected split is what reads."""
    for thr in sorted(sub["bitscore_per_residue"].unique()):
        s = sub[sub["bitscore_per_residue"] == thr].sort_values("year")
        color = THRESHOLD_COLORS.get(thr, "#888888")
        ax.plot(s["year"], s["spearman_rho"], color=color, linewidth=1.5,
                alpha=0.3, zorder=2, label=f"{thr:g} bits/res")
        sel = s[s["clears_depth_floor"]]
        rej = s[~s["clears_depth_floor"]]
        ax.scatter(rej["year"], rej["spearman_rho"], s=30, facecolor="none",
                   edgecolor=color, linewidth=1.2, alpha=0.6, zorder=3)
        ax.scatter(sel["year"], sel["spearman_rho"], s=62, color=color,
                   edgecolor="#1a1a1a", linewidth=1.3, zorder=4)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("UniRef100 snapshot year")
    ax.set_xticks(sorted(sub["year"].unique()))
    ax.tick_params(axis="x", rotation=45)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#e1e0d9", linewidth=1, zorder=0)
    ax.set_axisbelow(True)


def _grid(df, assays, plot_fn):
    fig, axes = plt.subplots(1, len(assays), figsize=(6.2 * len(assays), 5),
                             sharey=True)
    axes = [axes] if len(assays) == 1 else axes
    for ax, assay in zip(axes, assays):
        plot_fn(ax, df[df["dms_id"] == assay], assay)
    axes[0].set_ylabel("Spearman's ρ (PSSM vs. DMS)")
    return fig, axes


def main():
    df = pd.read_csv(IN_CSV)
    df = df[(df["status"] == "DONE") & (df["protein"] == PROTEIN)]
    if df.empty:
        raise SystemExit(f"no DONE rows for protein={PROTEIN!r} in {IN_CSV}")

    assays = sorted(df["dms_id"].unique())

    # Figure 1: the plain threshold sweep (lines only).
    fig, axes = _grid(df, assays, plot_assay)
    axes[-1].legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0),
                    title="bit-score threshold")
    fig.suptitle(f"{PROTEIN}: PSSM accuracy vs. snapshot year, by bit-score threshold",
                 fontsize=13)
    fig.tight_layout(rect=(0, 0, 0.9, 1))
    fig.savefig(OUT_PNG, dpi=150)
    print(f"Wrote {OUT_PNG}")

    # Figure 2: same curves, with EVEREST's selection floor applied.
    fig, axes = _grid(df, assays, plot_assay_selected)
    thr_handles = [Line2D([], [], color=THRESHOLD_COLORS[t], lw=2, label=f"{t:g} bits/res")
                   for t in sorted(df["bitscore_per_residue"].unique())]
    sel_handles = [
        Line2D([], [], marker="o", linestyle="none", markersize=8, color="#555",
               markeredgecolor="#1a1a1a", label="selected (Neff/L ≥ 1.0)"),
        Line2D([], [], marker="o", linestyle="none", markersize=6, markerfacecolor="none",
               markeredgecolor="#555", label="rejected (too shallow)"),
    ]
    leg1 = axes[-1].legend(handles=thr_handles, frameon=False, loc="upper left",
                           bbox_to_anchor=(1.02, 1.0), title="bit-score threshold")
    axes[-1].add_artist(leg1)
    axes[-1].legend(handles=sel_handles, frameon=False, loc="upper left",
                    bbox_to_anchor=(1.02, 0.42), title="EVEREST selection")
    fig.suptitle(f"{PROTEIN}: which (year × threshold) alignments EVEREST would select",
                 fontsize=13)
    fig.tight_layout(rect=(0, 0, 0.9, 1))
    fig.savefig(OUT_PNG_SELECTED, dpi=150)
    print(f"Wrote {OUT_PNG_SELECTED}")


if __name__ == "__main__":
    main()
