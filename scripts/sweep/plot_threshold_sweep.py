#!/usr/bin/env python3
"""Plot the bit-score threshold sweep: rho vs. snapshot year, one line per threshold.

Companion to `plot.py` (which draws the single-threshold headline curve). Writes
two PNGs per protein under `plots/`: `*_threshold_sweep.png` overlays every swept
bit-score threshold so the effect of hit stringency can be read off against
database growth, and `*_threshold_sweep_everest_pick.png` traces the single
alignment EVEREST's DMS-blind rule would select each year. Run:

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
OUT_PNG_PICK = REPO_ROOT / "plots" / f"{PROTEIN}_threshold_sweep_everest_pick.png"

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


def everest_picks(sub):
    """The single alignment EVEREST would run per year, selected DMS-blind by the
    paper's within-protein rule (Methods A.6.1): among alignments that clear the
    depth floor (Neff/L > 1), the one with the highest *proportion* of sequences
    within 90% identity of the query (prop90 = Neff@90%ID / Neff). The pick never
    looks at rho. Note this is NOT max-depth: the paper falls back to max Neff/L
    only for proteins where no alignment reaches the floor. The ``Neff@90%ID >= 30``
    reliability bar is deliberately not applied here -- that is EVEREST's
    cross-protein *confidence* flag, not part of within-protein selection. Years
    where nothing clears the floor are absent (EVEREST would decline to model them)."""
    elig = sub[sub["clears_depth_floor"]]
    picks = [g.loc[g["prop90"].idxmax()] for _, g in elig.groupby("year")]
    return pd.DataFrame(picks).sort_values("year")


def plot_assay_pick(ax, sub, title, label_thr=False):
    for thr in sorted(sub["bitscore_per_residue"].unique()):
        s = sub[sub["bitscore_per_residue"] == thr].sort_values("year")
        ax.plot(s["year"], s["spearman_rho"], color=THRESHOLD_COLORS.get(thr, "#888"),
                linewidth=1.2, alpha=0.22, zorder=2)
    p = everest_picks(sub)
    ax.plot(p["year"], p["spearman_rho"], color="#1a1a1a", linewidth=2.2, zorder=4)
    for _, r in p.iterrows():
        ax.scatter(r["year"], r["spearman_rho"], s=72,
                   color=THRESHOLD_COLORS.get(r["bitscore_per_residue"], "#888"),
                   edgecolor="#1a1a1a", linewidth=1.4, zorder=5)
        if label_thr:
            ax.annotate(f"{r['bitscore_per_residue']:g}", (r["year"], r["spearman_rho"]),
                        textcoords="offset points", xytext=(0, 9), ha="center",
                        fontsize=7.5, color="#444")
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
    df = df[(df["status"] == "DONE") & (df["protein"] == PROTEIN)].copy()
    if df.empty:
        raise SystemExit(f"no DONE rows for protein={PROTEIN!r} in {IN_CSV}")

    # Alignment relevance: fraction of the alignment's effective sequences that
    # sit within 90% identity of the query. EVEREST's within-protein selector.
    df["prop90"] = df["Neff_at_90pct_identity"] / df["Neff"]

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

    # Figure 2: the single alignment EVEREST would select per year, traced as one
    # bold curve; thresholds still colour-code the markers so the shifting pick is
    # visible.
    fig, axes = plt.subplots(1, len(assays), figsize=(6.2 * len(assays), 5), sharey=True)
    axes = [axes] if len(assays) == 1 else axes
    for i, (ax, assay) in enumerate(zip(axes, assays)):
        plot_assay_pick(ax, df[df["dms_id"] == assay], assay, label_thr=(i == 0))
    axes[0].set_ylabel("Spearman's ρ (PSSM vs. DMS)")
    pick_handles = [Line2D([], [], color="#1a1a1a", lw=2.2, marker="o",
                           markeredgecolor="#1a1a1a", markerfacecolor="#888", markersize=8,
                           label="EVEREST pick")]
    thr_handles = [Line2D([], [], marker="o", linestyle="none", markersize=8,
                          color=THRESHOLD_COLORS[t], markeredgecolor="#1a1a1a",
                          label=f"{t:g} bits/res")
                   for t in sorted(df["bitscore_per_residue"].unique())]
    leg1 = axes[-1].legend(handles=pick_handles, frameon=False, loc="upper left",
                           bbox_to_anchor=(1.02, 1.0))
    axes[-1].add_artist(leg1)
    axes[-1].legend(handles=thr_handles, frameon=False, loc="upper left",
                    bbox_to_anchor=(1.02, 0.82), title="picked threshold")
    fig.suptitle(f"{PROTEIN}: the alignment EVEREST would select each year "
                 "— Neff/L > 1, then max fraction ≥90% ID (DMS-blind; "
                 "2010–2011 unmodellable)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 0.9, 1))
    fig.savefig(OUT_PNG_PICK, dpi=150)
    print(f"Wrote {OUT_PNG_PICK}")


if __name__ == "__main__":
    main()
