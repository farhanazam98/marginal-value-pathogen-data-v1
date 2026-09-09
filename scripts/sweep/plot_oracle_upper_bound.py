#!/usr/bin/env python3
"""Plot the theoretical upper bound: max Spearman rho over every swept bit-score
threshold, per (protein, DMS assay, year).

This is a hindsight-optimal ceiling ("best achievable rho if the ideal threshold
were picked after seeing the DMS labels"), not a deployable selector's output --
no single-threshold selector (e.g. the within-protein EVEREST pick in
`plot_threshold_sweep.py`) is used, so the `Neff_at_90pct_identity` formulation
bug that affects that selector is irrelevant here. Every axis/caption in this
script says "theoretical upper bound" for that reason; don't let a plain "PSSM
accuracy" reading of these lines stand in for pipeline performance elsewhere.

Covers all five proteins -- protease, hiv_env, dengue_polg, flu_h1_ha and spike
-- each multi-assay protein (flu_h1_ha, spike) plotted as two separate lines,
since its two DMS assays measure different phenotypes on a shared alignment. Run:

    python scripts/sweep/plot_oracle_upper_bound.py
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
IN_CSV = REPO_ROOT / "data" / "sweep_results.csv"
OUT_PNG = REPO_ROOT / "plots" / "oracle_threshold_upper_bound.png"
OUT_PNG_BASELINE = REPO_ROOT / "plots" / "oracle_threshold_upper_bound_everest_baseline.png"

# (protein, dms_id) -> (legend label, color, linestyle, marker). A multi-assay
# protein's two assays share a hue and are told apart by line style/marker
# instead of a second color -- a light tint of either hue fails the palette
# validator's lightness/chroma checks.
SERIES = {
    ("dengue_polg", "suphatrakul_fitness"): ("dengue_polg (Suphatrakul)", "#2a78d6", "-", "o"),
    ("hiv_env", "haddox_fitness"): ("hiv_env (Haddox)", "#eb6834", "-", "o"),
    ("protease", "flynn_fitness"): ("protease (Flynn)", "#2e9e5b", "-", "o"),
    ("flu_h1_ha", "doud_fitness"): ("flu_h1_ha (Doud)", "#9a5cd0", "-", "o"),
    ("flu_h1_ha", "wu_fitness"): ("flu_h1_ha (Wu)", "#9a5cd0", "--", "^"),
    ("spike", "starr_binding"): ("spike (Starr binding)", "#1baf7a", "-", "o"),
    ("spike", "starr_expression"): ("spike (Starr expression)", "#1baf7a", "--", "^"),
}

# EVEREST's own PSSM_Spearman baseline per assay (their alignment/threshold
# choice, not tied to a UniRef year) -- from
# github.com/debbiemarkslab/priority-viruses results/summary/hybrid_summary.csv.
# Confirmed with the user before use (2026-09-09).
EVEREST_BASELINE = {
    ("dengue_polg", "suphatrakul_fitness"): 0.4268238377770937,
    ("hiv_env", "haddox_fitness"): 0.4858169201621786,
    ("protease", "flynn_fitness"): 0.5804935452338189,
    ("flu_h1_ha", "doud_fitness"): 0.3979960234111402,
    ("flu_h1_ha", "wu_fitness"): 0.3795681525987909,
    ("spike", "starr_binding"): 0.150197432112948,
    ("spike", "starr_expression"): 0.1785275307513163,
}

def oracle_max(df):
    """Max spearman_rho per (protein, dms_id, year) across all swept thresholds."""
    return df.groupby(["protein", "dms_id", "year"])["spearman_rho"].max() \
        .rename("oracle_rho").reset_index()


def draw(ax, oracle, with_baseline):
    for key, (label, color, ls, marker) in SERIES.items():
        s = oracle[(oracle["protein"] == key[0]) & (oracle["dms_id"] == key[1])].sort_values("year")
        ax.plot(s["year"], s["oracle_rho"], color=color, linestyle=ls, linewidth=2,
                zorder=2, label=label)
        ax.scatter(s["year"], s["oracle_rho"], s=50, color=color, marker=marker,
                   edgecolor="#fcfcfb", linewidth=1.2, zorder=3)
        if with_baseline:
            years = s["year"]
            if len(years):
                ax.hlines(EVEREST_BASELINE[key], years.min(), years.max(),
                          color=color, linestyle=ls, linewidth=1, alpha=0.45, zorder=1)

    ax.set_xlabel("UniRef100 snapshot year")
    ax.set_ylabel("Spearman's ρ — theoretical upper bound")
    ax.set_xticks(sorted(oracle["year"].unique()))
    ax.tick_params(axis="x", rotation=45)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#e1e0d9", linewidth=1, zorder=0)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0),
              title="protein (assay)")


def main():
    df = pd.read_csv(IN_CSV)
    df = df[df["status"] == "DONE"].copy()

    oracle = oracle_max(df)

    fig, ax = plt.subplots(figsize=(10, 5.5))
    draw(ax, oracle, with_baseline=False)
    fig.suptitle("Max Spearman ρ over all swept bit-score thresholds, per protein-year",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 0.82, 0.95))
    fig.savefig(OUT_PNG, dpi=150)
    print(f"Wrote {OUT_PNG}")

    fig, ax = plt.subplots(figsize=(10, 5.5))
    draw(ax, oracle, with_baseline=True)
    fig.suptitle("Max Spearman ρ vs. EVEREST PSSM baseline (flat lines; dashed = second assay)",
                fontsize=12)
    fig.tight_layout(rect=(0, 0, 0.82, 0.95))
    fig.savefig(OUT_PNG_BASELINE, dpi=150)
    print(f"Wrote {OUT_PNG_BASELINE}")

    return oracle


if __name__ == "__main__":
    main()
