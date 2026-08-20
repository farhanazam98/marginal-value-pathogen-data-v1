#!/usr/bin/env python3
"""Plot the headline rho-vs-snapshot-year curve from `data/sweep_results.csv`.

Reads the collected sweep table (`scripts/sweep/collect.py` produces it) and
renders the trend as a PNG at the repo root, so it's the first thing visible
next to the README. Re-run any time `data/sweep_results.csv` changes:

    python scripts/sweep/plot.py
"""

import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
IN_CSV = REPO_ROOT / "data" / "sweep_results.csv"

# Same PROTEIN_CONFIG env var used to select a protein everywhere else in the
# pipeline (default config/spike.yaml) -- its filename stem picks which
# protein's rows to plot and matches the sandbox naming under $SWEEP_ROOT.
PROTEIN = Path(os.environ.get("PROTEIN_CONFIG", "config/spike.yaml")).stem
# Spike keeps the original, un-prefixed filename the README already links to;
# any other protein gets its own file so the two don't overwrite each other.
OUT_PNG = REPO_ROOT / (
    "pssm_accuracy_vs_snapshot_year.png" if PROTEIN == "spike"
    else f"{PROTEIN}_accuracy_vs_snapshot_year.png"
)

# Tier A (2010-2018) was searched annually; Tier B (2020-2026) picks up
# biennially after a redownload/rerun -- same pipeline, same bit-score threshold
# throughout, so it's one continuous series per assay, not two measurements.
# Palette skill's fixed categorical slots (validated CVD-safe): BLUE/ORANGE mark
# the newly added Tier B years in the single-assay view; when several assays are
# present each gets its own color from ASSAY_COLORS.
BLUE = "#2a78d6"
ORANGE = "#eb6834"
TIER_B_YEARS = {2020, 2022, 2024, 2026}
ASSAY_COLORS = ["#2a78d6", "#eb6834", "#2e9e5b", "#9a5cd0"]


def plot_single_assay(ax, df):
    """One protein, one assay: the original Tier A/B two-color view, where color
    marks the newly added later years rather than a different quantity."""
    ax.fill_between(
        df["year"], df["bootstrap_ci_95_lo"], df["bootstrap_ci_95_hi"],
        color=BLUE, alpha=0.12, linewidth=0,
    )
    ax.plot(df["year"], df["spearman_rho"], color=BLUE, linewidth=2, zorder=2)

    tier_a = df[~df["year"].isin(TIER_B_YEARS)]
    tier_b = df[df["year"].isin(TIER_B_YEARS)]
    ax.scatter(tier_a["year"], tier_a["spearman_rho"], s=64, color=BLUE,
               edgecolor="#fcfcfb", linewidth=1.5, zorder=3,
               label="Tier A -- annual, 2010-2018")
    ax.scatter(tier_b["year"], tier_b["spearman_rho"], s=64, color=ORANGE,
               edgecolor="#fcfcfb", linewidth=1.5, zorder=3,
               label="Tier B -- biennial, 2020-2026 (new)")


def plot_by_assay(ax, df):
    """One line per DMS assay, each its own color, with its own CI band."""
    for color, (assay_id, sub) in zip(ASSAY_COLORS, df.groupby("dms_id")):
        sub = sub.sort_values("year")
        ax.fill_between(sub["year"], sub["bootstrap_ci_95_lo"], sub["bootstrap_ci_95_hi"],
                        color=color, alpha=0.12, linewidth=0)
        ax.plot(sub["year"], sub["spearman_rho"], color=color, linewidth=2, zorder=2)
        ax.scatter(sub["year"], sub["spearman_rho"], s=64, color=color,
                   edgecolor="#fcfcfb", linewidth=1.5, zorder=3, label=assay_id)


def main():
    df = pd.read_csv(IN_CSV)
    df = df[(df["status"] == "DONE") & (df["protein"] == PROTEIN)].sort_values("year")
    if df.empty:
        raise SystemExit(f"no DONE rows for protein={PROTEIN!r} in {IN_CSV}")

    fig, ax = plt.subplots(figsize=(8, 5))

    # Split into a line per assay only when the table actually carries several;
    # the single-assay path reproduces the original committed figure exactly.
    multi_assay = "dms_id" in df.columns and df["dms_id"].nunique() > 1
    if multi_assay:
        plot_by_assay(ax, df)
    else:
        plot_single_assay(ax, df)

    ax.set_xticks(sorted(df["year"].unique()))
    ax.tick_params(axis="x", rotation=0)
    ax.set_xlabel("UniRef100 snapshot year")
    ax.set_ylabel("Spearman's ρ (PSSM vs. DMS)")
    title = "PSSM accuracy vs. UniRef100 snapshot year, 2010–2026"
    if PROTEIN != "spike":
        title = f"{PROTEIN}: {title}"
    ax.set_title(title)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#e1e0d9", linewidth=1, zorder=0)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="upper right")

    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=150)
    print(f"Wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
