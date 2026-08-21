#!/usr/bin/env python3
"""Plot the headline rho-vs-snapshot-year curve from `data/sweep_results.csv`.

Reads the collected sweep table (`scripts/sweep/collect.py` produces it) and
renders the trend as a PNG at the repo root, so it's the first thing visible
next to the README. Re-run any time `data/sweep_results.csv` changes:

    python scripts/sweep/plot.py
"""

import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
IN_CSV = REPO_ROOT / "data" / "sweep_results.csv"

sys.path.insert(0, str(REPO_ROOT / "scripts" / "pssm_pipeline"))
from config import load_config  # noqa: E402

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

# Palette skill's fixed categorical slots (validated CVD-safe), cycled one per
# DMS assay. A single-assay protein just uses the first color.
ASSAY_COLORS = ["#2a78d6", "#eb6834", "#2e9e5b", "#9a5cd0"]


def plot_by_assay(ax, df, label):
    """One line per DMS assay (each its own color, own CI band); a single
    assay draws one unlabeled line so the legend stays empty."""
    for color, (assay_id, sub) in zip(ASSAY_COLORS, df.groupby("dms_id")):
        sub = sub.sort_values("year")
        ax.fill_between(sub["year"], sub["bootstrap_ci_95_lo"], sub["bootstrap_ci_95_hi"],
                        color=color, alpha=0.12, linewidth=0)
        ax.plot(sub["year"], sub["spearman_rho"], color=color, linewidth=2, zorder=2)
        ax.scatter(sub["year"], sub["spearman_rho"], s=64, color=color,
                   edgecolor="#fcfcfb", linewidth=1.5, zorder=3,
                   label=assay_id if label else None)


def main():
    # This is the rho-vs-year curve at ONE threshold, so it must select a single
    # threshold or the threshold sweep's <year>_t<thr> rows give several rho points
    # per year and the line doubles back. Plot the active config's baseline
    # threshold (the _t0.3 column for spike); the threshold-vs-rho comparison is a
    # separate visualization. Read the baseline from the YAML, ignoring any stray
    # BITSCORE_PER_RESIDUE override so plotting isn't re-thresholded by a leftover.
    os.environ.pop("BITSCORE_PER_RESIDUE", None)
    baseline = load_config()["bitscore_per_residue"]

    df = pd.read_csv(IN_CSV)
    df = df[(df["status"] == "DONE") & (df["protein"] == PROTEIN)
            & ((df["bitscore_per_residue"] - baseline).abs() < 1e-9)].sort_values("year")
    if df.empty:
        raise SystemExit(
            f"no DONE rows for protein={PROTEIN!r} at bitscore_per_residue={baseline} in {IN_CSV}")

    fig, ax = plt.subplots(figsize=(8, 5))

    multi_assay = df["dms_id"].nunique() > 1
    plot_by_assay(ax, df, label=multi_assay)

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
    if multi_assay:
        ax.legend(frameon=False, loc="upper right")

    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=150)
    print(f"Wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
