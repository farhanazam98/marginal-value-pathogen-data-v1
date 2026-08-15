#!/usr/bin/env python3
"""Plot the headline rho-vs-snapshot-year curve from `data/sweep_results.csv`.

Reads the collected sweep table (`scripts/sweep/collect.py` produces it) and
renders the trend as a PNG at the repo root, so it's the first thing visible
next to the README. Re-run any time `data/sweep_results.csv` changes:

    python scripts/sweep/plot.py
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
IN_CSV = REPO_ROOT / "data" / "sweep_results.csv"
OUT_PNG = REPO_ROOT / "pssm_accuracy_vs_snapshot_year.png"

# Tier A (2010-2018) was searched annually; Tier B (2020-2026) picks up
# biennially after a redownload/rerun -- same pipeline, same DMS, same
# bit-score threshold throughout, so it's one continuous series, not two
# different measurements. Colors are the palette skill's fixed categorical
# slots 1 and 2 (validated CVD-safe adjacent pair), used only to mark which
# years are the newly added ones -- not to imply a different quantity.
BLUE = "#2a78d6"
ORANGE = "#eb6834"
TIER_B_YEARS = {2020, 2022, 2024, 2026}


def main():
    df = pd.read_csv(IN_CSV)
    df = df[df["status"] == "DONE"].sort_values("year")

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.fill_between(
        df["year"], df["bootstrap_ci_95_lo"], df["bootstrap_ci_95_hi"],
        color=BLUE, alpha=0.12, linewidth=0,
    )
    ax.plot(df["year"], df["spearman_rho"], color=BLUE, linewidth=2, zorder=2)

    tier_a = df[~df["year"].isin(TIER_B_YEARS)]
    tier_b = df[df["year"].isin(TIER_B_YEARS)]
    ax.scatter(
        tier_a["year"], tier_a["spearman_rho"],
        s=64, color=BLUE, edgecolor="#fcfcfb", linewidth=1.5, zorder=3,
        label="Tier A -- annual, 2010-2018",
    )
    ax.scatter(
        tier_b["year"], tier_b["spearman_rho"],
        s=64, color=ORANGE, edgecolor="#fcfcfb", linewidth=1.5, zorder=3,
        label="Tier B -- biennial, 2020-2026 (new)",
    )

    ax.set_xticks(sorted(df["year"].unique()))
    ax.tick_params(axis="x", rotation=0)
    ax.set_xlabel("UniRef100 snapshot year")
    ax.set_ylabel("Spearman's ρ (PSSM vs. DMS)")
    ax.set_title("PSSM accuracy vs. UniRef100 snapshot year, 2010–2026")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#e1e0d9", linewidth=1, zorder=0)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="upper right")

    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=150)
    print(f"Wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
