#!/usr/bin/env python3
"""Step 6: join predictions back to the DMS and compute Spearman rho.
Writes data/pssm_pipeline/scatter.png.

Spearman rho is Pearson correlation computed on ranks rather than raw
values: it asks whether the two orderings (by predicted_score, by
DMS_score) move together, ignoring magnitude entirely. That's the right
tool here because our log-odds scores have no claim to being on the same
scale as the DMS's binding measurement -- only their relative order is
meant to carry information. Ties (e.g. Step 5's pseudocount-floor ties, many
variants at a column sharing an identical score) are handled by assigning
the tied variants the average of the ranks they span, scipy's default.
"""

import json
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

DMS_FILE = "data/SARS2_RBD_Starr_binding_dms.csv"
PREDICTIONS = "data/pssm_pipeline/predictions.csv"
OUT_SCATTER = "data/pssm_pipeline/scatter.png"
OUT_META = "data/pssm_pipeline/evaluate_meta.json"

MUTANT_RE = re.compile(r"^([A-Z])(\d+)([A-Z])$")
N_BOOTSTRAP = 10_000
SEED = 0


def main():
    dms = pd.read_csv(DMS_FILE)[["mutant", "DMS_score"]]
    parsed = dms["mutant"].str.extract(MUTANT_RE)
    parsed.columns = ["wt_aa", "position", "mut_aa"]
    parsed["position"] = parsed["position"].astype(int)
    dms = pd.concat([dms, parsed], axis=1)
    print(f"DMS variants: {len(dms)}")

    predictions = pd.read_csv(PREDICTIONS)
    print(f"Prediction rows: {len(predictions)}")

    # --- Inner join on (position, wt_aa, mut_aa), independent of Step 5's own bookkeeping ---
    joined = pd.merge(
        dms,
        predictions[["position", "wt_aa", "mut_aa", "predicted_score", "imputed"]],
        on=["position", "wt_aa", "mut_aa"],
        how="inner",
    )
    n_dropped_from_dms = len(dms) - len(joined)
    n_dropped_from_predictions = len(predictions) - len(joined)
    print(f"\nJoined: {len(joined)} variants")
    print(f"DMS variants dropped by the join: {n_dropped_from_dms} "
          f"(no matching prediction row for that position/wt/mut)")
    print(f"Prediction rows dropped by the join: {n_dropped_from_predictions} "
          f"(no matching DMS row -- shouldn't happen, predictions.csv was built from this same DMS file)")

    # ---------------- Spearman rho ----------------
    rho, pval = spearmanr(joined["predicted_score"], joined["DMS_score"])
    print(f"\nSpearman rho = {rho:.4f} (p={pval:.2e}), n={len(joined)}")

    rng = np.random.default_rng(SEED)
    n = len(joined)
    pred_vals = joined["predicted_score"].to_numpy()
    dms_vals = joined["DMS_score"].to_numpy()
    boot_rhos = np.empty(N_BOOTSTRAP)
    for b in range(N_BOOTSTRAP):
        idx = rng.integers(0, n, size=n)
        boot_rhos[b] = spearmanr(pred_vals[idx], dms_vals[idx]).statistic
    ci_lo, ci_hi = np.percentile(boot_rhos, [2.5, 97.5])
    print(f"Bootstrap 95% CI ({N_BOOTSTRAP} resamples): [{ci_lo:.4f}, {ci_hi:.4f}]")

    # --- Recompute excluding imputed variants, to see how much they're propping up rho ---
    non_imputed = joined[~joined["imputed"]]
    rho_no_impute, pval_no_impute = spearmanr(non_imputed["predicted_score"], non_imputed["DMS_score"])
    print(f"\nSpearman rho excluding imputed variants = {rho_no_impute:.4f} "
          f"(p={pval_no_impute:.2e}), n={len(non_imputed)}")
    print(f"({len(joined) - len(non_imputed)} imputed variants excluded)")

    print("\nNote: rho is a property of the whole variant set evaluated here, not of any "
          "individual mutation -- a single variant's predicted vs. DMS score can disagree "
          "sharply while rho over the full set is still strongly positive, or vice versa.")

    # ---------------- Scatter plot ----------------
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(
        non_imputed["predicted_score"], non_imputed["DMS_score"],
        s=8, alpha=0.4, color="#1f77b4", label=f"scored (n={len(non_imputed)})",
    )
    imputed_rows = joined[joined["imputed"]]
    ax.scatter(
        imputed_rows["predicted_score"], imputed_rows["DMS_score"],
        s=8, alpha=0.4, color="#d62728", label=f"imputed (n={len(imputed_rows)})",
    )
    ax.set_xlabel("Predicted score: log f(mut, pos) - log f(wt, pos)")
    ax.set_ylabel("DMS score (Starr 2020 ACE2 binding)")
    ax.set_title(f"PSSM vs DMS -- Spearman rho = {rho:.3f} (95% CI [{ci_lo:.3f}, {ci_hi:.3f}])")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_SCATTER, dpi=150)
    print(f"\nWrote {OUT_SCATTER}")

    out_meta = {
        "n_dms_variants": int(len(dms)),
        "n_prediction_rows": int(len(predictions)),
        "n_joined": int(len(joined)),
        "n_dropped_from_dms": int(n_dropped_from_dms),
        "n_dropped_from_predictions": int(n_dropped_from_predictions),
        "spearman_rho": float(rho),
        "spearman_pvalue": float(pval),
        "bootstrap_n": N_BOOTSTRAP,
        "bootstrap_ci_95_lo": float(ci_lo),
        "bootstrap_ci_95_hi": float(ci_hi),
        "spearman_rho_excluding_imputed": float(rho_no_impute),
        "n_excluding_imputed": int(len(non_imputed)),
    }
    with open(OUT_META, "w") as f:
        json.dump(out_meta, f, indent=2)
    print(f"Wrote {OUT_META}")


if __name__ == "__main__":
    main()
