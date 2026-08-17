#!/usr/bin/env python3
"""Step 6: join each assay's predictions back to its DMS and compute Spearman rho.

Loops over the assays named by the active PROTEIN_CONFIG, reading the
predictions_<assay_id>.csv that step 05 wrote and producing scatter_<assay_id>.png
and evaluate_meta_<assay_id>.json for each.

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
from config import load_config

CKPT = "data/pssm_pipeline"

MUTANT_RE = re.compile(r"^([A-Z])(\d+)([A-Z])$")
N_BOOTSTRAP = 10_000
SEED = 0


def evaluate_assay(assay):
    dms_file = assay["csv"]
    predictions_file = f"{CKPT}/predictions_{assay['id']}.csv"
    out_scatter = f"{CKPT}/scatter_{assay['id']}.png"
    out_meta_path = f"{CKPT}/evaluate_meta_{assay['id']}.json"

    print(f"\n{'='*70}\nAssay '{assay['id']}' ({assay['label']})")
    dms = pd.read_csv(dms_file)[["mutant", "DMS_score"]]
    parsed = dms["mutant"].str.extract(MUTANT_RE)
    parsed.columns = ["wt_aa", "position", "mut_aa"]
    parsed["position"] = parsed["position"].astype(int)
    dms = pd.concat([dms, parsed], axis=1)
    print(f"DMS variants: {len(dms)}")

    predictions = pd.read_csv(predictions_file)
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
    ax.set_ylabel(f"DMS score ({assay['label']})")
    ax.set_title(f"PSSM vs {assay['id']} -- Spearman rho = {rho:.3f} (95% CI [{ci_lo:.3f}, {ci_hi:.3f}])")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_scatter, dpi=150)
    print(f"\nWrote {out_scatter}")

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
    with open(out_meta_path, "w") as f:
        json.dump(out_meta, f, indent=2)
    print(f"Wrote {out_meta_path}")


def main():
    cfg = load_config()
    for assay in cfg["assays"]:
        evaluate_assay(assay)


if __name__ == "__main__":
    main()
