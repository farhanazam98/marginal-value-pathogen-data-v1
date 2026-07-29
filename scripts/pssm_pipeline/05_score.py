#!/usr/bin/env python3
"""Step 5: score every DMS variant with the Step 4 PSSM.
score(wt->mut, pos) = log f(mut, pos) - log f(wt, pos).
Writes data/pssm_pipeline/predictions.csv.

Numbering: the DMS's `mutant` column (e.g. "N331C") and query.fasta both use
full-spike, 1-indexed coordinates -- verified below before any scoring
happens. Positions with no surviving MSA column (dropped by Step 2's gap
filter) are imputed with the mean predicted score across the rest of the
protein (Methods A.3.5): this keeps the evaluated variant set identical
across alignment-based models and contributes no rank signal instead of a
biased one.
"""

import json
import re

import numpy as np
import pandas as pd

DMS_FILE = "data/SARS2_RBD_Starr_binding_dms.csv"
QUERY_FASTA = "data/pssm_pipeline/query.fasta"
MSA_META = "data/pssm_pipeline/msa_clean_meta.json"
PSSM = "data/pssm_pipeline/pssm.npy"
OUT_PREDICTIONS = "data/pssm_pipeline/predictions.csv"
OUT_META = "data/pssm_pipeline/predictions_meta.json"

ALPHABET = "ACDEFGHIKLMNPQRSTVWY-"
AA_TO_CODE = {c: i for i, c in enumerate(ALPHABET[:-1])}
MUTANT_RE = re.compile(r"^([A-Z])(\d+)([A-Z])$")

MATCH_RATE_FLOOR = 0.99  # stop if WT-residue reconciliation isn't ~100%


def main():
    dms = pd.read_csv(DMS_FILE)
    print(f"Loaded {len(dms)} DMS variants from {DMS_FILE}")

    parsed = dms["mutant"].str.extract(MUTANT_RE)
    parsed.columns = ["wt_aa", "position", "mut_aa"]
    if parsed.isna().any(axis=1).any():
        raise AssertionError("Some `mutant` entries didn't parse as WT+position+MUT -- stopping.")
    parsed["position"] = parsed["position"].astype(int)
    dms = pd.concat([dms[["mutant", "DMS_score"]], parsed], axis=1)

    query = open(QUERY_FASTA).read().split("\n", 1)[1].replace("\n", "")
    print(f"Query length: {len(query)}")

    # --- Numbering reconciliation: DMS WT residue must match query.fasta at that position ---
    query_wt_at_pos = dms["position"].apply(lambda p: query[p - 1])
    match = query_wt_at_pos == dms["wt_aa"]
    match_rate = match.mean()
    print(f"\nWT-residue reconciliation: {match.sum()}/{len(match)} match ({match_rate:.4f})")
    if match_rate < MATCH_RATE_FLOOR:
        mismatches = dms.loc[~match, ["mutant", "position", "wt_aa"]].head(10)
        print("Mismatch sample:")
        print(mismatches.to_string(index=False))
        raise AssertionError(
            f"WT-residue match rate {match_rate:.4f} is below the {MATCH_RATE_FLOOR} floor -- "
            "likely a coordinate offset between the DMS and query.fasta. Stopping rather than "
            "scoring against the wrong columns."
        )

    meta = json.load(open(MSA_META))
    kept_cols_1indexed = meta["kept_query_columns_1indexed"]
    pos_to_col = {pos: idx for idx, pos in enumerate(kept_cols_1indexed)}

    pssm = np.load(PSSM)
    log_pssm = np.log(pssm)
    L, n_aa = pssm.shape
    print(f"\nLoaded PSSM: {pssm.shape}")

    # --- Score every variant with a surviving MSA column; leave the rest NaN for now ---
    scores = np.full(len(dms), np.nan)
    imputed = np.zeros(len(dms), dtype=bool)
    for i, row in dms.iterrows():
        col = pos_to_col.get(row["position"])
        if col is None:
            imputed[i] = True
            continue
        wt_code = AA_TO_CODE[row["wt_aa"]]
        mut_code = AA_TO_CODE[row["mut_aa"]]
        scores[i] = log_pssm[col, mut_code] - log_pssm[col, wt_code]

    n_imputed = int(imputed.sum())
    mean_score = np.nanmean(scores)
    scores[imputed] = mean_score
    print(f"\nVariants scored directly: {len(dms) - n_imputed}")
    print(f"Variants imputed (position outside surviving MSA columns): {n_imputed}")
    print(f"Imputed value (mean predicted score across scored variants): {mean_score:.4f}")

    dms["predicted_score"] = scores
    dms["imputed"] = imputed
    dms[["mutant", "position", "wt_aa", "mut_aa", "DMS_score", "predicted_score", "imputed"]].to_csv(
        OUT_PREDICTIONS, index=False
    )
    print(f"\nWrote {OUT_PREDICTIONS}")

    # ---------------- Sanity checks ----------------
    print("\n--- Sanity checks ---")

    # WT->WT scores are exactly 0, for every position with a surviving column.
    covered_positions = sorted(set(dms.loc[~imputed, "position"]))
    wt_wt_scores = []
    for pos in covered_positions:
        col = pos_to_col[pos]
        wt_code = AA_TO_CODE[query[pos - 1]]
        wt_wt_scores.append(log_pssm[col, wt_code] - log_pssm[col, wt_code])
    wt_wt_scores = np.array(wt_wt_scores)
    print(f"WT->WT scores across {len(covered_positions)} covered positions: "
          f"all exactly 0.0: {bool(np.all(wt_wt_scores == 0.0))}")

    scored = dms.loc[~imputed, "predicted_score"]
    print(f"\nScore distribution (non-imputed, n={len(scored)}):")
    print(scored.describe())
    print(f"Fraction of non-imputed scores < 0 (deleterious-leaning): {(scored < 0).mean():.3f}")

    print(f"\n10 most deleterious predictions:")
    most_del = dms.nsmallest(10, "predicted_score")[["mutant", "position", "predicted_score", "imputed"]]
    print(most_del.to_string(index=False))

    print(f"\n10 least deleterious (most tolerated) predictions:")
    least_del = dms.nlargest(10, "predicted_score")[["mutant", "position", "predicted_score", "imputed"]]
    print(least_del.to_string(index=False))

    # Cross-reference with per-column entropy to see if extremes sit at conserved positions.
    freqs = pssm
    entropy = -np.sum(freqs * np.log2(freqs), axis=1)
    pos_to_entropy = {pos: entropy[col] for pos, col in pos_to_col.items()}
    print(f"\nEntropy (bits) at the 10 most-deleterious positions "
          f"(low = conserved; alignment-wide mean = {entropy.mean():.3f}):")
    for _, r in most_del.iterrows():
        e = pos_to_entropy.get(r["position"])
        print(f"  {r['mutant']:>6}  entropy={'imputed/no column' if e is None else f'{e:.3f}'}")

    print(f"\nEntropy (bits) at the 10 least-deleterious positions "
          f"(high = variable; alignment-wide mean = {entropy.mean():.3f}):")
    for _, r in least_del.iterrows():
        e = pos_to_entropy.get(r["position"])
        print(f"  {r['mutant']:>6}  entropy={'imputed/no column' if e is None else f'{e:.3f}'}")

    out_meta = {
        "n_variants": int(len(dms)),
        "n_scored_directly": int(len(dms) - n_imputed),
        "n_imputed": n_imputed,
        "imputed_value": float(mean_score),
        "wt_wt_all_zero": bool(np.all(wt_wt_scores == 0.0)),
        "frac_negative_scores_non_imputed": float((scored < 0).mean()),
        "predicted_score_mean": float(scored.mean()),
        "predicted_score_std": float(scored.std()),
    }
    with open(OUT_META, "w") as f:
        json.dump(out_meta, f, indent=2)
    print(f"\nWrote {OUT_META}")


if __name__ == "__main__":
    main()
