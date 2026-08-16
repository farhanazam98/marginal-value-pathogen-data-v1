#!/usr/bin/env python3
"""Step 5: score every DMS variant with the Step 4 PSSM.
score(wt->mut, pos) = log f(mut, pos) - log f(wt, pos).

Scores each assay named by the active PROTEIN_CONFIG against the single PSSM
built for this protein, writing predictions_<assay_id>.csv (+ a
predictions_meta_<assay_id>.json) per assay. Building the PSSM (steps 00-04) is
the expensive per-protein work; scoring is a cheap fan-out over the protein's
assays.

Numbering: each DMS's `mutant` column (e.g. "N331C") and query.fasta must share
1-indexed coordinates -- verified per assay below before any scoring happens.
Positions with no surviving MSA column (dropped by Step 2's gap filter) are
imputed with the mean predicted score across the rest of the protein (Methods
A.3.5): this keeps the evaluated variant set identical across alignment-based
models and contributes no rank signal instead of a biased one.
"""

import json
import re

import numpy as np
import pandas as pd
from config import load_config

QUERY_FASTA = "data/pssm_pipeline/query.fasta"
MSA_META = "data/pssm_pipeline/msa_clean_meta.json"
PSSM = "data/pssm_pipeline/pssm.npy"
CKPT = "data/pssm_pipeline"

ALPHABET = "ACDEFGHIKLMNPQRSTVWY-"
AA_TO_CODE = {c: i for i, c in enumerate(ALPHABET[:-1])}
MUTANT_RE = re.compile(r"^([A-Z])(\d+)([A-Z])$")

MATCH_RATE_FLOOR = 0.99  # stop if WT-residue reconciliation isn't ~100%


def score_assay(assay, query, pssm, log_pssm, pos_to_col, entropy):
    dms_file = assay["csv"]
    out_predictions = f"{CKPT}/predictions_{assay['id']}.csv"
    out_meta_path = f"{CKPT}/predictions_meta_{assay['id']}.json"

    dms = pd.read_csv(dms_file)
    print(f"\n{'='*70}\nAssay '{assay['id']}': loaded {len(dms)} DMS variants from {dms_file}")

    parsed = dms["mutant"].str.extract(MUTANT_RE)
    parsed.columns = ["wt_aa", "position", "mut_aa"]
    if parsed.isna().any(axis=1).any():
        raise AssertionError("Some `mutant` entries didn't parse as WT+position+MUT -- stopping.")
    parsed["position"] = parsed["position"].astype(int)
    dms = pd.concat([dms[["mutant", "DMS_score"]], parsed], axis=1)

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
        out_predictions, index=False
    )
    print(f"\nWrote {out_predictions}")

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
    with open(out_meta_path, "w") as f:
        json.dump(out_meta, f, indent=2)
    print(f"\nWrote {out_meta_path}")


def main():
    cfg = load_config()

    # Shared across all assays: build once, score many.
    query = open(QUERY_FASTA).read().split("\n", 1)[1].replace("\n", "")
    print(f"Query length: {len(query)}")

    meta = json.load(open(MSA_META))
    kept_cols_1indexed = meta["kept_query_columns_1indexed"]
    pos_to_col = {pos: idx for idx, pos in enumerate(kept_cols_1indexed)}

    pssm = np.load(PSSM)
    log_pssm = np.log(pssm)
    entropy = -np.sum(pssm * np.log2(pssm), axis=1)
    print(f"Loaded PSSM: {pssm.shape}")
    print(f"Scoring {len(cfg['assays'])} assay(s): {', '.join(a['id'] for a in cfg['assays'])}")

    for assay in cfg["assays"]:
        score_assay(assay, query, pssm, log_pssm, pos_to_col, entropy)


if __name__ == "__main__":
    main()
