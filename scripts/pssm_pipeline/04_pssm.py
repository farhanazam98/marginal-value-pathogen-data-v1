#!/usr/bin/env python3
"""Step 4: fit a site-independent model (PSSM) from the reweighted alignment.
For each column, compute weighted amino acid frequencies with a pseudocount.
Writes data/pssm_pipeline/pssm.npy (L x 20, rows sum to 1).

This approximates the site-wise maximum-entropy model of Hopf et al. 2017:
their per-column fields h_i(a) come out of a regularized pseudo-likelihood
fit of the full pairwise (coupling) model, so they implicitly absorb some of
what the pairwise couplings would otherwise explain. Here we skip that fit
and just read off regularized weighted frequencies column by column. The gap
is not a length in the pseudocount denominator.
"""

import json

import numpy as np

MSA_MATRIX = "data/pssm_pipeline/msa_clean.npy"
MSA_META = "data/pssm_pipeline/msa_clean_meta.json"
WEIGHTS = "data/pssm_pipeline/weights.npy"
OUT_PSSM = "data/pssm_pipeline/pssm.npy"
OUT_META = "data/pssm_pipeline/pssm_meta.json"

ALPHABET = "ACDEFGHIKLMNPQRSTVWY-"  # matches Step 2's encoding
GAP_CODE = len(ALPHABET) - 1
N_AA = 20

PSEUDOCOUNT = 1.0  # one pseudo-observation of each amino acid per column


def main():
    matrix = np.load(MSA_MATRIX)  # N x L, int8 codes 0..20 (20 = gap)
    weights = np.load(WEIGHTS)  # N
    meta = json.load(open(MSA_META))
    N, L = matrix.shape
    print(f"Loaded MSA matrix: N={N}, L={L}; weights: {weights.shape}")

    # --- Weighted counts per column, one amino acid at a time (20 passes over N x L) ---
    counts = np.zeros((N_AA, L))
    for a in range(N_AA):
        is_a = matrix == a  # N x L boolean
        counts[a] = weights @ is_a  # length-L weighted count of amino acid a per column

    non_gap = matrix != GAP_CODE
    neff_col = weights @ non_gap  # length L: weighted, non-gap effective count per column

    # --- Regularized frequencies: add-one-per-amino-acid pseudocount ---
    pssm = (counts + PSEUDOCOUNT) / (neff_col + N_AA * PSEUDOCOUNT)  # 20 x L
    pssm = pssm.T  # L x 20

    np.save(OUT_PSSM, pssm)

    print(f"\nWrote {OUT_PSSM}, shape {pssm.shape}")

    # ---------------- Sanity checks ----------------
    print("\n--- Sanity checks ---")

    row_sums = pssm.sum(axis=1)
    print(f"Row sums: min={row_sums.min():.6f}, max={row_sums.max():.6f} "
          f"(should all be 1.0): {np.allclose(row_sums, 1.0)}")

    query_row_idx = meta["query_row_index_in_final_matrix"]
    query_codes = matrix[query_row_idx]  # length L, should be all 0..19 (no gaps)
    n_query_gaps = int(np.sum(query_codes == GAP_CODE))
    print(f"Gaps in query row (should be 0): {n_query_gaps}")

    argmax_aa = np.argmax(pssm, axis=1)
    is_wt_top = argmax_aa == query_codes
    frac_wt_not_top = 1.0 - is_wt_top.mean()
    print(f"Fraction of columns where WT residue is NOT the highest-frequency residue: "
          f"{frac_wt_not_top:.3f} ({np.sum(~is_wt_top)} / {L})")

    # --- Per-position conservation (Shannon entropy in bits; low = conserved) ---
    entropy = -np.sum(pssm * np.log2(pssm), axis=1)  # length L
    print(f"\nEntropy (bits): min={entropy.min():.3f}, max={entropy.max():.3f}, "
          f"mean={entropy.mean():.3f} (max possible = log2(20) = {np.log2(20):.3f})")
    for p in (0, 10, 25, 50, 75, 90, 100):
        print(f"  p{p:>3}: {np.percentile(entropy, p):.3f}")

    kept_cols_1indexed = meta["kept_query_columns_1indexed"]  # length L, maps to full spike numbering
    order = np.argsort(entropy)
    print("\n10 most conserved columns (lowest entropy) -- (spike position, WT aa, entropy, top aa, top freq):")
    for j in order[:10]:
        wt_aa = ALPHABET[query_codes[j]]
        top_aa = ALPHABET[argmax_aa[j]]
        print(f"  pos {kept_cols_1indexed[j]:>4}  WT={wt_aa}  entropy={entropy[j]:.3f}  "
              f"top={top_aa} (freq={pssm[j, argmax_aa[j]]:.3f})")

    print("\n10 least conserved columns (highest entropy) -- (spike position, WT aa, entropy, top aa, top freq):")
    for j in order[-10:][::-1]:
        wt_aa = ALPHABET[query_codes[j]]
        top_aa = ALPHABET[argmax_aa[j]]
        print(f"  pos {kept_cols_1indexed[j]:>4}  WT={wt_aa}  entropy={entropy[j]:.3f}  "
              f"top={top_aa} (freq={pssm[j, argmax_aa[j]]:.3f})")

    out_meta = {
        "alphabet_aa_only": ALPHABET[:-1],
        "pseudocount": PSEUDOCOUNT,
        "L": L,
        "row_sums_ok": bool(np.allclose(row_sums, 1.0)),
        "n_query_gaps": n_query_gaps,
        "frac_columns_wt_not_top": float(frac_wt_not_top),
        "entropy_bits_min": float(entropy.min()),
        "entropy_bits_max": float(entropy.max()),
        "entropy_bits_mean": float(entropy.mean()),
    }
    with open(OUT_META, "w") as f:
        json.dump(out_meta, f, indent=2)
    print(f"\nWrote {OUT_META}")


if __name__ == "__main__":
    main()
