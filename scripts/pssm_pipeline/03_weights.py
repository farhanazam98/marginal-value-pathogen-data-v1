#!/usr/bin/env python3
"""Step 3: compute per-sequence weights via 99%-identity clustering, then Neff
and depth. Writes data/pssm_pipeline/weights.npy (N-length float64, aligned
to the rows of msa_clean.npy) and data/pssm_pipeline/weights_meta.json.

Weighting scheme (Hopf et al. 2017 / EVcouplings-style): for sequence i, S_i
is the number of sequences (including itself) whose pairwise identity to i is
>= 1 - theta. Weight pi_i = 1 / S_i, so a cluster of k near-duplicate
sequences collectively contributes 1 total unit of evidence instead of k.
Identity is the fraction of the L aligned columns where two sequences have
the same character (gap counts as a character like any other -- two
sequences gapped at the same position are "agreeing" there, which is the
standard convention for this kind of reweighting).
"""

import json

import numpy as np
from scipy.spatial.distance import pdist, squareform

MSA_MATRIX = "data/pssm_pipeline/msa_clean.npy"
MSA_META = "data/pssm_pipeline/msa_clean_meta.json"
OUT_WEIGHTS = "data/pssm_pipeline/weights.npy"
OUT_META = "data/pssm_pipeline/weights_meta.json"

THETA = 0.01  # cluster at 1 - theta = 99% identity
RELIABILITY_ID_CUTOFF = 0.90  # Methods A.6.1 reliability metric: Neff @ 90% ID
RELIABILITY_NEFF_THRESHOLD = 30
DEPTH_FLOOR = 1.0


def main():
    matrix = np.load(MSA_MATRIX)
    meta = json.load(open(MSA_META))
    N, L = matrix.shape
    print(f"Loaded MSA matrix: N={N} sequences, L={L} columns")

    print("\nComputing all-vs-all pairwise identity (this is the slow step)...")
    hamming_dist = pdist(matrix, metric="hamming")  # fraction of columns that differ, per pair
    identity = 1.0 - squareform(hamming_dist)
    print(f"Identity matrix shape: {identity.shape}")

    # --- theta = 0.01 -> cluster at 99% identity ---
    cutoff = 1.0 - THETA
    in_cluster = identity >= cutoff  # N x N boolean; diagonal is always True (self)
    cluster_size = in_cluster.sum(axis=1)  # S_i for each sequence
    weights = 1.0 / cluster_size
    Neff = weights.sum()
    depth = Neff / L

    np.save(OUT_WEIGHTS, weights)

    print(f"\n--- Weights at theta={THETA} (99% identity clustering) ---")
    print(f"N = {N}")
    print(f"Neff = {Neff:.2f}")
    print(f"Neff / L = {depth:.3f}")
    if depth < DEPTH_FLOOR:
        print(f"  *** WARNING: Neff/L = {depth:.3f} is below EVEREST's selection floor of {DEPTH_FLOOR}. ***")
        print("  *** EVEREST would not select this alignment for downstream modeling. ***")
    else:
        print(f"  Neff/L clears EVEREST's selection floor of {DEPTH_FLOOR}.")

    print("\nCluster size (S_i) distribution:")
    for p in (0, 10, 25, 50, 75, 90, 99, 100):
        print(f"  p{p:>3}: {np.percentile(cluster_size, p):.1f}")
    n_singletons = int(np.sum(cluster_size == 1))
    print(f"Sequences that are singletons (S_i=1, unique at 99% identity): {n_singletons} / {N}")
    largest_clusters = sorted(set(cluster_size.tolist()), reverse=True)[:5]
    print(f"5 largest cluster sizes observed: {largest_clusters}")

    # --- Reliability metric (Methods A.6.1): Neff @ 90% identity ---
    cutoff_90 = RELIABILITY_ID_CUTOFF
    in_cluster_90 = identity >= cutoff_90
    cluster_size_90 = in_cluster_90.sum(axis=1)
    weights_90 = 1.0 / cluster_size_90
    Neff_90 = weights_90.sum()

    print(f"\n--- Reliability metric: Neff @ {RELIABILITY_ID_CUTOFF:.0%} identity ---")
    print(f"Neff @ 90% ID = {Neff_90:.2f}")
    if Neff_90 >= RELIABILITY_NEFF_THRESHOLD:
        print(f"  Clears the paper's reliability threshold of {RELIABILITY_NEFF_THRESHOLD}.")
    else:
        print(f"  *** Does NOT clear the paper's reliability threshold of {RELIABILITY_NEFF_THRESHOLD}. ***")

    print("\nSanity checks:")
    print(f"  0 < Neff <= N: {0 < Neff <= N}")

    out_meta = {
        "theta": THETA,
        "identity_cutoff": cutoff,
        "N": N,
        "L": L,
        "Neff": Neff,
        "depth_Neff_over_L": depth,
        "depth_floor": DEPTH_FLOOR,
        "clears_depth_floor": bool(depth >= DEPTH_FLOOR),
        "reliability_id_cutoff": RELIABILITY_ID_CUTOFF,
        "Neff_at_90pct_identity": Neff_90,
        "reliability_neff_threshold": RELIABILITY_NEFF_THRESHOLD,
        "clears_reliability_threshold": bool(Neff_90 >= RELIABILITY_NEFF_THRESHOLD),
        "n_singleton_sequences": n_singletons,
    }
    with open(OUT_META, "w") as f:
        json.dump(out_meta, f, indent=2)
    print(f"\nWrote {OUT_WEIGHTS} and {OUT_META}")


if __name__ == "__main__":
    main()
