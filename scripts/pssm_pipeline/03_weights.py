#!/usr/bin/env python3
"""Step 3 (revised): weight every candidate, then select one alignment.

For each cleaned candidate from Step 2 this computes per-sequence weights via
99%-identity clustering (Hopf et al. 2017 / EVcouplings-style), Neff, depth
Neff/L, and the alignment-relevance fraction, writing labeled weights_t{X}.npy
+ weights_meta_t{X}.json. Then it applies EVEREST's Methods A.4.1 selection:

  1. keep candidates with Neff/L > 1.0 (the depth floor),
  2. among survivors, choose the one maximizing alignment relevance,
  3. fallback: if none clear the floor, choose max raw Neff/L and flag it loudly.

The winner's trio (matrix, meta, weights) is copied onto the canonical filenames
so the unchanged Steps 4-6 fit/score/evaluate the selected alignment.

CAVEAT (this run, deliberately not worked around): EVEREST's selection is a
whole-protein criterion, but this DMS covers only a sub-region (RBD 331-531).
When no candidate clears the depth floor, the max-Neff/L fallback prefers the
loosest threshold -- whose divergent homologs align poorly through the variable
RBD, dropping every assayed column and yielding an undefined rho. The stricter
thresholds (0.3, 0.5) keep the RBD; among candidates, 0.5 scores best (rho 0.266
vs 0.3's 0.248). The selection is kept faithful to the paper rather than patched,
because the EBI web service (live 'uniprot', single-pass search, the hits/rows
discrepancy) makes different assumptions than the paper's local UniRef100
pipeline -- so this degeneracy may be an artifact of the substitute and should be
re-judged there, not compensated for here.

Two distinct 90%-identity quantities, do not conflate:
  - alignment relevance = fraction of sequences with identity to the query
    STRICTLY > 0.90 (the selection tie-breaker; a fraction in [0,1]).
  - Neff @ 90% ID = effective count when clustering at >= 0.90 (Methods A.6.1
    reliability metric, threshold ~30; an absolute count). Reported, not used
    for selection.
"""

import filecmp
import json
import shutil

import numpy as np
from scipy.spatial.distance import pdist, squareform

OUT_DIR = "data/pssm_pipeline"
ALL_THRESHOLDS = [0.5, 0.3, 0.1, 0.05, 0.03, 0.01]

THETA = 0.01               # cluster at 1 - theta = 99% identity
DEPTH_FLOOR = 1.0          # Neff/L selection floor
RELEVANCE_ID_CUTOFF = 0.90  # relevance: fraction of sequences with identity to query > this
RELIABILITY_ID_CUTOFF = 0.90  # Neff @ 90% ID reliability metric (>=)
RELIABILITY_NEFF_THRESHOLD = 30

CANON_CLEAN = f"{OUT_DIR}/msa_clean.npy"
CANON_CLEAN_META = f"{OUT_DIR}/msa_clean_meta.json"
CANON_WEIGHTS = f"{OUT_DIR}/weights.npy"
CANON_WEIGHTS_META = f"{OUT_DIR}/weights_meta.json"


def clean_matrix_path(t):
    return f"{OUT_DIR}/msa_clean_t{t}.npy"


def clean_meta_path(t):
    return f"{OUT_DIR}/msa_clean_meta_t{t}.json"


def weights_path(t):
    return f"{OUT_DIR}/weights_t{t}.npy"


def weights_meta_path(t):
    return f"{OUT_DIR}/weights_meta_t{t}.json"


def compute_candidate(threshold):
    """Load one cleaned candidate, compute weights + selection stats, write labeled
    weights outputs, and return a stats dict."""
    matrix = np.load(clean_matrix_path(threshold))
    meta = json.load(open(clean_meta_path(threshold)))
    N, L = matrix.shape
    query_row_idx = meta["query_row_index_in_final_matrix"]

    # All-vs-all identity: fraction of the L columns where two rows agree
    # (gap counts as a character, standard reweighting convention).
    identity = 1.0 - squareform(pdist(matrix, metric="hamming"))

    # Weights at 99% identity: pi_i = 1 / (#sequences within 99% identity of i).
    in_cluster = identity >= (1.0 - THETA)
    cluster_size = in_cluster.sum(axis=1)
    weights = 1.0 / cluster_size
    Neff = float(weights.sum())
    depth = Neff / L

    # Alignment relevance: fraction of sequences strictly >90% identical to query
    # (the query self-row is one such sequence -- a +1/N effect, noted).
    query_identities = identity[query_row_idx]
    relevance = float(np.mean(query_identities > RELEVANCE_ID_CUTOFF))

    # Reliability metric (distinct from relevance): Neff clustering at >=90% ID.
    in_cluster_90 = identity >= RELIABILITY_ID_CUTOFF
    Neff_90 = float((1.0 / in_cluster_90.sum(axis=1)).sum())

    np.save(weights_path(threshold), weights)
    wmeta = {
        "threshold": threshold,
        "database": "uniprot",
        "theta": THETA,
        "N": N,
        "L": L,
        "Neff": Neff,
        "depth_Neff_over_L": depth,
        "depth_floor": DEPTH_FLOOR,
        "clears_depth_floor": bool(depth > DEPTH_FLOOR),
        "relevance_id_cutoff": RELEVANCE_ID_CUTOFF,
        "alignment_relevance_frac_gt_90id": relevance,
        "reliability_id_cutoff": RELIABILITY_ID_CUTOFF,
        "Neff_at_90pct_identity": Neff_90,
        "clears_reliability_threshold": bool(Neff_90 >= RELIABILITY_NEFF_THRESHOLD),
    }
    with open(weights_meta_path(threshold), "w") as f:
        json.dump(wmeta, f, indent=2)

    return {
        "threshold": threshold, "N": N, "L": L, "Neff": Neff, "depth": depth,
        "clears_floor": depth > DEPTH_FLOOR, "relevance": relevance, "Neff_90": Neff_90,
    }


def select(stats):
    """Two-step EVEREST selection. Returns (winner_stat, fallback_triggered)."""
    survivors = [s for s in stats if s["clears_floor"]]
    if survivors:
        winner = max(survivors, key=lambda s: s["relevance"])
        return winner, False
    winner = max(stats, key=lambda s: s["depth"])
    return winner, True


def copy_trio(threshold):
    """Copy the winner's matrix, meta, and weights onto the canonical filenames."""
    pairs = [
        (clean_matrix_path(threshold), CANON_CLEAN),
        (clean_meta_path(threshold), CANON_CLEAN_META),
        (weights_path(threshold), CANON_WEIGHTS),
        (weights_meta_path(threshold), CANON_WEIGHTS_META),
    ]
    for src, dst in pairs:
        shutil.copyfile(src, dst)
    return pairs


def main():
    print(f"Weighting {len(ALL_THRESHOLDS)} candidates (99% identity clustering); "
          f"depth floor Neff/L > {DEPTH_FLOOR}\n")

    stats = []
    for threshold in ALL_THRESHOLDS:
        print(f"Weighting t{threshold} ...")
        s = compute_candidate(threshold)
        stats.append(s)
        print(f"  N={s['N']} L={s['L']} Neff={s['Neff']:.1f} Neff/L={s['depth']:.3f} "
              f"relevance={s['relevance']:.3f} Neff@90={s['Neff_90']:.1f}")

    winner, fallback = select(stats)

    # ---------------- Audit table ----------------
    print("\n=== Step 3 audit: every candidate, selection basis ===")
    print(f"{'thresh':>7}  {'N':>6}  {'L':>5}  {'Neff':>8}  {'Neff/L':>7}  "
          f"{'floor>1':>7}  {'relevance':>9}  {'Neff@90':>8}  {'selected':>8}")
    for s in stats:
        mark = "  <== SELECTED" if s["threshold"] == winner["threshold"] else ""
        print(f"{s['threshold']:>7}  {s['N']:>6}  {s['L']:>5}  {s['Neff']:>8.1f}  {s['depth']:>7.3f}  "
              f"{str(s['clears_floor']):>7}  {s['relevance']:>9.3f}  {s['Neff_90']:>8.1f}  "
              f"{('YES' if s['threshold'] == winner['threshold'] else ''):>8}{mark}")

    print("\n--- Selection ---")
    n_survivors = sum(1 for s in stats if s["clears_floor"])
    if fallback:
        print("*** FALLBACK TRIGGERED: no candidate cleared Neff/L > 1.0. ***")
        print("*** Selected the maximum raw Neff/L, faithful to EVEREST's stated fallback. ***")
        print("*** WARNING: for a sub-region DMS (RBD 331-531) this can pick a globally-deep ***")
        print("*** alignment that drops the assayed region entirely -> undefined rho. See the ***")
        print("*** module docstring: kept faithful rather than patched; the EBI web service ***")
        print("*** makes different assumptions than the paper's UniRef100 pipeline, so this ***")
        print("*** should be re-judged there, not worked around on the web-service substitute. ***")
        print(f"Selected t{winner['threshold']} by max Neff/L = {winner['depth']:.3f} "
              f"(still below the {DEPTH_FLOOR} floor).")
    else:
        print(f"{n_survivors} candidate(s) cleared Neff/L > {DEPTH_FLOOR}.")
        print(f"Selected t{winner['threshold']} by max relevance = {winner['relevance']:.3f} "
              f"among floor-clearing candidates (Neff/L = {winner['depth']:.3f}).")

    # ---------------- Copy winner trio to canonical names, verify byte-for-byte ----------------
    print("\n--- Copying winner to canonical filenames ---")
    pairs = copy_trio(winner["threshold"])
    all_match = True
    for src, dst in pairs:
        same = filecmp.cmp(src, dst, shallow=False)
        all_match = all_match and same
        print(f"  {dst}  <=  {src}   byte-identical={same}")

    print()
    if not all_match:
        raise RuntimeError("A canonical copy does not match its source byte-for-byte. Stopping.")
    print(f"Canonical trio now points at t{winner['threshold']}. Steps 4-6 will run on it unchanged.")
    print(f"Selected Neff/L = {winner['depth']:.3f}"
          + ("  (FALLBACK: below floor)" if fallback else "  (cleared floor)"))


if __name__ == "__main__":
    main()
