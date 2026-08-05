#!/usr/bin/env python3
"""Step 2 (revised): clean every candidate alignment from the threshold sweep.

Applies EVEREST's Methods A.3.1 filters -- drop columns >50% gaps, drop
sequences covering <50% of query length -- independently to each candidate raw
alignment produced by Step 1. Writes a labeled msa_clean_t{X}.npy (N x L
integer-encoded) plus msa_clean_meta_t{X}.json per candidate. It does NOT write
the canonical msa_clean.npy: Step 3 selects a winner and copies that one to the
canonical name, so Steps 4-6 stay untouched.

Coordinate mapping: the raw Stockholm alignment's "#=GC RF" line marks each
column as a match column ('x') or an insert column ('.'). Every candidate
converged at jackhmmer iteration 1 (single-query profile), so each has exactly
L match columns -- one per query residue -- and match-column index i *is* query
position i, no offset arithmetic. Insert columns are discarded. If a candidate
ever fails to reconstruct the query from its match columns, that assumption
broke and the script stops loudly rather than guessing coordinates.
"""

import json

import numpy as np

QUERY_FASTA = "data/pssm_pipeline/query.fasta"
OUT_DIR = "data/pssm_pipeline"

# Same sweep as Step 1 (0.3 reused). Highest threshold first for a readable table.
ALL_THRESHOLDS = [0.5, 0.3, 0.1, 0.05, 0.03, 0.01]

ALPHABET = "ACDEFGHIKLMNPQRSTVWY-"  # 20 standard amino acids + gap; index 20 = gap
GAP_CODE = len(ALPHABET) - 1
STANDARD_AA = set(ALPHABET[:-1])

COLUMN_GAP_MAX = 0.5
SEQUENCE_COVERAGE_MIN = 0.5


def raw_path(threshold):
    return f"{OUT_DIR}/msa_raw_t{threshold}.sto"


def clean_matrix_path(threshold):
    return f"{OUT_DIR}/msa_clean_t{threshold}.npy"


def clean_meta_path(threshold):
    return f"{OUT_DIR}/msa_clean_meta_t{threshold}.json"


def parse_stockholm(path):
    """Return (rf_string, {name: concatenated raw row string}) preserving row order."""
    rf = ""
    seq_parts = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#=GC RF"):
                rf += line.split(None, 2)[2].strip()
            elif line.startswith("#") or not line.strip() or line.strip() == "//":
                continue
            else:
                parts = line.split(None, 1)
                if len(parts) != 2:
                    continue
                name, seq = parts
                seq_parts.setdefault(name, []).append(seq.strip())
    return rf, {name: "".join(chunks) for name, chunks in seq_parts.items()}


def clean_one(threshold, query_fasta_seq):
    """Clean one candidate alignment; write labeled matrix + meta; return stats."""
    rf, seq_parts = parse_stockholm(raw_path(threshold))
    match_cols = [i for i, c in enumerate(rf) if c == "x"]
    L = len(match_cols)
    if L != len(query_fasta_seq):
        raise AssertionError(
            f"t{threshold}: {L} match columns but query is {len(query_fasta_seq)} residues -- "
            "the single-query-profile coordinate assumption broke. Stopping."
        )

    names = list(seq_parts.keys())
    N = len(names)

    raw_chars = np.empty((N, L), dtype="<U1")
    for i, name in enumerate(names):
        row = seq_parts[name]
        raw_chars[i] = [row[c] for c in match_cols]

    is_standard = np.isin(raw_chars, list(STANDARD_AA))
    is_gap = ~is_standard  # '-' and any non-standard call (e.g. 'X') both treated as gap
    n_nonstandard = int(np.sum((raw_chars != "-") & is_gap))

    # Locate the query row (its match-column string equals the full query).
    query_row_idx = None
    for i in range(N):
        if "".join(raw_chars[i]) == query_fasta_seq:
            query_row_idx = i
            break
    if query_row_idx is None:
        raise AssertionError(f"t{threshold}: no raw row is identical to query.fasta -- stopping.")

    # Column filter: drop columns >50% gaps over all N raw sequences.
    col_gap_frac = is_gap.mean(axis=0)
    keep_cols = np.where(col_gap_frac <= COLUMN_GAP_MAX)[0]

    # Sequence filter: drop sequences covering <50% of full query length L.
    seq_coverage = 1.0 - is_gap.mean(axis=1)
    keep_rows = np.where(seq_coverage >= SEQUENCE_COVERAGE_MIN)[0]
    if query_row_idx not in set(keep_rows.tolist()):
        raise AssertionError(f"t{threshold}: query row dropped by coverage filter -- stopping.")

    final_chars = raw_chars[np.ix_(keep_rows, keep_cols)]
    final_names = [names[i] for i in keep_rows]
    new_query_row_idx = keep_rows.tolist().index(query_row_idx)

    encoded = np.full(final_chars.shape, GAP_CODE, dtype=np.int8)
    for code, char in enumerate(ALPHABET[:-1]):
        encoded[final_chars == char] = code
    np.save(clean_matrix_path(threshold), encoded)

    # Sanity: query row gap-free, and matches query.fasta at every kept position.
    query_row_chars = final_chars[new_query_row_idx]
    n_gaps_in_query_row = int(np.sum(query_row_chars == "-"))
    expected_at_kept = "".join(query_fasta_seq[c] for c in keep_cols)
    query_matches_kept = ("".join(query_row_chars) == expected_at_kept)

    meta = {
        "threshold": threshold,
        "database": "uniprot",
        "alphabet": ALPHABET,
        "gap_code": GAP_CODE,
        "N_raw": N,
        "L_raw": L,
        "N_final": len(keep_rows),
        "L_final": len(keep_cols),
        "column_gap_max": COLUMN_GAP_MAX,
        "sequence_coverage_min": SEQUENCE_COVERAGE_MIN,
        "n_nonstandard_folded_to_gap": n_nonstandard,
        "kept_query_columns_1indexed": (np.array(keep_cols) + 1).tolist(),
        "query_row_name": names[query_row_idx],
        "query_row_index_in_final_matrix": new_query_row_idx,
        "sequence_names_in_final_matrix": final_names,
    }
    with open(clean_meta_path(threshold), "w") as f:
        json.dump(meta, f, indent=2)

    return {
        "threshold": threshold,
        "N_raw": N,
        "N_final": len(keep_rows),
        "L_raw": L,
        "L_final": len(keep_cols),
        "query_gapfree": n_gaps_in_query_row == 0,
        "query_matches_kept": query_matches_kept,
    }


def main():
    query_fasta_seq = open(QUERY_FASTA).read().split("\n", 1)[1].replace("\n", "")
    print(f"Query length L = {len(query_fasta_seq)}")
    print(f"Filters: drop columns >{COLUMN_GAP_MAX:.0%} gaps, drop sequences <{SEQUENCE_COVERAGE_MIN:.0%} query coverage\n")

    stats = []
    for threshold in ALL_THRESHOLDS:
        print(f"Cleaning t{threshold} ...")
        s = clean_one(threshold, query_fasta_seq)
        stats.append(s)
        print(f"  N {s['N_raw']} -> {s['N_final']}   L {s['L_raw']} -> {s['L_final']}   "
              f"query_gapfree={s['query_gapfree']}  query_matches_kept={s['query_matches_kept']}")

    # ---------------- Sanity table: shrinkage across thresholds ----------------
    print("\n=== Step 2 sanity: cleaning results per candidate ===")
    print(f"{'thresh':>7}  {'N_raw':>7}  {'N_final':>7}  {'seq_drop':>8}  "
          f"{'L_raw':>6}  {'L_final':>7}  {'col_drop':>8}  {'q_gapfree':>9}  {'q_match':>7}")
    all_ok = True
    for s in stats:
        seq_drop = s["N_raw"] - s["N_final"]
        col_drop = s["L_raw"] - s["L_final"]
        ok = s["query_gapfree"] and s["query_matches_kept"]
        all_ok = all_ok and ok
        print(f"{s['threshold']:>7}  {s['N_raw']:>7}  {s['N_final']:>7}  {seq_drop:>8}  "
              f"{s['L_raw']:>6}  {s['L_final']:>7}  {col_drop:>8}  "
              f"{str(s['query_gapfree']):>9}  {str(s['query_matches_kept']):>7}")

    print()
    if not all_ok:
        raise RuntimeError("A candidate failed the query-row sanity check (see q_gapfree/q_match above). Stopping.")
    print("All candidates cleaned; query row gap-free and consistent with query.fasta in every one.")
    print("Ready for Step 3 (weights, Neff/L, relevance, selection).")


if __name__ == "__main__":
    main()
