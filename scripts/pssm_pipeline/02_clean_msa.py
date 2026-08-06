#!/usr/bin/env python3
"""Step 2: map the raw alignment to query coordinates and apply EVEREST's
Methods A.3.1 filters (drop columns >50% gaps, drop sequences covering <50%
of query length). Writes data/pssm_pipeline/msa_clean.npy (N x L
integer-encoded) and data/pssm_pipeline/msa_clean_meta.json.

Coordinate mapping: the raw Stockholm alignment's "#=GC RF" line marks each
column as a match column ('x') or an insert column ('.'). This jackhmmer
profile was built from a single query sequence (iteration 1), so it has
exactly L match columns -- one per query residue -- and match-column index i
*is* query position i, with no further offset arithmetic required. Insert
columns hold residues that some hits have but the query doesn't, and are
simply discarded.
"""

import json

import numpy as np

RAW_ALIGNMENT = "data/pssm_pipeline/msa_raw.sto"
QUERY_FASTA = "data/pssm_pipeline/query.fasta"
OUT_MATRIX = "data/pssm_pipeline/msa_clean.npy"
OUT_META = "data/pssm_pipeline/msa_clean_meta.json"

ALPHABET = "ACDEFGHIKLMNPQRSTVWY-"  # 20 standard amino acids + gap; index 20 = gap
GAP_CODE = len(ALPHABET) - 1
STANDARD_AA = set(ALPHABET[:-1])

COLUMN_GAP_MAX = 0.5
SEQUENCE_COVERAGE_MIN = 0.5


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


def main():
    rf, seq_parts = parse_stockholm(RAW_ALIGNMENT)
    match_cols = [i for i, c in enumerate(rf) if c == "x"]
    L = len(match_cols)
    print(f"Alignment width: {len(rf)} columns ({L} match columns, {len(rf) - L} insert columns)")

    names = list(seq_parts.keys())
    N = len(names)
    print(f"N sequences in raw alignment: {N}")

    raw_chars = np.empty((N, L), dtype="<U1")
    for i, name in enumerate(names):
        row = seq_parts[name]
        raw_chars[i] = [row[c] for c in match_cols]

    is_standard = np.isin(raw_chars, list(STANDARD_AA))
    is_gap = ~is_standard  # '-' and any non-standard call (e.g. 'X') both treated as gap
    n_nonstandard = int(np.sum((raw_chars != "-") & is_gap))
    print(f"Non-gap, non-standard residue characters found (folded into gap channel): {n_nonstandard}")

    query_fasta_seq = open(QUERY_FASTA).read().split("\n", 1)[1].replace("\n", "")
    query_row_idx = None
    for i, name in enumerate(names):
        if "".join(raw_chars[i]) == query_fasta_seq:
            query_row_idx = i
            break
    if query_row_idx is None:
        raise AssertionError("No row in the raw alignment is identical to data/query.fasta -- stopping.")
    print(f"Query row located: {names[query_row_idx]} (raw row index {query_row_idx})")

    # --- Column filter: drop columns >50% gaps, over all N raw sequences ---
    col_gap_frac = is_gap.mean(axis=0)
    keep_cols = np.where(col_gap_frac <= COLUMN_GAP_MAX)[0]
    print(f"\nColumns before filter: {L}")
    print(f"Columns with >{COLUMN_GAP_MAX:.0%} gaps: {L - len(keep_cols)}")
    print(f"Columns after filter: {len(keep_cols)}")

    # --- Sequence filter: drop sequences covering <50% of query length ---
    # Coverage is measured against the full query length L, not the
    # post-column-filter column count -- query length is a fixed property of
    # the protein, independent of which columns this particular alignment
    # run happened to fill in.
    seq_coverage = 1.0 - is_gap.mean(axis=1)
    keep_rows = np.where(seq_coverage >= SEQUENCE_COVERAGE_MIN)[0]
    print(f"\nSequences before filter: {N}")
    print(f"Sequences covering <{SEQUENCE_COVERAGE_MIN:.0%} of query length: {N - len(keep_rows)}")
    print(f"Sequences after filter: {len(keep_rows)}")

    if query_row_idx not in set(keep_rows.tolist()):
        raise AssertionError("Query row was dropped by the sequence coverage filter -- something is wrong.")

    final_chars = raw_chars[np.ix_(keep_rows, keep_cols)]
    final_names = [names[i] for i in keep_rows]
    new_query_row_idx = keep_rows.tolist().index(query_row_idx)

    print(f"\nFinal matrix shape: {final_chars.shape}")

    encoded = np.full(final_chars.shape, GAP_CODE, dtype=np.int8)
    for code, char in enumerate(ALPHABET[:-1]):
        encoded[final_chars == char] = code

    np.save(OUT_MATRIX, encoded)

    meta = {
        "alphabet": ALPHABET,
        "gap_code": GAP_CODE,
        "N_raw": N,
        "L_raw": L,
        "N_final": len(keep_rows),
        "L_final": len(keep_cols),
        "column_gap_max": COLUMN_GAP_MAX,
        "sequence_coverage_min": SEQUENCE_COVERAGE_MIN,
        "kept_query_columns_1indexed": (np.array(keep_cols) + 1).tolist(),
        "query_row_name": names[query_row_idx],
        "query_row_index_in_final_matrix": new_query_row_idx,
        "sequence_names_in_final_matrix": final_names,
    }
    with open(OUT_META, "w") as f:
        json.dump(meta, f, indent=2)

    print(f"\nWrote {OUT_MATRIX} and {OUT_META}")

    print("\nSanity checks:")
    query_row_chars = final_chars[new_query_row_idx]
    query_row_str = "".join(query_row_chars)
    n_gaps_in_query_row = int(np.sum(query_row_chars == "-"))
    print(f"  Gaps in query row: {n_gaps_in_query_row} (should be 0)")

    if len(keep_cols) == L:
        matches_full_query = query_row_str == query_fasta_seq
        print(f"  Query row identical to data/query.fasta: {matches_full_query} (no columns dropped)")
    else:
        print(f"  Columns were dropped ({L - len(keep_cols)} of {L}), so the query row is a "
              f"{len(keep_cols)}-residue subsequence of the full {L}-residue query and cannot be "
              f"literally identical to data/query.fasta. Checking it matches at the kept positions instead:")
        expected = "".join(query_fasta_seq[c] for c in keep_cols)
        print(f"    Matches query.fasta at all kept positions: {query_row_str == expected}")

    print(f"\n  Column gap-fraction distribution (pre-filter, over all {L} columns):")
    for p in (0, 10, 25, 50, 75, 90, 100):
        print(f"    p{p:>3}: {np.percentile(col_gap_frac, p):.4f}")

    print(f"\n  Sequence coverage distribution (pre-filter, over all {N} sequences):")
    for p in (0, 10, 25, 50, 75, 90, 100):
        print(f"    p{p:>3}: {np.percentile(seq_coverage, p):.4f}")


if __name__ == "__main__":
    main()
