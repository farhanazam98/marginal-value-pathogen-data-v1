#!/usr/bin/env python3
"""Diagnose the RBD gap-column-filter artifact.

02_clean_msa.py drops any alignment column that is >50% gaps across all
sequences (see CLAUDE.md's "Key methodology"). The hypothesis this script
checks: RBD columns are getting dropped because distantly-related homologs
align poorly across that specific region (it's the most variable part of the
protein) while still covering the rest of the query well enough to pass the
per-sequence coverage filter -- i.e. a handful of divergent sequences carry
most of the RBD gaps.

For every sequence in a jackhmmer alignment (.sto or .a3m), this computes:
  - sequence_identity_pct: % identity to the wild-type query, over all match
    columns (denominator = query length L, so a sequence with gaps scores
    lower here too -- this is "% of the query correctly covered", not just
    identity among aligned positions).
  - rbd_gap_count: number of '-' characters within the query's RBD window
    (positions --rbd-start..--rbd-end, 1-indexed, inclusive) among that
    sequence's match-column characters.
Plotting one against the other should show low-identity sequences clustering
at high RBD gap counts if the hypothesis holds.

Coordinate mapping follows the same convention as 02_clean_msa.py: for .sto,
match columns are the ones marked 'x' on the "#=GC RF" line, and the i-th
match column (0-indexed) is query position i -- this jackhmmer profile was
built from a single query, so there's no further offset. For .a3m, match
columns are the uppercase/'-' characters; lowercase characters are insert
states relative to the query and are ignored, per the same convention.

Memory: a Stockholm alignment interleaves each sequence's characters across
many blocks (see any .sto file -- blocks of ~70 columns, one per sequence,
repeated), so no single pass can finalize a sequence's row before EOF.
Rather than buffering full per-sequence strings the way 02_clean_msa.py does
(O(N x L) characters, fine for a few thousand rows but wasteful at the hit
counts a 59 GB snapshot can produce), this script keeps only two running
integers per sequence (identity matches, RBD gaps) and one block's worth of
lines at a time -- O(N x block_width), independent of total alignment
length. .a3m has no interleaving, so it's read with Bio.SeqIO.parse's
generator interface and each row is scored and written immediately.
"""

import argparse
import csv
import sys
from pathlib import Path

from Bio import SeqIO


def load_query(path):
    record = next(SeqIO.parse(path, "fasta"))
    return str(record.seq).upper()


def iter_stockholm_blocks(path):
    """Yield (rf_string, [(name, chars), ...]) once per interleaved block."""
    block_lines = []
    block_rf = None
    with open(path) as f:
        for line in f:
            if line.startswith("#=GC RF"):
                block_rf = line.split(None, 2)[2].strip()
                continue
            if line.startswith("#"):
                continue
            if not line.strip() or line.strip() == "//":
                if block_lines or block_rf:
                    yield block_rf, block_lines
                block_lines, block_rf = [], None
                continue
            parts = line.split(None, 1)
            if len(parts) != 2:
                continue
            name, chars = parts
            block_lines.append((name, chars.strip()))
    if block_lines or block_rf:
        yield block_rf, block_lines


def process_stockholm(path, query_seq, rbd_start, rbd_end):
    """Yield (sequence_id, sequence_identity_pct, rbd_gap_count) per sequence."""
    rbd_lo, rbd_hi = rbd_start - 1, rbd_end - 1  # 0-indexed, inclusive

    match_counter = 0  # running count of match columns seen so far == query position
    order = []
    acc = {}  # name -> [identity_matches, rbd_gaps]

    for block_rf, block_lines in iter_stockholm_blocks(path):
        if not block_rf:
            continue
        match_positions = []
        for col, c in enumerate(block_rf):
            if c == "x":
                match_positions.append((col, match_counter))
                match_counter += 1
        if not match_positions:
            continue

        for name, chars in block_lines:
            row = acc.get(name)
            if row is None:
                row = acc[name] = [0, 0]
                order.append(name)
            for col, qpos in match_positions:
                if col >= len(chars) or qpos >= len(query_seq):
                    continue
                c = chars[col]
                if c != "-" and c.upper() == query_seq[qpos]:
                    row[0] += 1
                if rbd_lo <= qpos <= rbd_hi and c == "-":
                    row[1] += 1

    L = match_counter
    if L != len(query_seq):
        print(
            f"WARNING: alignment has {L} match columns but the query is {len(query_seq)} aa; "
            "identity percentages are relative to the alignment's match-column count, not the "
            "query length.",
            file=sys.stderr,
        )

    for name in order:
        identity_matches, rbd_gaps = acc[name]
        yield name, 100.0 * identity_matches / L, rbd_gaps


def process_a3m(path, query_seq, rbd_start, rbd_end):
    """Yield (sequence_id, sequence_identity_pct, rbd_gap_count) per sequence."""
    rbd_lo, rbd_hi = rbd_start - 1, rbd_end - 1  # 0-indexed, inclusive
    L = len(query_seq)

    for record in SeqIO.parse(path, "fasta"):
        raw = str(record.seq)
        match_chars = [c for c in raw if not c.islower() and c != "."]
        if len(match_chars) != L:
            print(
                f"WARNING: {record.id} has {len(match_chars)} match-column characters, "
                f"expected {L} (query length); skipping.",
                file=sys.stderr,
            )
            continue
        identity_matches = sum(
            1 for qpos, c in enumerate(match_chars) if c != "-" and c.upper() == query_seq[qpos]
        )
        rbd_gaps = sum(1 for c in match_chars[rbd_lo:rbd_hi + 1] if c == "-")
        yield record.id, 100.0 * identity_matches / L, rbd_gaps


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--alignment", required=True, type=Path, help="Jackhmmer alignment file (.sto or .a3m)")
    parser.add_argument("--format", choices=["auto", "sto", "a3m"], default="auto",
                         help="Alignment format; auto infers from --alignment's extension (default: auto)")
    parser.add_argument("--query", required=True, type=Path, help="Wild-type query FASTA (single record)")
    parser.add_argument("--rbd-start", type=int, default=361,
                         help="RBD window start, 1-indexed query position, inclusive (default: 361)")
    parser.add_argument("--rbd-end", type=int, default=413,
                         help="RBD window end, 1-indexed query position, inclusive (default: 413)")
    parser.add_argument("--output", type=Path, default=Path("data/pssm_pipeline/rbd_gap_diagnostic.csv"),
                         help="Output CSV path (default: data/pssm_pipeline/rbd_gap_diagnostic.csv)")
    args = parser.parse_args()

    fmt = args.format
    if fmt == "auto":
        suffix = args.alignment.suffix.lower()
        if suffix == ".sto":
            fmt = "sto"
        elif suffix == ".a3m":
            fmt = "a3m"
        else:
            raise SystemExit(f"Can't infer format from extension {suffix!r}; pass --format sto|a3m explicitly.")

    query_seq = load_query(args.query)
    print(f"Query length L = {len(query_seq)}")
    print(f"RBD window: query positions {args.rbd_start}-{args.rbd_end} (1-indexed, inclusive)")
    print(f"Alignment: {args.alignment} (format={fmt})")

    rows = (
        process_stockholm(args.alignment, query_seq, args.rbd_start, args.rbd_end)
        if fmt == "sto"
        else process_a3m(args.alignment, query_seq, args.rbd_start, args.rbd_end)
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sequence_id", "sequence_identity_pct", "rbd_gap_count"])
        for seq_id, identity_pct, rbd_gaps in rows:
            writer.writerow([seq_id, f"{identity_pct:.2f}", rbd_gaps])
            n += 1

    print(f"Wrote {n} rows to {args.output}")


if __name__ == "__main__":
    main()
