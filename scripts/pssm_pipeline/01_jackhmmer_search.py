#!/usr/bin/env python3
"""Step 1: run jackhmmer locally to iteratively search the query against a
local UniRef100 FASTA snapshot, write the alignment of significant hits, and
print sanity checks.

This runs the *same* algorithm the EBI HMMER web service ran -- jackhmmer from
the HMMER suite -- but pointed at our own on-disk UniRef100 snapshot instead of
EBI's hosted databases. Doing it locally is both simpler (one subprocess call,
no submit/poll/download HTTP dance) and more faithful to the paper: the EBI API
has no UniRef100 option, so the web-service version fell back to `uniprot`;
here we search the actual UniRef100 release.

Threshold: EVEREST length-normalizes the inclusion cutoff as
bits/residue x L, passed here as a raw bit-score threshold (-T), not an E-value
threshold. A bit score is an absolute, database-size-independent measure of
alignment quality; an E-value is not -- it scales with database size, so the
*same* E-value cutoff would silently get stricter as the 2010->2026 snapshots
grow. A fixed bits/residue cutoff keeps "how good must a hit be" constant
across snapshots, which is the whole point of comparing across them.

Requires the HMMER suite on PATH: `conda activate marginal-value-pathogen-data`.
"""

import json
import os
import shutil
import subprocess
import time
from pathlib import Path

from Bio import SeqIO
from config import load_config

QUERY_FASTA = "data/pssm_pipeline/query.fasta"
# The local UniRef100 snapshot to search. Override with the SEQ_DB env var
# (e.g. `SEQ_DB=data/snapshots/uniref100_2015_01.fasta python ...`) instead
# of editing this default -- the sweep driver does exactly that per year.
SEQ_DB = os.environ.get("SEQ_DB", "data/uniref100_2010.fasta")
OUT_ALIGNMENT = "data/pssm_pipeline/msa_raw.sto"  # consumed by Step 2
OUT_TABLE = "data/pssm_pipeline/msa_raw_hits.tbl"  # per-sequence hit table
OUT_RUN_META = "data/pssm_pipeline/msa_raw_run_meta.json"  # run summary / debug

# Per-protein jackhmmer inclusion cutoff, bits/residue (lower = more permissive:
# deeper, noisier MSA). Set in the active PROTEIN_CONFIG, not here.
BITSCORE_PER_RESIDUE = load_config()["bitscore_per_residue"]
MAX_ITERATIONS = 5  # jackhmmer's default; it stops early once the hit set converges
CPU = 4  # per README benchmark, jackhmmer doesn't scale past ~2 cores on this workload


def parse_tblout(path):
    """Parse jackhmmer --tblout rows into the fields we report on.

    The format is whitespace-delimited with 18 fixed columns followed by a
    free-text description. For the full-sequence hit we want the target name
    (col 0), E-value (col 4) and bit score (col 5), plus the description.
    """
    hits = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split(None, 18)
            hits.append({
                "name": parts[0],
                "evalue": float(parts[4]),
                "score": float(parts[5]),
                "description": parts[18].strip() if len(parts) > 18 else "",
            })
    return hits


def main():
    record = next(SeqIO.parse(QUERY_FASTA, "fasta"))
    seq = str(record.seq)
    L = len(seq)
    T = BITSCORE_PER_RESIDUE * L

    if shutil.which("jackhmmer") is None:
        raise SystemExit(
            "jackhmmer not found on PATH. Activate the env first:\n"
            "  conda activate marginal-value-pathogen-data"
        )
    if not Path(SEQ_DB).exists():
        raise SystemExit(f"Sequence database not found: {SEQ_DB}")

    print(f"Query length L = {L}")
    print(f"BITSCORE_PER_RESIDUE = {BITSCORE_PER_RESIDUE}")
    print(f"Absolute inclusion threshold T = {BITSCORE_PER_RESIDUE} x {L} = {T} bits")
    print(f"Database = {SEQ_DB}")

    # jackhmmer <opts> <query.fasta> <seqdb.fasta>
    #   -T/--incT/--domT/--incdomT : the four bit-score thresholds the EBI
    #       payload set -- report-seq, include-seq, report-domain, include-domain.
    #   -A       : write the hit alignment (Stockholm) -- this is what Step 2 reads.
    #   --tblout : per-sequence hit table, parsed for the score/e-value prints.
    #   --noali  : keep stdout small (the alignment still goes to the -A file).
    cmd = [
        "jackhmmer",
        "-N", str(MAX_ITERATIONS),
        "--cpu", str(CPU),
        "-T", str(T), "--incT", str(T),
        "--domT", str(T), "--incdomT", str(T),
        "-A", OUT_ALIGNMENT,
        "--tblout", OUT_TABLE,
        "--noali",
        QUERY_FASTA, SEQ_DB,
    ]
    print("\nRunning:", " ".join(cmd))
    t0 = time.monotonic()
    result = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.monotonic() - t0
    if result.returncode != 0:
        raise SystemExit(f"jackhmmer failed (exit {result.returncode}):\n{result.stderr}")
    print(f"jackhmmer finished in {elapsed:.0f}s")

    # jackhmmer prints one "@@"-prefixed block per round summarizing what was
    # included and whether it converged; that's the useful signal in the log.
    # Round 1 (the initial phmmer pass) has no "@@ Round:" header -- only rounds
    # 2..K do -- so the number of rounds actually run is that count plus one.
    round_log = [ln.rstrip() for ln in result.stdout.splitlines() if ln.startswith("@@")]
    n_rounds = sum(1 for ln in round_log if ln.startswith("@@ Round:")) + 1
    converged = any("CONVERGED" in ln for ln in round_log)
    print(f"\nRounds run: {n_rounds}   Converged: {converged}")
    if not converged:
        print(f"  WARNING: did not converge within -N {MAX_ITERATIONS}; the hit set was still "
              "growing when it stopped. Consider a stricter threshold or more iterations.")

    hits = parse_tblout(OUT_TABLE)
    hits.sort(key=lambda h: h["score"], reverse=True)
    print(f"Number of significant hits: {len(hits)}")

    print("\nTop 10 hit descriptions:")
    for h in hits[:10]:
        print(f"  score={h['score']:8.1f}  evalue={h['evalue']:.2e}  {h['name']}  {h['description']}")

    aligned_records = list(SeqIO.parse(OUT_ALIGNMENT, "stockholm"))

    with open(OUT_RUN_META, "w") as f:
        json.dump(
            {
                "command": cmd,
                "database": SEQ_DB,
                "bitscore_per_residue": BITSCORE_PER_RESIDUE,
                "query_length": L,
                "threshold_bits": T,
                "elapsed_seconds": round(elapsed, 1),
                "n_hits": len(hits),
                "n_alignment_rows": len(aligned_records),
                "rounds": n_rounds,
                "converged": converged,
                "round_log": round_log,
            },
            f,
            indent=2,
        )
    print(f"\nWrote alignment to {OUT_ALIGNMENT}")
    print(f"Wrote hit table to {OUT_TABLE}")
    print(f"Wrote run metadata to {OUT_RUN_META}")

    print("\nSanity checks:")
    print(f"  Sequences in downloaded alignment: {len(aligned_records)}")

    exact_self_hits = [r for r in aligned_records if str(r.seq).replace("-", "").replace(".", "").upper() == seq.upper()]
    print(f"  Rows in alignment with sequence identical to query: {len(exact_self_hits)}")
    if exact_self_hits:
        print(f"    e.g. {exact_self_hits[0].id}")
    else:
        print("    WARNING: query sequence itself was not found verbatim among the hits.")
        print("    (Expected when the query post-dates the snapshot, e.g. SARS-CoV-2 vs a 2010 DB.)")
        if hits:
            top = hits[0]
            print(f"    Top hit: score={top['score']:.1f} evalue={top['evalue']:.2e} name={top['name']}")


if __name__ == "__main__":
    main()
