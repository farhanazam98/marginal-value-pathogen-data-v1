#!/usr/bin/env python3
"""Step 0: canonicalize the query sequence into data/query.fasta and sanity-check it.

Source: data/protein.fasta (SARS-CoV-2 Spike glycoprotein, precursor/full-length
numbering starting at the initiator Met). This is the same numbering used by the
'mutant' and 'mutated_sequence' columns in the Starr 2020 DMS file, so no coordinate
offset is needed downstream.
"""

from collections import Counter

from Bio import SeqIO

SOURCE_FASTA = "data/protein.fasta"
QUERY_FASTA = "data/pssm_pipeline/query.fasta"

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def main():
    record = next(SeqIO.parse(SOURCE_FASTA, "fasta"))
    seq = str(record.seq).upper()

    SeqIO.write(record, QUERY_FASTA, "fasta")

    L = len(seq)
    composition = Counter(seq)
    non_standard = sorted(set(seq) - STANDARD_AA)

    print(f"Wrote {QUERY_FASTA} from {SOURCE_FASTA}")
    print(f"Record ID: {record.id}")
    print(f"Length L = {L}")
    print("Residue composition (count, fraction):")
    for aa in sorted(composition):
        count = composition[aa]
        print(f"  {aa}: {count:4d}  ({count / L:.3%})")

    print(f"\nNon-standard characters found: {non_standard if non_standard else 'none'}")
    assert not non_standard, f"Non-standard amino acid characters present: {non_standard}"

    print(f"\nSequence starts: {seq[:20]}...")
    print(f"Sequence ends:   ...{seq[-20:]}")

    print("\nSanity checks:")
    print(f"  L == 1273 (full SARS-CoV-2 Spike length)? {'PASS' if L == 1273 else 'FAIL'} (L={L})")
    starts_met = seq.startswith("M")
    print(f"  Starts with initiator Met? {'PASS' if starts_met else 'FAIL'}")
    assert L == 1273, f"Expected full Spike length 1273, got {L}"
    assert starts_met, "Expected sequence to start with initiator Met"


if __name__ == "__main__":
    main()
