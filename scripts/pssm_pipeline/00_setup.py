#!/usr/bin/env python3
"""Step 0: canonicalize the query sequence into data/query.fasta and sanity-check it.

Source is the `query_fasta` named by the active PROTEIN_CONFIG (default Spike:
data/protein.fasta, SARS-CoV-2 Spike glycoprotein, precursor/full-length numbering
starting at the initiator Met). For Spike that's the same numbering used by the
'mutant' column in the Starr 2020 DMS files, so no coordinate offset is needed
downstream; step 05 reconciles this per assay.
"""

from collections import Counter

from Bio import SeqIO
from config import load_config

CONFIG = load_config()
SOURCE_FASTA = CONFIG["query_fasta"]
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


if __name__ == "__main__":
    main()
