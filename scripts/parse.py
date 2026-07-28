#!/usr/bin/env python3
"""Parse the Starr 2020 SARS-CoV-2 RBD DMS CSV and print the N501Y row."""

import csv
import sys

CSV_PATH = "data/SARS2_RBD_Starr_binding_dms.csv"
TARGET_MUTANT = "N501Y"

COLUMN_DESCRIPTIONS = {
    "mutant": "Mutation identifier (wild-type residue, position, mutant residue)",
    "mutated_sequence": "Full RBD/spike protein sequence with the substitution applied",
    "DMS_score": "Deep mutational scanning ACE2 binding affinity score",
    "DMS_score_bin": "Binarized binding score (1 = binds/tolerated, 0 = does not bind)",
}


def main():
    with open(CSV_PATH, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["mutant"] == TARGET_MUTANT:
                print(f"Found mutation {TARGET_MUTANT}:\n")
                for column, value in row.items():
                    description = COLUMN_DESCRIPTIONS.get(column, "No description available")
                    if column == "mutated_sequence" and len(value) > 10:
                        value = f"{value[:10]}... ({len(value)} aa)"
                    print(f"{column}: {value}\n  -> {description}\n")
                return

    print(f"Mutation {TARGET_MUTANT} not found in {CSV_PATH}", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    main()
