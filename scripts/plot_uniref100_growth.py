#!/usr/bin/env python3
"""Plot the growth of the combined UniRef50+90+100 archive across yearly UniProt releases."""

import csv

import matplotlib.pyplot as plt

CSV_PATH = "data/uniref100_yearly_release_sizes.csv"
OUTPUT_PATH = "data/uniref100_growth.png"

SUFFIX_TO_GB = {"K": 1 / 1024 / 1024, "M": 1 / 1024, "G": 1, "T": 1024}


def size_to_gb(raw):
    suffix = raw[-1]
    if suffix in SUFFIX_TO_GB:
        return float(raw[:-1]) * SUFFIX_TO_GB[suffix]
    return float(raw) / 1024 / 1024 / 1024


def main():
    with open(CSV_PATH, newline="") as f:
        rows = list(csv.DictReader(f))

    years = [int(row["year"]) for row in rows]
    sizes_gb = [size_to_gb(row["uniref_archive_size_raw"]) for row in rows]

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(years, sizes_gb, color="#4C72B0")
    ax.set_xlabel("Year (first release)")
    ax.set_ylabel("Archive size (GB)")
    ax.set_title(
        "UniProt UniRef archive growth, 2010–2026\n"
        "(combined UniRef50+90+100 download; not UniRef100 alone)"
    )
    ax.set_xticks(years)
    ax.set_xticklabels(years, rotation=45, ha="right")
    fig.tight_layout()
    fig.savefig(OUTPUT_PATH, dpi=150)
    print(f"Wrote {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
