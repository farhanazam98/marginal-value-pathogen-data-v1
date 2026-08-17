"""Load the per-protein pipeline config selected by the PROTEIN_CONFIG env var.

One small YAML file per protein names the query FASTA, the jackhmmer bit-score
threshold, and the list of DMS assays to score against it (see config/spike.yaml).
Steps 00/01 depend only on the query + threshold; steps 05/06 fan out over the
assay list. Kept deliberately minimal -- when this grows to many proteins it
becomes a single manifest CSV instead (see README).

Run as a script to print the input files a config references, one per line:
    python scripts/pssm_pipeline/config.py --input-files
The sweep driver uses that to symlink the right inputs into each sandbox.
"""

import os
import sys

import yaml

DEFAULT_CONFIG = "config/spike.yaml"
REQUIRED_KEYS = ("name", "query_fasta", "bitscore_per_residue", "assays")


def load_config(path=None):
    path = path or os.environ.get("PROTEIN_CONFIG", DEFAULT_CONFIG)
    try:
        with open(path) as fh:
            cfg = yaml.safe_load(fh)
    except FileNotFoundError:
        raise SystemExit(f"PROTEIN_CONFIG not found: {path}")

    missing = [k for k in REQUIRED_KEYS if k not in cfg]
    if missing:
        raise SystemExit(f"{path}: missing required key(s): {', '.join(missing)}")
    if not cfg["assays"]:
        raise SystemExit(f"{path}: 'assays' is empty -- at least one DMS assay is required")
    for a in cfg["assays"]:
        if "id" not in a or "csv" not in a:
            raise SystemExit(f"{path}: each assay needs an 'id' and a 'csv'")
        a.setdefault("label", a["id"])
    return cfg


def input_files(cfg):
    """Every on-disk input the config points at: the query plus each assay CSV.
    These are what the sweep sandbox must symlink in."""
    return [cfg["query_fasta"]] + [a["csv"] for a in cfg["assays"]]


if __name__ == "__main__":
    cfg = load_config()
    if "--input-files" in sys.argv:
        for f in input_files(cfg):
            print(f)
    else:
        print(f"{cfg['name']}: {len(cfg['assays'])} assay(s) -> "
              f"{', '.join(a['id'] for a in cfg['assays'])}")
