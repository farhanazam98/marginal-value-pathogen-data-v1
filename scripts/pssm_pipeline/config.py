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

import json
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

    # Opt-in override of the bit-score threshold without editing the YAML, used
    # by the threshold sweep (scripts/sweep/run_threshold_sweep.sh). Because both
    # the jackhmmer search and the PSSM-reuse fingerprint (matches_build) read the
    # threshold through this function, overriding it here keeps the search and the
    # reuse check consistent: a run built at a different threshold fails the
    # fingerprint and rebuilds. Unset -> unchanged behavior (year sweep, manual runs).
    # Treat empty/whitespace as unset (a common leftover `export BITSCORE_PER_RESIDUE=`
    # must not crash every step, and a plain year sweep must keep the YAML value),
    # and fail loudly on a non-numeric value instead of an opaque float() traceback.
    override = os.environ.get("BITSCORE_PER_RESIDUE")
    if override is not None and override.strip():
        try:
            cfg["bitscore_per_residue"] = float(override)
        except ValueError:
            raise SystemExit(f"BITSCORE_PER_RESIDUE must be numeric, got: {override!r}")
    return cfg


def input_files(cfg):
    """Every on-disk input the config points at: the query plus each assay CSV.
    These are what the sweep sandbox must symlink in."""
    return [cfg["query_fasta"]] + [a["csv"] for a in cfg["assays"]]


def matches_build(cfg, meta_path):
    """True if an existing MSA/PSSM build (its msa_raw_run_meta.json) was made
    with this config's query and threshold, so its PSSM is safe to reuse instead
    of re-running the jackhmmer search."""
    try:
        meta = json.load(open(meta_path))
    except (FileNotFoundError, json.JSONDecodeError):
        return False
    query_len = len(open(cfg["query_fasta"]).read().split("\n", 1)[1].replace("\n", ""))
    return (meta.get("query_length") == query_len
            and meta.get("bitscore_per_residue") == cfg["bitscore_per_residue"])


if __name__ == "__main__":
    cfg = load_config()
    if "--input-files" in sys.argv:
        for f in input_files(cfg):
            print(f)
    elif "--matches-build" in sys.argv:
        # Exit 0 (reuse the PSSM) only if the build fingerprint matches.
        sys.exit(0 if matches_build(cfg, sys.argv[sys.argv.index("--matches-build") + 1]) else 1)
    else:
        print(f"{cfg['name']}: {len(cfg['assays'])} assay(s) -> "
              f"{', '.join(a['id'] for a in cfg['assays'])}")
