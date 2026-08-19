#!/usr/bin/env python3
"""Collect the snapshot-year sweep into a single tidy table.

Reads every finished (or partly finished) sandbox under the sweep root and joins
the meta files each pipeline step leaves behind into `data/sweep_results.csv`,
one row per run.

This is a *stateless re-derivation* from on-disk checkpoints rather than an
append-as-you-go log. That makes it correct no matter when it runs -- mid-sweep,
repeatedly, or after a session crash -- and means a partially completed sweep
still yields a valid table of whatever has landed so far. Runs still in flight
appear with their completed columns filled and the rest blank.

Every available metric is carried through, not just rho, because the
visualization is curated separately and shouldn't have to go back to the raw
checkpoints to pick its axes.

Usage: python scripts/sweep/collect.py
"""

import csv
import json
import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SWEEP_ROOT = Path(os.environ.get("SWEEP_ROOT", REPO_ROOT / "data" / "sweep"))
SNAPSHOT_DIR = REPO_ROOT / "data" / "snapshots"
OUT_CSV = REPO_ROOT / "data" / "sweep_results.csv"

# Sandboxes are keyed by (protein, year) as <protein>/<year>. Runs from before
# that layout sit flat as <year> directly under the sweep root; those predate
# multi-protein support and are all spike (config/spike.yaml), so a flat run's
# protein defaults to this.
DEFAULT_PROTEIN = "spike"

# Steps 00-04 are per-(protein, year): one MSA/PSSM, shared by every assay. Their
# metas contribute the same values to each of a run's assay rows.
# (output column, meta file, key in that file)
SHARED_FIELDS = [
    ("tag",                     "sweep_run.json",        "tag"),
    ("year",                    "sweep_run.json",        "year"),
    ("snapshot_bytes",          "sweep_run.json",        "snapshot_bytes"),

    ("bitscore_per_residue",    "msa_raw_run_meta.json", "bitscore_per_residue"),
    ("query_length",            "msa_raw_run_meta.json", "query_length"),
    ("threshold_bits",          "msa_raw_run_meta.json", "threshold_bits"),
    ("jackhmmer_elapsed_s",     "msa_raw_run_meta.json", "elapsed_seconds"),
    ("jackhmmer_rounds",        "msa_raw_run_meta.json", "rounds"),
    ("jackhmmer_converged",     "msa_raw_run_meta.json", "converged"),
    ("n_hits",                  "msa_raw_run_meta.json", "n_hits"),
    ("n_alignment_rows",        "msa_raw_run_meta.json", "n_alignment_rows"),

    ("N_raw",                   "msa_clean_meta.json",   "N_raw"),
    ("L_raw",                   "msa_clean_meta.json",   "L_raw"),
    ("N_final",                 "msa_clean_meta.json",   "N_final"),
    ("L_final",                 "msa_clean_meta.json",   "L_final"),

    ("theta",                   "weights_meta.json",     "theta"),
    ("Neff",                    "weights_meta.json",     "Neff"),
    ("Neff_over_L",             "weights_meta.json",     "depth_Neff_over_L"),
    ("clears_depth_floor",      "weights_meta.json",     "clears_depth_floor"),
    ("Neff_at_90pct_identity",  "weights_meta.json",     "Neff_at_90pct_identity"),
    ("clears_reliability",      "weights_meta.json",     "clears_reliability_threshold"),
    ("n_singleton_sequences",   "weights_meta.json",     "n_singleton_sequences"),
]

# Steps 05-06 run once per assay. The meta filename carries the assay id
# (predictions_meta_<id>.json, evaluate_meta_<id>.json); "base" here is that name
# without the id, filled in per assay in collect_run.
# (output column, meta base, key in that file)
ASSAY_FIELDS = [
    ("n_variants",              "predictions_meta", "n_variants"),
    ("n_scored_directly",       "predictions_meta", "n_scored_directly"),
    ("n_imputed",               "predictions_meta", "n_imputed"),
    ("imputed_value",           "predictions_meta", "imputed_value"),
    ("wt_wt_all_zero",          "predictions_meta", "wt_wt_all_zero"),
    ("predicted_score_mean",    "predictions_meta", "predicted_score_mean"),
    ("predicted_score_std",     "predictions_meta", "predicted_score_std"),

    ("n_joined",                "evaluate_meta",    "n_joined"),
    ("n_dropped_from_dms",      "evaluate_meta",    "n_dropped_from_dms"),
    ("spearman_rho",            "evaluate_meta",    "spearman_rho"),
    ("spearman_pvalue",         "evaluate_meta",    "spearman_pvalue"),
    ("bootstrap_ci_95_lo",      "evaluate_meta",    "bootstrap_ci_95_lo"),
    ("bootstrap_ci_95_hi",      "evaluate_meta",    "bootstrap_ci_95_hi"),
    ("spearman_rho_excl_imputed", "evaluate_meta",  "spearman_rho_excluding_imputed"),
    ("n_excl_imputed",          "evaluate_meta",    "n_excluding_imputed"),
]

DERIVED = [
    "snapshot_gb", "db_n_seqs", "db_n_residues", "imputed_frac", "status",
]


def load_json(path):
    try:
        with open(path) as fh:
            return json.load(fh)
    except (FileNotFoundError, json.JSONDecodeError):
        return {}


def discover_assays(ckpt):
    """Assay ids present in this run, from the per-assay evaluate metas. Falls
    back to [""] for old single-assay checkpoints (evaluate_meta.json, no id)."""
    ids = sorted(p.stem[len("evaluate_meta_"):] for p in ckpt.glob("evaluate_meta_*.json"))
    return ids if ids else [""]


def discover_runs(sweep_root):
    """Yield (protein, run_dir) for every sandbox, tolerating both the flat
    legacy layout (<sweep_root>/<year>) and the per-protein layout
    (<sweep_root>/<protein>/<year>). A directory holding data/pssm_pipeline is
    itself a flat run; otherwise it's a protein dir whose children are runs."""
    def is_run(d):
        return d.is_dir() and (d / "data" / "pssm_pipeline").exists()

    for child in sorted(sweep_root.iterdir()):
        if is_run(child):
            yield DEFAULT_PROTEIN, child
        elif child.is_dir():
            for run_dir in sorted(child.iterdir()):
                if is_run(run_dir):
                    yield child.name, run_dir


def collect_run(run_dir, protein):
    """One row per assay in this run: the shared step 00-04 metrics repeated,
    plus that assay's step 05-06 metrics."""
    ckpt = run_dir / "data" / "pssm_pipeline"
    shared_metas = {name: load_json(ckpt / name) for name in {f[1] for f in SHARED_FIELDS}}

    base = {col: shared_metas.get(fname, {}).get(key, "")
            for col, fname, key in SHARED_FIELDS}
    base["protein"] = protein
    base["tag"] = base["tag"] or run_dir.name

    status_file = run_dir / "STATUS"
    base["status"] = status_file.read_text().strip() if status_file.exists() else "UNKNOWN"

    if base["snapshot_bytes"] != "":
        base["snapshot_gb"] = round(base["snapshot_bytes"] / 1e9, 2)

    # Sequence counts come from the download step's stats file, which is the only
    # place the database's own size in sequences/residues is recorded.
    year = base.get("year")
    if year:
        stats = load_json(SNAPSHOT_DIR / f"uniref100_{year}_01.stats.json")
        base["db_n_seqs"] = stats.get("n_seqs", "")
        base["db_n_residues"] = stats.get("n_residues", "")

    rows = []
    for assay_id in discover_assays(ckpt):
        suffix = f"_{assay_id}" if assay_id else ""
        assay_metas = {b: load_json(ckpt / f"{b}{suffix}.json")
                       for b in {f[1] for f in ASSAY_FIELDS}}
        row = dict(base)
        row["dms_id"] = assay_id
        for col, mbase, key in ASSAY_FIELDS:
            row[col] = assay_metas.get(mbase, {}).get(key, "")
        if row["n_imputed"] != "" and row["n_variants"] not in ("", 0):
            row["imputed_frac"] = round(row["n_imputed"] / row["n_variants"], 4)
        rows.append(row)

    return rows


def main():
    if not SWEEP_ROOT.exists():
        raise SystemExit(f"sweep root not found: {SWEEP_ROOT}")

    runs = list(discover_runs(SWEEP_ROOT))
    if not runs:
        raise SystemExit(f"no runs found under {SWEEP_ROOT}")

    rows = [row for protein, d in runs for row in collect_run(d, protein)]
    rows.sort(key=lambda r: (r.get("protein") or "", r.get("year") or 0,
                             r.get("tag") or "", r.get("dms_id") or ""))

    shared_cols = [f[0] for f in SHARED_FIELDS]
    columns = ["protein", "tag", "dms_id"] + [c for c in shared_cols if c != "tag"] \
        + [f[0] for f in ASSAY_FIELDS] + DERIVED
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_CSV, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({c: row.get(c, "") for c in columns})

    print(f"Wrote {OUT_CSV}  ({len(rows)} rows)\n")
    hdr = f"{'protein':<12}{'tag':<14}{'dms_id':<18}{'status':<10}{'GB':>7}{'hits':>7}{'N':>7}{'L':>7}{'Neff/L':>9}{'imp%':>7}{'rho':>8}{'rho_ni':>8}{'jh_s':>8}"
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        def fmt(key, spec, scale=1):
            v = r.get(key, "")
            if isinstance(v, (int, float)):
                return format(v * scale, spec)
            # Pad the placeholder to the column's width so in-flight runs, which
            # have most fields still empty, stay aligned with finished ones.
            return format("-", f">{spec.split('.')[0]}")
        print(
            f"{str(r.get('protein', '')):<12}{str(r['tag']):<14}{str(r.get('dms_id', '')):<18}{str(r['status']):<10}"
            f"{fmt('snapshot_gb', '7.1f')}"
            f"{fmt('n_hits', '7.0f')}"
            f"{fmt('N_final', '7.0f')}"
            f"{fmt('L_final', '7.0f')}"
            f"{fmt('Neff_over_L', '9.3f')}"
            f"{fmt('imputed_frac', '7.1f', 100)}"
            f"{fmt('spearman_rho', '8.4f')}"
            f"{fmt('spearman_rho_excl_imputed', '8.4f')}"
            f"{fmt('jackhmmer_elapsed_s', '8.0f')}"
        )


if __name__ == "__main__":
    main()
