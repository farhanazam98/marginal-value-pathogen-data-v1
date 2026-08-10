#!/usr/bin/env bash
# Runs the full PSSM pipeline (steps 00-06) against ONE UniRef100 snapshot year,
# inside an isolated sandbox directory.
#
# Why a sandbox: every script under scripts/pssm_pipeline/ hardcodes its paths
# (no argparse, no env vars) and resolves them relative to the working directory.
# So a sandbox whose data/ holds the right symlinks retargets the whole pipeline
# without editing a single tracked script -- which keeps BITSCORE_PER_RESIDUE at
# 0.3 for every year and keeps this run comparable to the existing 2010 result.
#
# data/uniref100_2010.fasta keeps that literal name in *every* sandbox, whatever
# the year, because SEQ_DB in 01_jackhmmer_search.py hardcodes it. Only the
# symlink target varies. This also means each run owns its own copy of that name,
# so years can run concurrently -- impossible with the repoint-the-shared-symlink
# workflow in CLAUDE.local.md, which is a single global mutable.
#
# Usage: run_year.sh <year> [tag]
#   tag defaults to <year>; pass one explicitly to run the same year more than
#   once (e.g. the concurrency probe runs 2010 six times as probe2010_1..6).
set -euo pipefail

YEAR="${1:?usage: run_year.sh <year> [tag]}"
TAG="${2:-$YEAR}"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SWEEP_ROOT="${SWEEP_ROOT:-/mnt/scratch/sweep}"
RUN_DIR="$SWEEP_ROOT/$TAG"
SNAPSHOT="$REPO_ROOT/data/snapshots/uniref100_${YEAR}_01.fasta"
STATUS_FILE="$RUN_DIR/STATUS"

STEPS=(00_setup 01_jackhmmer_search 02_clean_msa 03_weights 04_pssm 05_score 06_evaluate)

# Idempotent: a completed run is never redone, so re-running the driver after a
# crash resumes the sweep instead of restarting it.
if [ -f "$STATUS_FILE" ] && [ "$(cat "$STATUS_FILE")" = "DONE" ]; then
  echo "[$TAG] already DONE, skipping."
  exit 0
fi

# A 0-byte snapshot is a failed download, not a database -- data/snapshots holds
# one (2020, from the killed Tier B run). 01's Path(SEQ_DB).exists() check would
# happily pass on it and produce a silently empty result, so guard on size.
if [ ! -s "$SNAPSHOT" ]; then
  echo "[$TAG] snapshot missing or empty: $SNAPSHOT" >&2
  mkdir -p "$RUN_DIR"
  echo "FAILED:snapshot" > "$STATUS_FILE"
  exit 1
fi

mkdir -p "$RUN_DIR/data/pssm_pipeline"
ln -sfn "$REPO_ROOT/scripts"                              "$RUN_DIR/scripts"
ln -sfn "$REPO_ROOT/data/protein.fasta"                   "$RUN_DIR/data/protein.fasta"
ln -sfn "$REPO_ROOT/data/SARS2_RBD_Starr_binding_dms.csv" "$RUN_DIR/data/SARS2_RBD_Starr_binding_dms.csv"
ln -sfn "$SNAPSHOT"                                       "$RUN_DIR/data/uniref100_2010.fasta"

# Provenance: the pipeline's own meta files record the database only as the
# generic "data/uniref100_2010.fasta", so without this a finished sandbox can't
# be traced back to the year it actually searched.
cat > "$RUN_DIR/data/pssm_pipeline/sweep_run.json" <<EOF
{
  "tag": "$TAG",
  "year": $YEAR,
  "snapshot": "$SNAPSHOT",
  "snapshot_bytes": $(stat -c %s "$SNAPSHOT"),
  "started": "$(date -Is)"
}
EOF

echo "RUNNING" > "$STATUS_FILE"

CONDA_SH="$HOME/miniconda3/etc/profile.d/conda.sh"
if [ -f "$CONDA_SH" ]; then
  # shellcheck disable=SC1090
  . "$CONDA_SH"
  conda activate marginal-value-pathogen-data
elif command -v conda >/dev/null 2>&1; then
  eval "$(conda shell.bash hook)"
  conda activate marginal-value-pathogen-data
else
  echo "[$TAG] conda env not found; jackhmmer and scipy will be missing." >&2
  exit 1
fi

cd "$RUN_DIR"
START_EPOCH=$(date +%s)
for step in "${STEPS[@]}"; do
  echo "=== [$TAG] $step  $(date -Is) ==="
  if ! python "scripts/pssm_pipeline/${step}.py"; then
    echo "FAILED:$step" > "$STATUS_FILE"
    echo "[$TAG] FAILED at $step" >&2
    exit 1
  fi
done
ELAPSED=$(( $(date +%s) - START_EPOCH ))

echo "DONE" > "$STATUS_FILE"
echo "=== [$TAG] DONE in ${ELAPSED}s  $(date -Is) ==="
