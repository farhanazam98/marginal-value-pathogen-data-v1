#!/usr/bin/env bash
# Runs the PSSM pipeline against ONE UniRef100 snapshot year, inside an isolated
# sandbox directory. If the sandbox already holds a PSSM built for this same
# protein, the expensive search (steps 00-04) is skipped and only scoring
# (05-06) re-runs -- so adding or changing DMS assays costs seconds, not a
# re-search, and needs no snapshot on disk.
#
# Why a sandbox: every script under scripts/pssm_pipeline/ resolves its
# checkpoint paths (data/pssm_pipeline/...) relative to the working directory,
# with no year/tag parameter -- so each concurrent run needs its own working
# directory. A sandbox whose data/ holds the right symlinks retargets the
# whole pipeline without editing a tracked script, and keeps runs comparable
# (query + threshold come from the active PROTEIN_CONFIG).
#
# The database itself doesn't need a symlink: SEQ_DB (exported below) points
# 01_jackhmmer_search.py straight at this run's snapshot, so concurrent runs
# never share that mutable the way the old repoint-the-shared-symlink
# workflow in CLAUDE.local.md did.
#
# Usage: run_year.sh <year> [tag]
#   tag defaults to <year>; pass one explicitly to run the same year more than
#   once (e.g. the concurrency probe runs 2010 six times as probe2010_1..6).
set -euo pipefail

YEAR="${1:?usage: run_year.sh <year> [tag]}"
TAG="${2:-$YEAR}"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SWEEP_ROOT="${SWEEP_ROOT:-$REPO_ROOT/data/sweep}"
RUN_DIR="$SWEEP_ROOT/$TAG"
SNAPSHOT="$REPO_ROOT/data/snapshots/uniref100_${YEAR}_01.fasta"
STATUS_FILE="$RUN_DIR/STATUS"

# Which protein to run. Points steps 00/01/05/06 at a query, threshold, and the
# list of DMS assays to score (see config/spike.yaml). One protein per sweep --
# sandboxes are keyed by year only, so running a second protein reuses these
# dirs and overwrites the first. Run proteins serially, collecting results
# between runs.
PROTEIN_CONFIG="${PROTEIN_CONFIG:-config/spike.yaml}"

BUILD_STEPS=(00_setup 01_jackhmmer_search 02_clean_msa 03_weights 04_pssm)  # need the snapshot
SCORE_STEPS=(05_score 06_evaluate)                                          # need only the PSSM

# Activate conda first: reading the protein config below needs the env's python
# (pyyaml), and the pipeline steps need jackhmmer + scipy.
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

# Preflight: confirm the active env can actually load the config. Every step and
# both config.py calls below do `import yaml`, but those calls sit inside an
# `if ... &&` and a `$(...)` assignment, so under `set -e` a missing dep is
# swallowed -- the run then silently produces nothing and the summary still
# prints a stale STATUS. Fail loudly here instead. (A missing pyyaml in the env,
# despite environment.yml listing it, is exactly what wedged the 2026-08-17 rerun.)
if ! python -c "import yaml" 2>/dev/null; then
  echo "[$TAG] active env ($CONDA_DEFAULT_ENV) cannot 'import yaml' -- config.py will fail." >&2
  echo "[$TAG] fix: conda env update -n marginal-value-pathogen-data -f environment.yml" >&2
  mkdir -p "$RUN_DIR"
  echo "FAILED:preflight" > "$STATUS_FILE"
  exit 1
fi

# Reuse an existing PSSM built for this same protein (query + threshold),
# skipping the search; scoring always re-runs. A different protein/threshold
# fails the fingerprint and rebuilds from scratch.
BUILD_META="$RUN_DIR/data/pssm_pipeline/msa_raw_run_meta.json"
if [ -f "$RUN_DIR/data/pssm_pipeline/pssm.npy" ] && \
   (cd "$REPO_ROOT" && PROTEIN_CONFIG="$PROTEIN_CONFIG" \
      python scripts/pssm_pipeline/config.py --matches-build "$BUILD_META"); then
  REUSE_PSSM=1
  STEPS=("${SCORE_STEPS[@]}")
  echo "[$TAG] reusing existing PSSM; skipping search (00-04), re-scoring only."
else
  REUSE_PSSM=0
  STEPS=("${BUILD_STEPS[@]}" "${SCORE_STEPS[@]}")
  # A 0-byte snapshot is a failed download, not a database -- data/snapshots holds
  # one (2020, from the killed Tier B run). 01's Path(SEQ_DB).exists() check would
  # happily pass on it and produce a silently empty result, so guard on size.
  # Only matters when building; a reuse run needs no snapshot on disk.
  if [ ! -s "$SNAPSHOT" ]; then
    echo "[$TAG] snapshot missing or empty: $SNAPSHOT" >&2
    mkdir -p "$RUN_DIR"
    echo "FAILED:snapshot" > "$STATUS_FILE"
    exit 1
  fi
fi

mkdir -p "$RUN_DIR/data/pssm_pipeline"
ln -sfn "$REPO_ROOT/scripts" "$RUN_DIR/scripts"
ln -sfn "$REPO_ROOT/config"  "$RUN_DIR/config"
# Symlink exactly the inputs the active protein config names -- its query FASTA
# and each assay CSV -- instead of hardcoding one protein/DMS. config.py prints
# them repo-relative; recreate that layout inside the sandbox so the pipeline's
# relative paths resolve here.
INPUT_FILES=$(cd "$REPO_ROOT" && PROTEIN_CONFIG="$PROTEIN_CONFIG" \
  python scripts/pssm_pipeline/config.py --input-files)
while IFS= read -r f; do
  mkdir -p "$RUN_DIR/$(dirname "$f")"
  ln -sfn "$REPO_ROOT/$f" "$RUN_DIR/$f"
done <<< "$INPUT_FILES"
export SEQ_DB="$SNAPSHOT"
export PROTEIN_CONFIG

# Record build provenance only when we actually build -- a reuse run keeps the
# existing sweep_run.json (its real snapshot bytes) intact.
if [ "$REUSE_PSSM" -eq 0 ]; then
  cat > "$RUN_DIR/data/pssm_pipeline/sweep_run.json" <<EOF
{
  "tag": "$TAG",
  "year": $YEAR,
  "snapshot": "$SNAPSHOT",
  "snapshot_bytes": $(stat -c %s "$SNAPSHOT"),
  "protein_config": "$PROTEIN_CONFIG",
  "started": "$(date +%Y-%m-%dT%H:%M:%S%z)"
}
EOF
fi

echo "RUNNING" > "$STATUS_FILE"

cd "$RUN_DIR"
START_EPOCH=$(date +%s)
for step in "${STEPS[@]}"; do
  echo "=== [$TAG] $step  $(date +%Y-%m-%dT%H:%M:%S%z) ==="
  if ! python "scripts/pssm_pipeline/${step}.py"; then
    echo "FAILED:$step" > "$STATUS_FILE"
    echo "[$TAG] FAILED at $step" >&2
    exit 1
  fi
done
ELAPSED=$(( $(date +%s) - START_EPOCH ))

echo "DONE" > "$STATUS_FILE"
echo "=== [$TAG] DONE in ${ELAPSED}s  $(date +%Y-%m-%dT%H:%M:%S%z) ==="
