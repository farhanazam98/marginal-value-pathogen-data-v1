#!/usr/bin/env bash
# Driver for the (year x bit-score-threshold) grid sweep. A thin wrapper over the
# unchanged run_sweep.sh: for each threshold it exports BITSCORE_PER_RESIDUE and
# runs one year-sweep whose specs carry a per-cell tag <year>_t<thr>, so every
# cell gets its own sandbox (data/sweep/<protein>/<year>_t<thr>) and no cell
# overwrites another.
#
# The threshold override flows through config.load_config() (which reads
# BITSCORE_PER_RESIDUE), so both the jackhmmer search and the PSSM-reuse
# fingerprint see it; the env var set on the run_sweep.sh call is inherited by
# its run_year.sh children and their python steps -- no re-export needed.
#
# Thresholds run SEQUENTIALLY: run_sweep.sh holds a per-protein PID lock that
# forbids two concurrent sweeps under one root. Years run CONCURRENTLY within a
# threshold (-j), which already saturates the box since jackhmmer pins ~2 cores.
#
# Do NOT set SWEEP_ROOT: every cell must land in run_sweep.sh's default
# data/sweep root so collect.py can rebuild one merged sweep_results.csv.
#
# Usage: run_threshold_sweep.sh [-j N] [-t "0.01 0.03 0.05 0.1 0.3 0.5"] <year> [<year> ...]
#   -t "..."  space-separated bit-score-per-residue thresholds. Optional: the
#             default is EVEREST's exact grid (Methods A.6.1), which includes the
#             config baseline 0.3 (so its _t0.3 cell re-derives the plain year
#             sweep) and the low end where a sparse family clears the Neff/L depth
#             floor. COST: the default runs 6 thresholds SEQUENTIALLY (the PID
#             lock forbids concurrent sweeps), so an ordinary run is ~6x a single
#             year sweep -- pass -t to narrow it. NOTE this is 6 independent DB
#             scans; the "threshold changes are free" gotcha is per-search only.
#   -j N      max concurrent pipelines per threshold (default 6; passed to run_sweep.sh)
#
# Runs in the foreground; launch it detached (setsid nohup ... &) to survive the
# session ending. Progress: tail -f logs/sweep/<protein>/<year>_t<thr>.log
set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

JOBS=6
THRESHOLDS="0.01 0.03 0.05 0.1 0.3 0.5"   # EVEREST-exact default grid; -t overrides
usage() { echo 'usage: run_threshold_sweep.sh [-j N] [-t "0.01 0.03 0.05 0.1 0.3 0.5"] <year> [<year> ...]' >&2; exit 2; }

while getopts ":j:t:" opt; do
  case "$opt" in
    j) JOBS="$OPTARG" ;;
    t) THRESHOLDS="$OPTARG" ;;
    *) usage ;;
  esac
done
shift $((OPTIND - 1))

[ -n "$THRESHOLDS" ] && [ "$#" -gt 0 ] || usage   # default is non-empty; rejects only an explicit -t ""

for thr in $THRESHOLDS; do
  specs=()
  for y in "$@"; do
    specs+=("$y:${y}_t${thr}")
  done
  echo "=== threshold $thr : ${#specs[@]} year(s) ==="
  BITSCORE_PER_RESIDUE="$thr" bash "$DIR/run_sweep.sh" -j "$JOBS" "${specs[@]}"
done
