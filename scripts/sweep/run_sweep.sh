#!/usr/bin/env bash
# Driver for the snapshot-year sweep: runs run_year.sh for each requested year
# with a concurrency cap, longest job first.
#
# Longest-first matters: the jobs are independent and jackhmmer pins only ~2
# effective cores each, so wall time collapses to roughly the single longest job
# -- but only if that job starts immediately. Started last, it becomes a tail.
#
# Concurrency is capped rather than unbounded on purpose: run_tier_b_download.sh
# previously wedged this instance by oversubscribing it.
#
# Usage: run_sweep.sh [-j N] <spec> [<spec> ...]
#   spec is YEAR or YEAR:TAG   (YEAR:TAG lets one year run several times, as the
#                               concurrency probe does with 2010)
#   -j N   max concurrent pipelines (default 6: 16 vCPUs / ~2 effective cores)
#
# Runs in the foreground; launch it detached (setsid nohup ... &) to survive the
# session ending. Progress: tail -f logs/sweep/<tag>.log
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SWEEP_ROOT="${SWEEP_ROOT:-$REPO_ROOT/data/sweep}"

# Key the lock, logs, and sandboxes by protein (config filename stem) so a second
# protein's sweep can't collide with the first. Export the config so run_year.sh
# derives the identical tag.
export PROTEIN_CONFIG="${PROTEIN_CONFIG:-config/spike.yaml}"
PROTEIN_TAG="$(basename "${PROTEIN_CONFIG%.*}")"   # config/spike.yaml -> spike
LOG_DIR="$REPO_ROOT/logs/sweep/$PROTEIN_TAG"
PID_FILE="$SWEEP_ROOT/$PROTEIN_TAG/.run_sweep.pid"
JOBS=6

while getopts ":j:" opt; do
  case "$opt" in
    j) JOBS="$OPTARG" ;;
    *) echo "usage: run_sweep.sh [-j N] <year|year:tag> ..." >&2; exit 2 ;;
  esac
done
shift $((OPTIND - 1))

[ "$#" -gt 0 ] || { echo "usage: run_sweep.sh [-j N] <year|year:tag> ..." >&2; exit 2; }

mkdir -p "$LOG_DIR" "$(dirname "$PID_FILE")"

if [ -f "$PID_FILE" ] && kill -0 "$(cat "$PID_FILE")" 2>/dev/null; then
  echo "Sweep already running (PID $(cat "$PID_FILE")). Not starting a second copy." >&2
  exit 1
fi
echo $$ > "$PID_FILE"
trap 'rm -f "$PID_FILE"' EXIT

# Order specs by snapshot size, descending, so the long pole starts first.
SORTED=$(
  for spec in "$@"; do
    year="${spec%%:*}"
    snap="$REPO_ROOT/data/snapshots/uniref100_${year}_01.fasta"
    size=$(stat -c %s "$snap" 2>/dev/null || echo 0)
    echo "$size $spec"
  done | sort -rn | awk '{print $2}'
)

echo "Sweep starting $(date +%Y-%m-%dT%H:%M:%S%z)"
echo "  concurrency: $JOBS"
echo "  order (largest DB first): $(echo "$SORTED" | tr '\n' ' ')"
echo

for spec in $SORTED; do
  year="${spec%%:*}"
  tag="${spec##*:}"   # equals year when spec has no colon
  # Throttle by polling rather than `wait -n`: the latter is bash 4.3+ only and
  # aborts under `set -e` on the macOS system bash (3.2).
  while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do
    sleep 0.5
  done
  echo "launching $tag (year $year) -> $LOG_DIR/$tag.log"
  bash "$REPO_ROOT/scripts/sweep/run_year.sh" "$year" "$tag" \
    > "$LOG_DIR/$tag.log" 2>&1 &
done

wait
echo
echo "Sweep finished $(date +%Y-%m-%dT%H:%M:%S%z)"
for spec in $SORTED; do
  tag="${spec##*:}"
  printf '  %-16s %s\n' "$tag" "$(cat "$SWEEP_ROOT/$PROTEIN_TAG/$tag/STATUS" 2>/dev/null || echo MISSING)"
done
