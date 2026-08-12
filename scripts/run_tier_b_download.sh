#!/usr/bin/env bash
# Launches Tier B (UniRef100, 2020/2022/2024/2026, ~375 GB / 16.1h per README Phase 3)
# as a detached background job so it survives the SSH session ending.
#
# Output goes to data/snapshots, which is symlinked to /mnt/scratch -- the
# root volume only has ~44G free and Tier B alone won't fit on it.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUTPUT_DIR="$REPO_ROOT/data/snapshots"
LOG_DIR="$REPO_ROOT/logs"
LOG_FILE="$LOG_DIR/tier_b_download_$(date +%Y%m%d_%H%M%S).log"
PID_FILE="$REPO_ROOT/data/snapshots/.tier_b_download.pid"

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

if [ -f "$PID_FILE" ] && kill -0 "$(cat "$PID_FILE")" 2>/dev/null; then
  echo "Tier B download already running (PID $(cat "$PID_FILE")). Not starting a second copy." >&2
  exit 1
fi

if command -v conda >/dev/null 2>&1; then
  eval "$(conda shell.bash hook)"
  conda activate marginal-value-pathogen-data
  PYTHON=python
else
  # No conda env set up on this box (checked: not installed) -- fall back to
  # the system interpreter, which already has `requests` available.
  PYTHON=python3
fi

cd "$REPO_ROOT"
nohup "$PYTHON" scripts/download_uniref100.py \
  --years 2020 2022 2024 2026 \
  --output-dir "$OUTPUT_DIR" \
  --download-workers 4 \
  --parse-workers 3 \
  > "$LOG_FILE" 2>&1 &

echo $! > "$PID_FILE"
disown

echo "Tier B download started in background (PID $(cat "$PID_FILE"))."
echo "Output dir: $OUTPUT_DIR"
echo "Log file:   $LOG_FILE"
echo
echo "Follow progress:   tail -f $LOG_FILE"
echo "Check it's alive:  kill -0 \$(cat $PID_FILE) && echo running"
echo "Stop it:           kill \$(cat $PID_FILE)"
