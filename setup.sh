#!/usr/bin/env bash
# Bootstraps everything this repo needs on a fresh machine (e.g. a new EC2 instance):
# installs Miniconda if missing, accepts the conda channel Terms of Service (required
# non-interactively before any install), then creates/updates the `marginal-value-pathogen-data`
# conda environment from environment.yml, which provides the `datasets` / `dataformat`
# CLI (ncbi-datasets-cli) used by scripts/build_temporal_snapshots.py and friends.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONDA_DIR="${CONDA_DIR:-$HOME/miniconda3}"
ENV_NAME="marginal-value-pathogen-data"

os="$(uname -s)"
arch="$(uname -m)"

case "$os" in
  Darwin) os_tag="MacOSX" ;;
  Linux)  os_tag="Linux" ;;
  *) echo "Unsupported OS: $os" >&2; exit 1 ;;
esac

case "$arch" in
  arm64)   arch_tag="arm64" ;;
  aarch64) arch_tag="aarch64" ;;
  x86_64)  arch_tag="x86_64" ;;
  *) echo "Unsupported architecture: $arch" >&2; exit 1 ;;
esac
# Miniconda's macOS arm64 installer is named "arm64"; its Linux arm64 installer is named "aarch64".
if [ "$os_tag" = "MacOSX" ] && [ "$arch_tag" = "aarch64" ]; then arch_tag="arm64"; fi
if [ "$os_tag" = "Linux" ] && [ "$arch_tag" = "arm64" ]; then arch_tag="aarch64"; fi

if command -v conda >/dev/null 2>&1; then
  CONDA_BIN="$(command -v conda)"
elif [ -x "$CONDA_DIR/bin/conda" ]; then
  CONDA_BIN="$CONDA_DIR/bin/conda"
else
  installer="Miniconda3-latest-${os_tag}-${arch_tag}.sh"
  url="https://repo.anaconda.com/miniconda/${installer}"
  echo "conda not found -- installing Miniconda from $url"
  tmp_dir="$(mktemp -d)"
  trap 'rm -rf "$tmp_dir"' EXIT
  curl -fsSL "$url" -o "$tmp_dir/$installer"
  bash "$tmp_dir/$installer" -b -p "$CONDA_DIR"
  CONDA_BIN="$CONDA_DIR/bin/conda"
fi

echo "Accepting conda default-channel Terms of Service (required before installing)..."
"$CONDA_BIN" tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
"$CONDA_BIN" tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

if "$CONDA_BIN" env list | grep -qE "^${ENV_NAME}\s"; then
  echo "Updating existing '${ENV_NAME}' environment..."
  "$CONDA_BIN" env update -n "$ENV_NAME" -f "$REPO_ROOT/environment.yml" --prune
else
  echo "Creating '${ENV_NAME}' environment..."
  "$CONDA_BIN" env create -f "$REPO_ROOT/environment.yml"
fi

echo
echo "Verifying installed CLI tools..."
"$CONDA_BIN" run -n "$ENV_NAME" datasets version
"$CONDA_BIN" run -n "$ENV_NAME" dataformat --help >/dev/null && echo "dataformat OK"

echo
echo "Setup complete. Activate with:"
echo "  conda activate ${ENV_NAME}"
