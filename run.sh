#!/usr/bin/env bash

folder="$1"
optional_csv_file="$2"

# Pipeline location
SCRIPT_DIR="/opt/app/R"
# Default config file
DEFAULT_CFG="/opt/app/config/scReady.config"

print_help() {
    cat <<EOF
scReady pipeline

Usage:
    docker run [docker args] <CONTAINER> input_dir [optional_metadata_file]
    docker run [docker args] <CONTAINER> --print-default-config
    docker run [docker args] <CONTAINER> init-config

Helpers:
    --print-default-config  Print built-in config to stdout.
    init-config             Write template config file to current directory.
                            Refuses to overwrite existing file.

Notes:
    Example [docker args] for local runs:
      -u "$(id -u)":"$(id -g)" -v "$PWD":/work    
EOF
}

# Helper: --print-default-config
if [[ "${1:-}" == "--print-default-config" ]]; then
  cat "$DEFAULT_CFG"
  exit 0
fi

# Helper: init-config
if [[ "${1:-}" == "init-config" ]]; then
  target="/work/scReady.config"
  mkdir -p "$(dirname "$target")"
  if [[ -e "$target" ]]; then
    echo "Refusing to overwrite existing file: $(basename "$target")" >&2
    exit 1
  fi
  cp "$DEFAULT_CFG" "$target"
  # ensure host user owns it when run with -u
  chown "$(id -u)":"$(id -g)" "$target" 2>/dev/null || true
  echo "Wrote $(basename "$target")"
  exit 0
fi

# Print help
if [[ $# -eq 0 || "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
  print_help
  exit 0
fi

# Run the pipeline
if [[ -f "R/SeuratGeneration.R" ]]; then
    SCRIPT="R/SeuratGeneration.R"
else
    SCRIPT="${SCRIPT_DIR}/SeuratGeneration.R"
fi

if [ -z "$optional_csv_file" ]; then
    Rscript "${SCRIPT}" $folder
else
    Rscript "${SCRIPT}" $folder $optional_csv_file
fi
