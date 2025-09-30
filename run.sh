#!/usr/bin/env bash

folder="$1"
optional_csv_file="$2"

# Container location
SCRIPT_DIR="/opt/app/R"

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
