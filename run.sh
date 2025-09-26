#!/usr/bin/env bash

folder="$1"
optional_csv_file="$2"

if [ -z "$optional_csv_file" ]; then
    Rscript SeuratGeneration.R $folder
else
    Rscript SeuratGeneration.R $folder $optional_csv_file
fi
