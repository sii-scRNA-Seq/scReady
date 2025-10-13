#!/bin/bash

############# SLURM SETTINGS #############
#SBATCH --job-name=scReady
#SBATCH --output=scReady-%j.out
#SBATCH --error=scReady-%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1

############# LOADING MODULES (optional) #############
# module load apps/apptainer

############# MY CODE #############
#
# RUN:
# sbatch this_script.sh folder_path <optional_metadata_file.csv>

IMAGE="docker://ghcr.io/sii-scrna-seq/scready:latest"

folder="$1"
optional_csv_file="$2"

if [ -z "$optional_csv_file" ]; then
    apptainer run "$IMAGE" $folder
else
    apptainer run "$IMAGE" $folder $optional_csv_file
fi
