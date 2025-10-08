#!/bin/bash

############# SLURM SETTINGS #############
#SBATCH --account=project0001
#SBATCH --job-name=R_pipeline
#SBATCH --output=R_pipeline-%j.out
#SBATCH --error=R_pipeline-%j.err
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --ntasks-per-node=1

############# LOADING MODULES (optional) #############
module load apps/apptainer

############# MY CODE #############
#
# RUN:
# sbatch this_script.sh folder_path <optional_metadata_file.csv>
#
# folder_path need to be a folder containing cellranger mappings

IMAGE="docker://ghcr.io/sii-scrna-seq/scready:latest"

folder="$1"
optional_csv_file="$2"

if [ -z "$optional_csv_file" ]; then
    apptainer run "$IMAGE" $folder
else
    apptainer run "$IMAGE" $folder $optional_csv_file
fi
