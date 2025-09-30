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
module load apps/R/4.4.1
module load apps/miniforge

############# MY CODE #############
#
# RUN:
# sbatch this_script.sh folder_path <optional_metadata_file.csv>
#
# R directory containing SeuratGeneration.R needs to be in the same folder of this script 
#
# folder_path need to be a folder containing cellranger mappings

folder="$1"
optional_csv_file="$2"

conda activate pipeline

# Check if path is provided
#if [ -z "$input_file" ]; then
#    echo "Usage: $0 <path> [<optional_metadata_file.csv>]"
#    exit 1
#fi

if [ -z "$optional_csv_file" ]; then
    echo Rscript $seurat_script $folder
    Rscript R/SeuratGeneration.R $folder
else
    echo Rscript $seurat_script $folder $optional_csv_file
    Rscript R/SeuratGeneration.R $folder $optional_csv_file
fi
