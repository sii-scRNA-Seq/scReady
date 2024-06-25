#!/bin/bash

############# SLURM SETTINGS #############
#SBATCH --account=project0001
#SBATCH --job-name=R_pipeline
#SBATCH --output=R_pipeline-%j.out
#SBATCH --error=R_pipeline-%j.err
#SBATCH --partition=nodes
#SBATCH --time=1-00:00:00
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1

############# LOADING MODULES (optional) #############
module load apps/cellranger
module load apps/R/4.3.0

############# MY CODE #############
#
# RUN:
# sbatch this_script.sh folder_path <optional_metadata_file.csv>
#
# SeuratGeneration.R needs to be in the same folder 
#
# folder_path need to be a folder containing cellranger mappings

folder="$1"
optional_csv_file="$2"

## Check if path is provided
#if [ -z "$input_file" ]; then
#    echo "Usage: $0 <path> [<optional_metadata_file.csv>]"
#    exit 1
#fi

if [ -z "$optional_csv_file" ]; then
    echo sbatch --account=project0001 --wrap="Rscript $seurat_script $folder"
    sbatch --account=project0001 --wrap="Rscript SeuratGeneration.R $folder"
else
    echo sbatch --account=project0001 --wrap="Rscript $seurat_script $folder $optional_csv_file"
    sbatch --account=project0001 --wrap="Rscript SeuratGeneration.R $folder $optional_csv_file"
fi

