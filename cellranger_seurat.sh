#!/bin/bash

############# SLURM SETTINGS #############
#SBATCH --account=project0001
#SBATCH --job-name=pipeline
#SBATCH --output=pipeline-%j.out
#SBATCH --error=pipeline-%j.err
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
# sbatch this_script.sh lib.csv <optional_metadata_file.csv>
#
# cellranger_job_template.sh needs to be in the same folder
# SeuratGeneration.R needs to be in the same folder 
#
# lib.csv has:
# fastqs,sample,count introns,generate bam,output folder
# /path/to/fastq/files/,FASTQID1,true,false,/path/to/output/folder1
# /path/to/fastq/files/,FASTQID2,false,true,/path/to/output/folder2
# /path/to/fastq/files/,FASTQID3,true,true,/path/to/output/folder3

# Assuming input file is passed as a parameter to the script
input_file="$1"
optional_csv_file="$2"

# Check if input file is provided
if [ -z "$input_file" ]; then
    echo "Usage: $0 <input_file> [<optional_metadata_file.csv>]"
    exit 1
fi

# Set field separator to comma
IFS=','

# Flag to skip the first line
first_line=1

# Array to hold job IDs
job_ids=()

# Array to hold unique output folders
declare -A unique_folders

# Read the input file line by line
while read -r line; do
    if [ $first_line -eq 1 ]; then
        first_line=0  # Skip the first line
        continue
    fi

    # Split line into fields
    fields=($line)

    # Check number of fields
    if [ ${#fields[@]} -eq 5 ]; then
        fastqs="${fields[0]}"
        ID="${fields[1]}"
        intron="${fields[2]}"
	bam="${fields[3]}"
	output="${fields[4]}"

        echo "FASTQS: $fastqs"
        echo "ID: $ID"
        echo "INTRON: $intron"
	echo "BAM: $bam"
        echo "OUTPUT: $output"

        echo "Submitting job for ID: $ID"

        ## Create a SLURM script from template
        script_file="cellranger_job_${ID}.sh"
        sed -e "s|#ID|${ID}|g; s|#FASTQS|${fastqs}|g; s|#INTRON|${intron}|g; s|#BAM|${bam}|g; s|#OUTPUT|${output}|g" cellranger_job_template.sh > "$script_file"
	
        # Convert the script to Unix line endings
        # dos2unix "$script_file"

        # Submit the SLURM job and capture the job ID
        job_id=$(sbatch "$script_file" | awk '{print $4}')
        job_ids+=("$job_id")

        # Collect unique output folders
        unique_folders["$output"]=1

    else
        echo "Error: This line does not have exactly 5 columns. Line: $line"
    fi
done < "$input_file"

# Submit R script jobs for each unique output folder with dependencies
dependency_list=$(IFS=:; echo "${job_ids[*]}")
for folder in "${!unique_folders[@]}"; do
    mkdir -p $folder
    if [ -z "$optional_csv_file" ]; then
        sbatch --account=project0001 --dependency=afterok:$dependency_list --wrap="Rscript SeuratGeneration.R $folder"
    else
        sbatch --account=project0001 --dependency=afterok:$dependency_list --wrap="Rscript SeuratGeneration.R $folder $optional_csv_file"
    fi
done
