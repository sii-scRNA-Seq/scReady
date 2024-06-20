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

# Assuming input file is passed as a parameter to the script
input_file="$1"

# Check if input file is provided
if [ -z "$input_file" ]; then
    echo "Usage: $0 <input_file>"
    exit 1
fi

# Set field separator to comma
IFS=','

# Flag to skip the first line
first_line=1

# Array to hold job IDs
job_ids=()

# Read the input file line by line
while read -r line; do
    if [ $first_line -eq 1 ]; then
        first_line=0  # Skip the first line
        continue
    fi

    # Split line into fields
    fields=($line)

    # Check number of fields
    if [ ${#fields[@]} -eq 3 ]; then
        fastqs="${fields[0]}"
        ID="${fields[1]}"
        intron="${fields[2]}"

        echo "Submitting job for ID: $ID"

        # Create a SLURM script from template
        script_file="cellranger_job_${ID}.sh"
        sed -e "s|#ID|${ID}|g; s|#FASTQS|${fastqs}|g; s|#INTRON|${intron}|g" cellranger_job_template.sh > "$script_file"

        # Convert the script to Unix line endings
        dos2unix "$script_file"

        # Submit the SLURM job and capture the job ID
        job_id=$(sbatch "$script_file" | awk '{print $4}')
        job_ids+=("$job_id")
    else
        echo "Error: This line does not have exactly 3 columns. Line: $line"
    fi
done < "$input_file"

# Submit R script job with dependencies
dependency_list=$(IFS=:; echo "${job_ids[*]}")
sbatch --account=project0001 --dependency=afterok:$dependency_list --wrap="Rscript /users/ds286q/project0001/Dom/pipeline/SeuratGeneration.R /users/ds286q/project0001/Dom/pipeline/Rcodeoutput.txt"
