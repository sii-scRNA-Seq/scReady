#!/bin/bash

############# SLURM SETTINGS #############
#SBATCH --account=project0001
#SBATCH --job-name=pipeline
#SBATCH --output=pipeline-%j.out
#SBATCH --error=pipeline-%j.err
#SBATCH --partition=short
#SBATCH --time=0-02:00:00
#SBATCH --mem=1G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1

############# LOADING MODULES #############
module load apps/R/4.4.1

############# INPUT VALIDATION #############
input_file="$1"
feature_ref_file="$2"  # Optional parameter

# Check if input file is provided
if [ -z "$input_file" ]; then
    echo "Usage: $0 <input_file> [<feature_ref_file>]"
    exit 1
fi

# Check if input file exists
if [ ! -f "$input_file" ]; then
    echo "Error: Input file $input_file does not exist."
    exit 1
fi

# Initialize flag for feature data presence
has_feature_data=0

############# PROCESS INPUT FILE #############
IFS=','

# Arrays to hold job IDs and unique output folders
job_ids=()
declare -A unique_folders

# Read the input file line by line
first_line=1
while read -r line; do
    if [ $first_line -eq 1 ]; then
        first_line=0  # Skip the header
        continue
    fi

    # Split line into fields
    fields=($line)
    num_fields=${#fields[@]}

    # Validate number of fields (5, 6, or 8 columns)
    if [ $num_fields -ne 5 ] && [ $num_fields -ne 6 ] && [ $num_fields -ne 8 ]; then
        echo "Error: Line must have 5, 6, or 8 columns. Line: $line"
        continue
    fi

    # Extract fields (common to 5, 6, and 8 columns)
    gex_fastqs="${fields[0]}"
    gex_fastq_id="${fields[1]}"
    intron="${fields[2]}"
    bam="${fields[3]}"
    output="${fields[4]}"

    # Initialize optional fields
    souporcell_param=""
    feature_fastqs=""
    feature_fastq_id=""

    # Extract feature fields if present (8 columns)
    if [ $num_fields -eq 8 ]; then
        feature_fastqs="${fields[6]}"
        feature_fastq_id="${fields[7]}"
        # Set flag if feature data exists
        if [ -n "$feature_fastqs" ]; then
            has_feature_data=1
        fi
    fi
done < "$input_file"

# Check if feature_ref_file is required and valid
if [ $has_feature_data -eq 1 ]; then
    if [ -z "$feature_ref_file" ]; then
        echo "Error: Feature reference file is required for samples with feature data."
        echo "Usage: $0 <input_file> [<feature_ref_file>]"
        exit 1
    fi
    
    if [ ! -f "$feature_ref_file" ]; then
        echo "Error: Feature reference file $feature_ref_file does not exist."
        exit 1
    fi
fi

# Reset file for main processing
first_line=1
while read -r line; do
    if [ $first_line -eq 1 ]; then
        first_line=0  # Skip the header
        continue
    fi

    # Split line into fields
    fields=($line)
    num_fields=${#fields[@]}

    # Skip invalid lines
    if [ $num_fields -ne 5 ] && [ $num_fields -ne 6 ] && [ $num_fields -ne 8 ]; then
        continue
    fi

    # Extract fields
    gex_fastqs="${fields[0]}"
    gex_fastq_id="${fields[1]}"
    intron="${fields[2]}"
    bam="${fields[3]}"
    output="${fields[4]}"
    souporcell_param=""
    feature_fastqs=""
    feature_fastq_id=""

    # Extract souporcell_param if present (6 or 8 columns)
    if [ $num_fields -eq 6 ] || [ $num_fields -eq 8 ]; then
        souporcell_param="${fields[5]}"
        # Validate souporcell_param is a number
        if ! [[ "$souporcell_param" =~ ^[0-9]+$ ]]; then
            echo "Error: souporcell_param must be a non-negative integer. Line: $line"
            continue
        fi
        # Force bam=true if souporcell_param is non-zero
        if [ "$souporcell_param" -ne 0 ]; then
            original_bam="$bam"
            bam="true"
            if [ "$original_bam" != "true" ]; then
                echo "Warning: Overriding bam to 'true' for sample $gex_fastq_id (souporcell_param=$souporcell_param)."
            fi
        fi
    fi

    # Validate bam field
    if [ "$bam" != "true" ] && [ "$bam" != "false" ]; then
        echo "Error: Invalid value for 'bam' (must be 'true' or 'false'). Line: $line"
        continue
    fi

    # Extract feature fields if present (8 columns)
    if [ $num_fields -eq 8 ]; then
        feature_fastqs="${fields[6]}"
        feature_fastq_id="${fields[7]}"
    fi

    # Validate paths
    if [ ! -d "$gex_fastqs" ]; then
        echo "Error: GEX FASTQ directory does not exist: $gex_fastqs"
        continue
    fi

    if [ -n "$feature_fastqs" ] && [ ! -d "$feature_fastqs" ]; then
        echo "Error: Feature Barcode FASTQ directory does not exist: $feature_fastqs"
        continue
    fi

    if [ ! -d "$output" ]; then
        echo "Creating output directory: $output"
        mkdir -p "$output" || { echo "Failed to create output directory: $output"; continue; }
    fi

    # Create libraries.csv if feature data exists
    if [ $num_fields -eq 8 ] && [ -n "$feature_fastqs" ] && [ -n "$feature_fastq_id" ]; then
        # Generate libraries.csv
        libraries_csv="$output/libraries_${gex_fastq_id}.csv"
        echo "fastqs,sample,library_type" > "$libraries_csv"
        echo "$gex_fastqs,$gex_fastq_id,Gene Expression" >> "$libraries_csv"
        echo "$feature_fastqs,$feature_fastq_id,Antibody Capture" >> "$libraries_csv"
    fi

    # Create and submit cellranger job
    script_file="cellranger_job_${gex_fastq_id}.sh"
    cat <<EOF > "$script_file"
#!/bin/bash -l

############# SLURM SETTINGS #############
#SBATCH --account=project0001
#SBATCH --job-name=cellranger_count_${gex_fastq_id}
#SBATCH --output=cellranger_count_${gex_fastq_id}-%j.out
#SBATCH --error=cellranger_count_${gex_fastq_id}-%j.err
#SBATCH --partition=nodes
#SBATCH --time=1-00:00:00
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --ntasks-per-node=1

############# LOADING MODULES #############

############# MY CODE #############

# Clean the output directory (if it exists)
if [ -d "${output}/${gex_fastq_id}" ]; then
    echo "Cleaning existing output directory: ${output}/${gex_fastq_id}"
    rm -rf "${output}/${gex_fastq_id}"
fi

# Run cellranger count
if [ $num_fields -eq 8 ] && [ -n "$feature_fastqs" ] && [ -n "$feature_fastq_id" ]; then
    cellranger count \\
        --id=${gex_fastq_id} \\
        --libraries=${libraries_csv} \\
        --feature-ref=${feature_ref_file} \\
        --transcriptome=/users/ds286q/project0001/10XReference/Human/refdata-gex-GRCh38-2020-A \\
        --include-introns=${intron} \\
        --r1-length=26 \\
        --localcores=20 \\
        --localmem=40 \\
        --create-bam=${bam} \\
        --output-dir=${output}/${gex_fastq_id}
else
    cellranger count \\
        --id=${gex_fastq_id} \\
        --fastqs=${gex_fastqs} \\
        --transcriptome=/users/ds286q/project0001/10XReference/Human/refdata-gex-GRCh38-2020-A \\
        --include-introns=${intron} \\
        --r1-length=26 \\
        --sample=${gex_fastq_id} \\
        --localcores=20 \\
        --localmem=40 \\
        --create-bam=${bam} \\
        --output-dir=${output}/${gex_fastq_id}
fi

if [ \$? -ne 0 ]; then
    echo "cellranger count failed for sample ${gex_fastq_id}"
    exit 1
fi

echo "cellranger count completed successfully for sample ${gex_fastq_id}"
EOF

    # Submit the job
    job_id=$(sbatch "$script_file" | awk '{print $4}')
    if [ -z "$job_id" ]; then
        echo "Error: Failed to submit cellranger job for $gex_fastq_id"
        continue
    fi
    job_ids+=("$job_id")

    # Submit souporcell job if applicable
    if [ -n "$souporcell_param" ] && [ "$souporcell_param" -ne 0 ]; then
        echo "Submitting souporcell job for: $souporcell_param"
        souporcell_job_id=$(sbatch --dependency=afterok:$job_id do.souporcellV2.sh "${output}/${gex_fastq_id}" "$souporcell_param" | awk '{print $4}')
        if [ -z "$souporcell_job_id" ]; then
            echo "Error: Failed to submit souporcell job for $gex_fastq_id"
        else
            job_ids+=("$souporcell_job_id")
        fi
    fi

    # Collect unique output folders
    unique_folders["$output"]=1
done < "$input_file"

############# SUBMIT SEURAT JOBS #############
dependency_list=$(IFS=:; echo "${job_ids[*]}")
for folder in "${!unique_folders[@]}"; do
    sbatch --dependency=afterok:$dependency_list seurat_run_only.sh "$folder" || { echo "Failed to submit Seurat job for folder: $folder"; continue; }
done

echo "All jobs submitted successfully."
