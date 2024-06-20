#!/bin/bash -l

############# SLURM SETTINGS #############
#SBATCH --account=project0001
#SBATCH --job-name=cellranger_job_#ID
#SBATCH --output=cellranger_job_#ID-%j.out
#SBATCH --error=cellranger_job_#ID-%j.err
#SBATCH --partition=nodes
#SBATCH --time=1-00:00:00
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --ntasks-per-node=1

############# LOADING MODULES (optional) #############
module load apps/cellranger
module load apps/R/4.3.0

############# MY CODE #############

# cellranger count command
cellranger count --localcores=20 --localmem=40 --id=#ID --fastqs=#FASTQS --transcriptome=/users/ds286q/project0001/10XReference/Human/refdata-gex-GRCh38-2020-A --sample=#ID --r1-length=26 --include-introns=#INTRON --output-dir=/users/ds286q/project0001/Dom/pipeline/mapped/#ID

# Optionally, run additional commands or scripts here
