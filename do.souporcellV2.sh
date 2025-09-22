#!/bin/bash -l

############# SLURM SETTINGS #############
#SBATCH --account=project0001  # account name (mandatory), if the job runs under a project then it'll be the project name, if not then it should =none
#SBATCH --job-name=souporcell        # some descriptive job name of your choice
#SBATCH --output=%x-%j.out      # output file name will contain job name + job ID
#SBATCH --error=%x-%j.err       # error file name will contain job name + job ID
#SBATCH --partition=nodes        # which partition to use, default on MARS is “nodes"
#SBATCH --time=1-01:00:00       # time limit for the whole run, in the form of d-hh:mm:ss, also accepts mm, mm:ss, hh:mm:ss, d-hh, d-hh:mm
#SBATCH --mem=50G                # memory required per node, in the form of [num][M|G|T]
#SBATCH --nodes=1               # number of nodes to allocate, default is 1
#SBATCH --ntasks=1              # number of Slurm tasks to be launched, increase for multi-process runs ex. MPI
#SBATCH --cpus-per-task=10       # number of processor cores to be assigned for each task, default is 1, increase for multi-threaded runs

############# LOADING MODULES (optional) #############
module load apps/apptainer

############# MY CODE #############
# RUN:
# sbatch thisfile.sh path n_samples
#
# $1 should be folder path containing "outs" folder of cellranger. Possorted_genome_bam.bam is needed inside the outs folder
# $2 should be n_samples (e.g. 2)

#echo $1
#echo $2

#extract barcode
#gzip <= 1.5
gzip -dc $1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > $1/outs/filtered_feature_bc_matrix/barcodes.tsv

#when they will update to gzip > 1.6
#gzip -dk $1/filtered_feature_bc_matrix/barcodes.tsv.gz

apptainer exec --bind /mnt/autofs/data:/mnt/autofs/data /users/ds286q/souporcell/souporcell_latest.sif souporcell_pipeline.py -i $1/outs/possorted_genome_bam.bam -b $1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz -f /users/ds286q/souporcell/10X_human_genome.fa -t 10 -o $1/souporcell_results -k $2
