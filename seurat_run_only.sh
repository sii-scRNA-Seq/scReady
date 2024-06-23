#!/bin/bash -l

############# SLURM SETTINGS #############
#SBATCH --account=project0001   # account name (mandatory), if the job runs under a project then it'll be the project name, if not then it should =none
#SBATCH --job-name=R_pipeline        # some descriptive job name of your choice
#SBATCH --output=%x-%j.out      # output file name will contain job name + job ID
#SBATCH --error=%x-%j.err       # error file name will contain job name + job ID
#SBATCH --partition=nodes        # which partition to use, default on MARS is “nodes"
#SBATCH --time=1-00:00:00       # time limit for the whole run, in the form of d-hh:mm:ss, also accepts mm, mm:ss, hh:mm:ss, d-hh, d-hh:mm
#SBATCH --mem=50G
#SBATCH --nodes=1               # number of nodes to allocate, default is 1
#SBATCH --ntasks=1              # number of Slurm tasks to be launched, increase for multi-process runs ex. MPI
#SBATCH --cpus-per-task=40       # number of processor cores to be assigned for each task, default is 1, increase for multi-threaded runs
#SBATCH --ntasks-per-node=1     # number of tasks to be launched on each allocated node

############# LOADING MODULES (optional) #############
module load apps/cellranger
module load apps/R/4.3.0

############# MY CODE #############

# Assuming input file is passed as a parameter to the script
optional_csv_file="$1"

if [ -z "$optional_csv_file" ]; then
    sbatch --account=project0001 --wrap="Rscript /users/ds286q/project0001/Dom/pipeline/SeuratGeneration.R"
else
    sbatch --account=project0001 --wrap="Rscript /users/ds286q/project0001/Dom/pipeline/SeuratGeneration.R $optional_csv_file"
fi
