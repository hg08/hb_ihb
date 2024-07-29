#!/bin/bash

#SBATCH --job-name=orientation                  # Job name
#SBATCH --output=%x_id_%j.out                # Output file based on job name and job ID
#SBATCH --error=%x_id_%j.err                 # Error file based on job name and job ID
#SBATCH --time=05-00:00:00                   # Time limit
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1             # Number of tasks per node
#SBATCH --mem=20G                          # Memory per CPU core

# Automatically find the conda base path
# Initialize conda in the shell
#CONDA_BASE=$(conda info --base)
#source "$CONDA_BASE/etc/profile.d/conda.sh"
#conda activate ihb
module load mamba gcc
source activate ihb

echo "Running job with name ${SLURM_JOB_NAME} and ID ${SLURM_JOB_ID}"
start_time=$(date +%s)
echo "Job started at $(date)"

# --- Your code here ---
./01_run_orientation.sh
# ----------------------

end_time=$(date +%s)
echo "Job ended at $(date)"
# Calculate and log the total time used in hours
total_time=$((end_time - start_time))
total_hours=$(echo "scale=2; $total_time / 3600" | bc)
echo "Total time used: $total_hours hours"
