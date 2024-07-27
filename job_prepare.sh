#!/bin/bash
#SBATCH --job-name=prepare                  # Job name
#SBATCH --output=%x_id_%j.out                # Output file based on job name and job ID
#SBATCH --error=%x_id_%j.err                 # Error file based on job name and job ID
#SBATCH --time=48:00:00                   # Time limit
#SBATCH --partition=debug                 # Partition name
#SBATCH --nodes=1                         # Number of nodes
#SBATCH --ntasks=1                       # Number of tasks
#SBATCH --cpus-per-task=4                 # Number of CPU cores per task
#SBATCH --mem=4G                          # Memory per CPU core

# Automatically find the conda base path
# Initialize conda in the shell
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate ihb

echo "Running job with name ${SLURM_JOB_NAME} and ID ${SLURM_JOB_ID}"
start_time=$(date +%s)
echo "Job started at $(date)"

# --- Your code here ---
./01_run_prepare.sh
# ----------------------

end_time=$(date +%s)
echo "Job ended at $(date)"
# Calculate and log the total time used in hours
total_time=$((end_time - start_time))
total_hours=$(echo "scale=2; $total_time / 3600" | bc)
echo "Total time used: $total_hours hours"
