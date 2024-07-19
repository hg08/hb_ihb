#!/bin/bash

#SBATCH --job-name=case2                  # Job name
#SBATCH --output=%x_id_%j.out                # Output file based on job name and job ID
#SBATCH --error=%x_id_%j.err                 # Error file based on job name and job ID
#SBATCH --time=INFINITE                   # Time limit
#SBATCH --partition=debug                 # Partition name
#SBATCH --ntasks=1                        # Number of tasks
#SBATCH --cpus-per-task=1                 # Number of CPU cores per task
#SBATCH --mem-per-cpu=4G                  # Memory per CPU core

echo "Running job with name ${SLURM_JOB_NAME} and ID ${SLURM_JOB_ID}"
start_time=$(date +%s)
echo "Job started at $(date)"

# --- Your code here ---
./01_run_case2.sh
# ----------------------

end_time=$(date +%s)
echo "Job ended at $(date)"
# Calculate and log the total time used in hours
total_time=$((end_time - start_time))
total_hours=$(echo "scale=2; $total_time / 3600" | bc)
echo "Total time used: $total_hours hours"
