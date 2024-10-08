#!/bin/bash
#SBATCH --job-name=125w                  # Job name
#SBATCH --output=%x_id_%j.out                # Output file based on job name and job ID
#SBATCH --error=%x_id_%j.err                 # Error file based on job name and job ID
#SBATCH --partition=debug                 # Partition name
#SBATCH --nodes=1                         # Number of nodes
#SBATCH --ntasks=4                        # Number of tasks (equal to logical cores)
#SBATCH --cpus-per-task=1                 # Number of CPU cores per task
#SBATCH --mem=4G                          # Total memory for the job

echo "Running job with name ${SLURM_JOB_NAME} and ID ${SLURM_JOB_ID}"
start_time=$(date +%s)
echo "Job started at $(date)"

# --- Your code here ---
mpirun -np ${SLURM_NTASKS} lmp -in lammps.in
# --- End of your code ---

end_time=$(date +%s)
echo "Job ended at $(date)"
# Calculate and log the total time used in hours
total_time=$((end_time - start_time))
total_hours=$(echo "scale=2; $total_time / 3600" | bc)
echo "Total time used: $total_hours hours"
